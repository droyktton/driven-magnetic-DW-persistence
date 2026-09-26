"""Step 2: persistence, correlation C(n,tau), correlation length xi(tau) and chi4(tau).

Usage: python step2_persistence.py DIR [--int] [--eps 0.5 0.75 1] [--fps F] [--um-per-px X]
Input: DIR/h_xt_sub.npy (or DIR/h_xt.npy with --int) and DIR/meta.json, written by step1.
  fps and µm/px are taken from meta.json unless given explicitly.
Outputs in DIR: persistence_<tag>eps*.npz, fig_Cn_tau_<tag>eps*.png, fig_xi_chi4_<tag>eps*.png
  (tag = "sub_" for subpixel h, "" for integer h) and a summary on stdout.

Definitions (averages over columns i and start times t):
  p_i(t,tau) = 1 if |h(x_i,t+tau) - h(x_i,t)| < eps
  Pi(t,tau)  = (1/L) sum_i p_i ;  Pi(tau) = <Pi(t,tau)>_t
  C(n,tau)   = <p_i p_{i+n}> - Pi(tau)^2
  chi4(tau)  = L [ <Pi(t,tau)^2>_t - Pi(tau)^2 ]
Exact identity (open boundaries): chi4 = sum_{|n|<L} (1 - |n|/L) C(n).
Since C(n) includes the variance of Pi across start times t, it tends to a plateau
B = chi4_glob/L at large n; the fit is C(n) = A exp(-n/xi) + B.
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

ap = argparse.ArgumentParser()
ap.add_argument("dir", help="step1 output folder of one movie")
ap.add_argument("--int", action="store_true", help="use integer h (h_xt.npy) instead of subpixel")
ap.add_argument("--fps", type=float, default=None, help="default: the value in meta.json")
ap.add_argument("--um-per-px", type=float, default=None,
                help="scale in µm/px (default: the value in meta.json; 0 = report in px only)")
ap.add_argument("--eps", type=float, nargs="+", default=None,
                help="thresholds in px (default: 1, 2 and 3*sigma_Dh; with --int: 1 and 2)")
ap.add_argument("--nmax", type=int, default=300, help="maximum n for C(n)")
ap.add_argument("--min-pers", type=int, default=200,
                help="minimum number of persistent events (sum p) to fit xi")
args = ap.parse_args()

with open(os.path.join(args.dir, "meta.json")) as f:
    meta = json.load(f)
if args.fps is None:
    args.fps = meta["fps"]
if args.um_per_px is None:
    args.um_per_px = meta.get("um_per_px") or 0.0
for w in meta.get("warnings", []):
    print("WARNING (step1):", w)
tag0 = "" if args.int else "sub_"
h = np.load(os.path.join(args.dir, "h_xt.npy" if args.int else "h_xt_sub.npy"))
T, L = h.shape
taus = np.arange(1, T // 2 + 1)

# segmentation noise: high-frequency roughness along x
sigma_noise = (h - gaussian_filter1d(h, 3, axis=1)).std()
# temporal noise: the wall advances (h decreases) and does not retreat, so
# Delta h(tau=1) > 0 is pure noise from still columns -> sigma of Delta h
dh1 = h[1:] - h[:-1]
back = dh1[dh1 > 0]
sigma_dh = np.sqrt((back ** 2).mean()) if back.size else np.nan
if back.size < 500:
    print(f"WARNING: only {back.size} backward Delta h values; sigma_Dh is unreliable "
          f"(does the wall never stay still?). Pass --eps explicitly.")
if args.eps:
    eps_list = args.eps
elif args.int:   # with integer h, sigma_Dh cannot resolve the noise (< 1 px)
    eps_list = [1.0, 2.0]
else:
    eps_list = [1.0, 2.0, round(3 * sigma_dh, 2)]
eps_list = sorted(set(eps_list))
print(f"L={L} columns, T={T} frames, sigma_noise(high freq. along x)={sigma_noise:.2f} px, "
      f"sigma_Dh(still columns, backward steps)={sigma_dh:.2f} px")


def expo(n, A, xi, B):
    return A * np.exp(-n / xi) + B


def analyse(eps):
    nmax = min(args.nmax, L - 1)
    ns = np.arange(nmax + 1)
    out = dict(tau=taus, Pi=np.zeros(len(taus)), chi4=np.zeros(len(taus)),
               chi4_sumC=np.zeros(len(taus)), chi4_model=np.full(len(taus), np.nan),
               xi=np.full(len(taus), np.nan), xi_err=np.full(len(taus), np.nan),
               A=np.full(len(taus), np.nan), B=np.full(len(taus), np.nan),
               C=np.zeros((len(taus), nmax + 1)), npers=np.zeros(len(taus)))
    for k, tau in enumerate(taus):
        p = (np.abs(h[tau:] - h[:-tau]) < eps).astype(float)   # (T-tau, L)
        Pit = p.mean(1)
        Pi = Pit.mean()
        out["Pi"][k] = Pi
        out["npers"][k] = p.sum()
        out["chi4"][k] = L * Pit.var()
        # C(n) for every n (via FFT) -> exact identity with chi4
        f = np.fft.rfft(p, n=2 * L, axis=1)
        ac = np.fft.irfft(f * np.conj(f), axis=1)[:, :L].mean(0)   # sum_i p_i p_{i+n}, averaged over t
        Cfull = ac / (L - np.arange(L)) - Pi ** 2
        w = 1 - np.arange(L) / L
        out["chi4_sumC"][k] = Cfull[0] + 2 * np.sum(w[1:] * Cfull[1:])
        out["C"][k] = Cfull[: nmax + 1]
        if out["npers"][k] < args.min_pers or Pi <= 0:
            continue
        C = Cfull[: nmax + 1]
        try:
            B0 = np.median(C[nmax // 2:])
            popt, pcov = curve_fit(expo, ns, C, p0=[C[0] - B0, 5.0, B0],
                                   bounds=([0, 0.05, -np.inf], [np.inf, L, np.inf]),
                                   maxfev=20000)
        except RuntimeError:
            continue
        A, xi, B = popt
        out["A"][k], out["xi"][k], out["B"][k] = A, xi, B
        out["xi_err"][k] = np.sqrt(pcov[1, 1])
        # discrete sum of the model: A coth(1/2xi) (exponential part) + B L (plateau)
        out["chi4_model"][k] = A / np.tanh(1 / (2 * xi)) + B * L
    # chi4 normalized by the single-site variance Pi(1-Pi): an effective correlated
    # length in px (~ 2 xi when the plateau B L is negligible)
    var1 = out["Pi"] * (1 - out["Pi"])
    with np.errstate(invalid="ignore", divide="ignore"):
        out["chi4_norm"] = np.where(var1 > 0, out["chi4"] / var1, np.nan)
        out["chi4_norm_local"] = np.where(var1 > 0, (out["chi4"] - out["B"] * L) / var1, np.nan)
    return out


def fmt_len(px, err=None):
    e = lambda v: f" ± {v:.2f}" if err is not None else ""
    s = f"{px:.2f}{e(err)} px"
    if args.um_per_px:
        s += f" = {px * args.um_per_px:.2f}" + (e(err * args.um_per_px) if err is not None else "") + " µm"
    return s


summary = []
for eps in eps_list:
    r = analyse(eps)
    tag = f"{tag0}eps{eps:g}"
    np.savez(os.path.join(args.dir, f"persistence_{tag}.npz"), eps=eps, sigma_dh=sigma_dh,
             sigma_noise=sigma_noise, fps=args.fps, um_per_px=args.um_per_px, **r)
    valid = np.isfinite(r["xi"])
    k_chi = int(np.argmax(r["chi4"]))
    k_xi = int(np.nanargmax(r["xi"])) if valid.any() else None
    tau_star = r["tau"][k_chi]

    # --- C(n,tau) for several tau ---
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    ns = np.arange(r["C"].shape[1])
    ks = [k for k in np.where(valid)[0]][:8]
    cmap = plt.cm.viridis(np.linspace(0, 0.9, max(len(ks), 1)))
    for c, k in zip(cmap, ks):
        C = r["C"][k]
        lab = rf"$\tau$={r['tau'][k]}  $\xi$={r['xi'][k]:.1f}"
        axes[0].plot(ns, C, "o", ms=2, color=c, label=lab)
        axes[0].plot(ns, expo(ns, r["A"][k], r["xi"][k], r["B"][k]), "-", lw=1, color=c)
        Cc = (C - r["B"][k]) / r["A"][k]
        axes[1].semilogy(ns, np.where(Cc > 0, Cc, np.nan), "o", ms=2, color=c)
        axes[1].semilogy(ns, np.exp(-ns / r["xi"][k]), "-", lw=1, color=c)
    axes[0].set_xlim(0, 80); axes[0].set_xlabel("n [px]"); axes[0].set_ylabel(r"$C(n,\tau)$")
    axes[0].legend(fontsize=7); axes[0].set_title(rf"$\varepsilon$={eps:g} px  (dots: data; lines: $Ae^{{-n/\xi}}+B$)", fontsize=9)
    axes[1].set_xlim(0, 80); axes[1].set_ylim(1e-3, 1.5)
    axes[1].set_xlabel("n [px]"); axes[1].set_ylabel(r"$(C-B)/A$")
    axes[1].set_title("Semilog scale", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(args.dir, f"fig_Cn_tau_{tag}.png"), dpi=150); plt.close(fig)

    # --- xi(tau) and chi4(tau), linear and log-log ---
    fig, axes = plt.subplots(2, 4, figsize=(18.5, 7.5))
    t = r["tau"]
    for row, logs in enumerate([False, True]):
        ax = axes[row, 0]
        ax.errorbar(t[valid], r["xi"][valid], yerr=r["xi_err"][valid], fmt="o-", ms=3)
        ax.set_ylabel(r"$\xi(\tau)$ [px]")
        ax = axes[row, 1]
        ax.plot(t, r["chi4"], "o-", ms=3, label=r"$L\,\mathrm{Var}_t(\Pi)$")
        ax.plot(t, r["chi4_sumC"], "x", ms=5, label=r"$\sum_n (1-|n|/L)\,C(n)$")
        ax.plot(t[valid], r["chi4_model"][valid], "s", mfc="none", ms=5,
                label=r"model: $A\coth(1/2\xi)+BL$")
        ax.axvline(tau_star, color="r", ls="--", lw=0.8)
        ax.set_ylabel(r"$\chi_4(\tau)$")
        ax = axes[row, 2]
        ax.plot(t, r["Pi"], "o-", ms=3)
        ax.set_ylabel(r"$\Pi(\tau)$")
        ax = axes[row, 3]
        ax.plot(t, r["chi4_norm"], "o-", ms=3, label=r"$\chi_4/[\Pi(1-\Pi)]$")
        ax.plot(t[valid], r["chi4_norm_local"][valid], "s", mfc="none", ms=5,
                label=r"$(\chi_4 - BL)/[\Pi(1-\Pi)]$")
        ax.plot(t[valid], 2 * r["xi"][valid], "--", color="0.4", lw=1, label=r"$2\xi$")
        ax.set_ylabel(r"$\chi_4/[\Pi(1-\Pi)]$ [px]")
        for ax in axes[row]:
            ax.set_xlabel(r"$\tau$ [frames]")
            if logs:
                ax.set_xscale("log"); ax.set_yscale("log")
        for ax in (axes[row, 0], axes[row, 3]):
            if logs:   # short range: plain tick labels instead of 3x10^0
                for a in (ax.xaxis, ax.yaxis):
                    a.set_major_formatter(FormatStrFormatter("%g"))
                    a.set_minor_formatter(FormatStrFormatter("%g"))
                ax.set_xlim(left=0.95)   # tau >= 1; avoids a 0.9 tick next to 1
        ax = axes[row, 0]
        if args.um_per_px:
            k = args.um_per_px
            sec = ax.secondary_yaxis("right", functions=(lambda y: y * k, lambda y: y / k))
            sec.set_ylabel(r"$\xi$ [µm]")
            if logs:
                sec.yaxis.set_major_formatter(FormatStrFormatter("%g"))
                sec.yaxis.set_minor_formatter(FormatStrFormatter("%g"))
    axes[0, 1].legend(fontsize=7)
    axes[0, 3].legend(fontsize=7)
    axes[0, 3].set_title("Normalized $\\chi_4$", fontsize=10)
    axes[0, 0].set_title(rf"$\varepsilon$={eps:g} px", fontsize=10)
    axes[0, 1].set_title(rf"$\tau^*$={tau_star} frames", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(args.dir, f"fig_xi_chi4_{tag}.png"), dpi=150); plt.close(fig)

    print(f"\n=== eps = {eps:g} px ===")
    print(" tau   Pi       chi4     sumC     model    xi[px]        A        B*L   chi4/[Pi(1-Pi)]  n_pers")
    for k in range(len(t)):
        if r["npers"][k] == 0 and k > 0 and r["npers"][k - 1] == 0:
            continue
        print(f"{t[k]:4d} {r['Pi'][k]:8.4f} {r['chi4'][k]:8.3f} {r['chi4_sumC'][k]:8.3f} "
              f"{r['chi4_model'][k]:8.3f} {r['xi'][k]:6.2f}±{r['xi_err'][k]:<5.2f} "
              f"{r['A'][k]:8.4f} {r['B'][k] * L:8.3f} {r['chi4_norm'][k]:12.2f}     {int(r['npers'][k]):6d}")
    s = (f"eps={eps:g} px: tau* (max chi4) = {tau_star} frames = {tau_star / args.fps * 1e3:.0f} ms, "
         f"chi4(tau*)={r['chi4'][k_chi]:.2f}, xi(tau*)={fmt_len(r['xi'][k_chi], r['xi_err'][k_chi])}")
    if k_xi is not None:
        s += (f"; max xi at tau={t[k_xi]} frames ({t[k_xi] / args.fps * 1e3:.0f} ms): "
              f"xi={fmt_len(r['xi'][k_xi])}")
    summary.append(s)

print("\n=== SUMMARY ===")
print(f"sigma_noise(x) = {sigma_noise:.2f} px; sigma_Dh = {sigma_dh:.2f} px; dt = {1e3 / args.fps:.0f} ms/frame")
v = -np.diff(h.mean(1)).mean()
print(f"mean wall velocity = {v:.2f} px/frame"
      + (f" = {v * args.um_per_px * args.fps:.1f} µm/s" if args.um_per_px else ""))
for s in summary:
    print(s)
