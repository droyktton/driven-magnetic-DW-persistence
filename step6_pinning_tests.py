"""Step 6: strong pinning vs collective creep (TIFF studies).

Usage: python step6_pinning_tests.py DIR [--eps 0.97] [--tau T] [--nsurr 200]
Inputs: DIR/h_xt_sub.npy, DIR/meta.json, DIR/tiff_extra.npz (step1); optionally the step2 results
        DIR/persistence_sub_rm_eps<eps>.npz, DIR/persistence_sub_rm_md<D>_eps<eps>.npz (far from the
        defects) and DIR/persistence_sub_rm_nd<D>_eps<eps>.npz (near them)
        (step2 --row-mean [--mask-defects D [--near-defects]]).
Outputs: DIR/fig_mask_defects_eps<eps>.png (if masked step2 results exist), DIR/fig_pinning_tests.png, summary.

Tests:
 1. Near vs far from the static defects: xi(tau), chi4(tau), chi4/[Pi(1-Pi)], tau* and xi_max from
    step2 --mask-defects D (far) and --mask-defects D --near-defects (near), overlaid with the
    whole wall. If the chi4 peak came from the defects it would be present near them and absent far
    from them.
 2. Location of the persistent clusters at tau*: a cluster is a run of consecutive rows with
    p_i(t, tau*) = 1 at a given t. For clusters longer than 2 xi(tau*), the distance from their
    centre (on the wall) to the nearest defect is compared with a surrogate that moves each cluster
    to a random position along the wall at the same t (--nsurr replicas). Also: mean cluster
    length vs distance to the nearest defect.
 3. Growth of xi(tau): fits of a power law a tau^b, a log law a [ln(1+tau)]^b and a saturating
    form a (1 - exp(-tau/c)) + d, compared with the AIC.
 4. Shape of the slow regions (possible hidden strong pins): waiting map w = 1/v (frames per px)
    from the arrival-time map. A strong pin is a point-like spot with very long waiting that holds a
    large share of the total time; collective creep gives extended slow regions. Reported: share of
    the total waiting time in the slowest 1 % and 5 % of the pixels, and correlation lengths of
    ln w along the wall (y) and along the motion (x), to compare with xi.
"""
import argparse
import glob
import json
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import binary_dilation, binary_erosion, distance_transform_edt, gaussian_filter
from scipy.optimize import curve_fit

ap = argparse.ArgumentParser()
ap.add_argument("dir", help="step1/step2 output folder of a TIFF study")
ap.add_argument("--eps", type=float, default=0.97, help="threshold in px")
ap.add_argument("--tau", type=int, default=None, help="lag for the cluster test (default tau*)")
ap.add_argument("--nsurr", type=int, default=200, help="surrogate replicas")
ap.add_argument("--seed", type=int, default=1)
args = ap.parse_args()
rng = np.random.default_rng(args.seed)

with open(os.path.join(args.dir, "meta.json")) as f:
    meta = json.load(f)
k_um, fps, Wc = meta.get("um_per_px") or 0.0, meta["fps"], meta["height"]
h = np.load(os.path.join(args.dir, "h_xt_sub.npy"))
T, L = h.shape
e = np.load(os.path.join(args.dir, "tiff_extra.npz"))
defects, ta = e["defects"], e["arrival"]
x = Wc - 1 - h                                             # wall position x_eff (subpixel), (T, L)
dmap = distance_transform_edt(~defects)
xi_ = np.clip(np.round(x).astype(int), 0, dmap.shape[1] - 1)
dwall = dmap[np.arange(L)[None, :], xi_]                    # distance wall point -> nearest defect
um = (lambda v: f" ({v * k_um:.2f} µm)") if k_um else (lambda v: "")
summary = []

# ---------------------------------------------------------------- test 1: masked step2 results
def load(tag):
    f = os.path.join(args.dir, f"persistence_sub_rm_{tag}eps{args.eps:g}.npz")
    return np.load(f) if os.path.exists(f) else None

base = load("")
far, near = {}, {}
for kind, store in (("md", far), ("nd", near)):
    for f in glob.glob(os.path.join(args.dir, f"persistence_sub_rm_{kind}*_eps{args.eps:g}.npz")):
        store[float(re.search(rf"_{kind}([\d.]+)_eps", f).group(1))] = np.load(f)


def char(r):
    k = int(np.argmax(r["chi4"]))
    kx = int(np.nanargmax(r["xi"])) if np.isfinite(r["xi"]).any() else k
    return dict(kept=r["kept"][k] if "kept" in r.files else 1.0, tau=r["tau"][k], chi4=r["chi4"][k], xi=r["xi"][k], dxi=r["xi_err"][k],
                ximax=r["xi"][kx], tauxi=r["tau"][kx], norm=r["chi4_norm"][k])


if far or near:
    print(f"\n[1] defect exclusion, eps = {args.eps:g} px (row-mean); far: wall > D px from every defect "
          f"during [t, t+tau]; near: the complement")
    print("   subset      kept@tau*  tau*  chi4(tau*)  chi4/[Pi(1-Pi)]   xi(tau*) [px]        xi_max [px] (tau)")
    for name, store in (("far", far), ("near", near)):
        for D in sorted(store):
            c = char(store[D])
            print(f"   {name:4s} D={D:<4g}  {c['kept']:6.0%}   {c['tau']:4d}  {c['chi4']:8.2f}   {c['norm']:10.1f}      "
                  f"{c['xi']:6.2f}±{c['dxi']:.2f}{um(c['xi'])}   {c['ximax']:6.2f}{um(c['ximax'])} ({c['tauxi']})")
    if base is not None and 0 in far:
        n0 = min(len(far[0]["tau"]), len(base["tau"]))
        print(f"   check D=0 vs unmasked run (tau <= {n0}): max |Δxi| = "
              f"{np.nanmax(np.abs(far[0]['xi'][:n0] - base['xi'][:n0])):.1e}, "
              f"max |Δchi4| = {np.max(np.abs(far[0]['chi4'][:n0] - base['chi4'][:n0])):.1e}")
    Dpair = sorted(set(far) & set(near) - {0})
    fig, ax = plt.subplots(1, 4, figsize=(22, 4.6))
    curves = [("all", base, "k", "-")] if base is not None else []
    for j, D in enumerate(Dpair):
        curves += [(f"far, D={D:g}", far[D], plt.cm.Blues(0.6 + 0.3 * j / max(len(Dpair) - 1, 1)), "-"),
                   (f"near, D={D:g}", near[D], plt.cm.Reds(0.6 + 0.3 * j / max(len(Dpair) - 1, 1)), "--")]
    taus_star = []
    for lab_, r, c, ls in curves:
        ok = np.isfinite(r["xi"])
        k = int(np.argmax(r["chi4"])); taus_star.append(r["tau"][k])
        ax[0].plot(r["tau"][ok], r["xi"][ok], ls, color=c, label=lab_)
        ax[1].plot(r["tau"], r["chi4"], ls, color=c, label=lab_)
        ax[1].plot(r["tau"][k], r["chi4"][k], "*", ms=10, color=c, mec="k")
        ax[2].plot(r["tau"], r["chi4_norm"], ls, color=c, label=lab_)
    ax[0].set_xscale("log"); ax[0].set_xlabel(r"$\tau$ [frames]"); ax[0].set_ylabel(r"$\xi(\tau)$ [px]")
    ax[0].legend(fontsize=8); ax[0].set_title(rf"$\xi(\tau)$ near / far from defects ($\varepsilon$={args.eps:g} px)", fontsize=10)
    for a_ in ax[1:3]:
        a_.set_xlim(0.5, 5 * max(taus_star)); a_.set_ylim(bottom=0); a_.set_xlabel(r"$\tau$ [frames]")
    ax[1].set_ylabel(r"$\chi_4(\tau)$ (per $L_{eff}$)"); ax[1].set_title(r"$\chi_4(\tau)$ (stars: $\tau^*$)", fontsize=10)
    ax[2].set_ylabel(r"$\chi_4/[\Pi(1-\Pi)]$ [px]"); ax[2].set_title(r"Normalized $\chi_4$", fontsize=10)
    if far:
        # D = 0 is the whole wall: take it from the full unmasked run (the D = 0 check may use fewer lags)
        src = {D: far[D] for D in far if D > 0}
        if base is not None:
            src[0.0] = base
        Ds = sorted(src)
        rows = np.array([[D, char(src[D])["kept"] if "kept" in src[D].files else 1.0, char(src[D])["tau"],
                          char(src[D])["xi"], char(src[D])["dxi"], char(src[D])["ximax"]] for D in Ds])
        ax[3].errorbar(rows[:, 0], rows[:, 3], yerr=rows[:, 4], fmt="o-", label=r"$\xi(\tau^*)$, far")
        ax[3].plot(rows[:, 0], rows[:, 5], "D-", mfc="none", label=r"$\xi_{max}$, far")
        if near:
            nr = np.array([[D, char(near[D])["xi"], char(near[D])["ximax"]] for D in sorted(near)])
            ax[3].plot(nr[:, 0], nr[:, 1], "o--", color="C3", mfc="none", label=r"$\xi(\tau^*)$, near")
            ax[3].plot(nr[:, 0], nr[:, 2], "D--", color="C3", mfc="none", alpha=0.6, label=r"$\xi_{max}$, near")
        ax2 = ax[3].twinx(); ax2.plot(rows[:, 0], rows[:, 2], "s:", color="0.4", mfc="none")
        ax2.set_ylabel(r"$\tau^*$ far [frames]", color="0.4"); ax2.set_ylim(0, 2 * rows[:, 2].max())
        ax[3].set_xlabel("distance D [px]"); ax[3].set_ylabel("length [px]"); ax[3].set_ylim(bottom=0)
        ax[3].legend(fontsize=7, loc="lower left"); ax[3].set_title("vs D (labels: share of data kept, far)", fontsize=10)
        for D, kp in zip(rows[:, 0], rows[:, 1]):
            ax[3].annotate(f"{kp:.0%}", (D, rows[:, 5].max() * 1.04), ha="center", fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(args.dir, f"fig_mask_defects_eps{args.eps:g}.png"), dpi=120); plt.close(fig)

# ---------------------------------------------------------------- test 2: cluster locations
if args.tau is None:
    args.tau = int(np.argmax(base["chi4"])) + 1 if base is not None else 7
tau = args.tau
xi_star = float(base["xi"][tau - 1]) if base is not None else 13.0
p = np.abs(h[tau:] - h[:-tau]) < args.eps                  # (T - tau, L)
# runs of consecutive persistent rows at each t
pad = np.zeros((p.shape[0], 1), bool)
dp = np.diff(np.hstack([pad, p, pad]).astype(np.int8), axis=1)
tt_s, ii_s = np.nonzero(dp == 1)
tt_e, ii_e = np.nonzero(dp == -1)                           # same order (row-major)
ell = ii_e - ii_s
ic = (ii_s + ii_e - 1) // 2
dcl = dwall[tt_s, ic]
big = ell > 2 * xi_star
d_th = (10, 20, 40)
obs = np.array([np.mean(dcl[big] <= d) for d in d_th])
sur = np.empty((args.nsurr, len(d_th)))
tb, lb = tt_s[big], ell[big]
for r_ in range(args.nsurr):
    c = (rng.random(tb.size) * (L - lb) + lb / 2).astype(int)   # random centre, cluster fits in [0, L)
    ds = dwall[tb, c]
    sur[r_] = [np.mean(ds <= d) for d in d_th]
pval = [(np.sum(sur[:, j] >= obs[j]) + 1) / (args.nsurr + 1) for j in range(len(d_th))]
print(f"\n[2] persistent clusters at tau = {tau} frames, eps = {args.eps:g} px: {ell.size} clusters, "
      f"{big.sum()} longer than 2 xi(tau*) = {2 * xi_star:.1f} px")
for j, d in enumerate(d_th):
    print(f"   centre within {d:3d} px of a defect: observed {obs[j]:.3f}, surrogate {sur[:, j].mean():.3f} "
          f"± {sur[:, j].std():.3f}, ratio {obs[j] / sur[:, j].mean():.2f}, p = {pval[j]:.3f}")
bins = np.array([0, 10, 20, 40, 80, 160, np.inf])
lmean = [ell[(dcl >= bins[i]) & (dcl < bins[i + 1])].mean() for i in range(len(bins) - 1)]
nb = [np.sum((dcl >= bins[i]) & (dcl < bins[i + 1])) for i in range(len(bins) - 1)]
print("   mean cluster length by distance to the nearest defect: "
      + ", ".join(f"[{bins[i]:g},{bins[i + 1]:g}) {lmean[i]:.1f} px (n={nb[i]})" for i in range(len(bins) - 1)))
# control of the method: clusters planted at the wall point closest to a defect (one per real
# cluster, same t) must give a clear excess over the surrogate
planted = dwall[tb].min(axis=1)
print("   control, clusters planted at the defects: "
      + ", ".join(f"within {d} px {np.mean(planted <= d):.3f} (ratio {np.mean(planted <= d) / sur[:, j].mean():.1f})"
                  for j, d in enumerate(d_th)))

# ---------------------------------------------------------------- test 3: growth of xi(tau)
def aic(y, yhat, k):
    n = y.size
    return n * np.log(np.sum((y - yhat) ** 2) / n) + 2 * k

fits = {}
if base is not None:
    t_, xi_ = base["tau"], base["xi"]
    kx = int(np.nanargmax(xi_))
    sel = np.isfinite(xi_) & (t_ <= t_[kx])                # growth part, up to the maximum
    t_, y_ = t_[sel].astype(float), xi_[sel]
    models = {
        "power  a τ^b": (lambda t, a, b: a * t ** b, [5, 0.3]),
        "log    a [ln(1+τ)]^b": (lambda t, a, b: a * np.log1p(t) ** b, [5, 0.5]),
        "satur. a(1-e^{-τ/c})+d": (lambda t, a, c, d: a * (1 - np.exp(-t / c)) + d, [15, 20, 5]),
    }
    print(f"\n[3] growth of xi(tau) for tau = 1..{int(t_[-1])} (eps = {args.eps:g}, row-mean):")
    for name, (fn, p0) in models.items():
        try:
            popt, _ = curve_fit(fn, t_, y_, p0=p0, maxfev=20000)
            yhat = fn(t_, *popt)
            fits[name] = (fn, popt, aic(y_, yhat, len(popt)))
            print(f"   {name:24s} params {np.round(popt, 3)}  AIC {fits[name][2]:8.1f}  rms {np.sqrt(np.mean((y_ - yhat) ** 2)):.2f} px")
        except RuntimeError:
            print(f"   {name:24s} fit failed")

# ---------------------------------------------------------------- test 4: shape of the slow regions
fin = np.isfinite(ta) & (ta >= 1)
w_ = gaussian_filter(np.where(fin, ta, 0.0), 2)
n_ = gaussian_filter(fin.astype(float), 2)
with np.errstate(divide="ignore", invalid="ignore"):
    tt = w_ / n_
    gy, gx = np.gradient(tt)
    wait = np.hypot(gx, gy)                                   # frames per px (= 1/v)
inner = binary_erosion(fin, iterations=6) & ~binary_dilation(defects, iterations=8) & np.isfinite(wait)
wv = np.sort(wait[inner])[::-1]
share1 = wv[: wv.size // 100].sum() / wv.sum()
share5 = wv[: wv.size // 20].sum() / wv.sum()
lw = np.where(inner, np.log(np.clip(wait, 1e-3, None)), np.nan)
lw = lw - np.nanmean(lw)


def acorr(a, axis, nmax=80):
    out = []
    va = np.nanmean(a * a)
    for n in range(nmax):
        if axis == 0:
            prod = a[: a.shape[0] - n] * a[n:]
        else:
            prod = a[:, : a.shape[1] - n] * a[:, n:]
        out.append(np.nanmean(prod) / va)
    return np.array(out)

acy, acx = acorr(lw, 0), acorr(lw, 1)
ly, lx = int(np.argmax(acy < 1 / np.e)), int(np.argmax(acx < 1 / np.e))
print(f"\n[4] waiting map 1/v (smoothing 2 px, defects excluded): slowest 1 % of pixels hold {share1:.1%} "
      f"of the waiting time, slowest 5 % hold {share5:.1%}; correlation length of ln(1/v): "
      f"along the wall (y) {ly} px{um(ly)}, along the motion (x) {lx} px{um(lx)}; xi(tau*) = {xi_star:.1f} px{um(xi_star)}")

# ---------------------------------------------------------------- figure
fig, ax = plt.subplots(1, 4, figsize=(21, 4.6))
dd = np.arange(0, 121, 2)
ax[0].plot(dd, [np.mean(dcl[big] <= d) for d in dd], "-", color="C3", label=rf"clusters $\ell > 2\xi$ (n={big.sum()})")
sc = np.array([[np.mean(dwall[tb, (rng.random(tb.size) * (L - lb) + lb / 2).astype(int)] <= d) for d in dd]
               for _ in range(30)])
ax[0].fill_between(dd, sc.mean(0) - 2 * sc.std(0), sc.mean(0) + 2 * sc.std(0), color="0.7", label="surrogate ±2σ")
ax[0].set_xlabel("distance to nearest defect d [px]"); ax[0].set_ylabel("fraction of cluster centres within d")
ax[0].legend(fontsize=8); ax[0].set_title(rf"Location of persistent clusters ($\tau$={tau}, $\varepsilon$={args.eps:g})", fontsize=10)
ax[1].plot(range(len(lmean)), lmean, "o-")
ax[1].set_xticks(range(len(lmean))); ax[1].set_xticklabels([f"{bins[i]:g}–{bins[i + 1]:g}" for i in range(len(lmean))], fontsize=8)
ax[1].set_xlabel("distance to nearest defect [px]"); ax[1].set_ylabel(r"mean cluster length $\ell$ [px]")
ax[1].set_title("Cluster length vs distance to defects", fontsize=10)
if base is not None:
    ax[2].plot(base["tau"], base["xi"], ".", color="k", ms=3, label="data")
    tt2 = np.linspace(1, base["tau"][int(np.nanargmax(base["xi"]))], 300)
    for c, (name, (fn, popt, A)) in zip(("C0", "C1", "C2"), fits.items()):
        ax[2].plot(tt2, fn(tt2, *popt), "-", color=c, label=f"{name.split()[0]} (AIC {A:.0f})")
    ax[2].set_xscale("log"); ax[2].set_xlabel(r"$\tau$ [frames]"); ax[2].set_ylabel(r"$\xi$ [px]")
    ax[2].legend(fontsize=8); ax[2].set_title(r"Growth of $\xi(\tau)$", fontsize=10)
ax[3].plot(acy, label=f"along the wall (y): {ly} px")
ax[3].plot(acx, label=f"along the motion (x): {lx} px")
ax[3].axhline(1 / np.e, color="k", ls="--", lw=0.6); ax[3].axvline(xi_star, color="C3", ls=":", label=rf"$\xi(\tau^*)$")
ax[3].set_xlabel("distance [px]"); ax[3].set_ylabel(r"autocorrelation of ln(1/v)"); ax[3].legend(fontsize=8)
ax[3].set_title(f"Slow regions (slowest 1 %: {share1:.0%} of waiting time)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(args.dir, "fig_pinning_tests.png"), dpi=120); plt.close(fig)
print(f"\noutputs: {args.dir}/fig_pinning_tests.png" + (f", fig_mask_defects_eps{args.eps:g}.png" if (far or near) else ""))
