"""Paso 2: persistencia, correlación C(n,tau), longitud xi(tau) y chi4(tau).

Uso:  python step2_persistence.py [--h h_xt_sub.npy --tag sub_] [--fps 25] [--um-per-px 0.17] [--eps 1 2]
Entrada: h_xt.npy (del paso 1).
Salidas: fig_Cn_tau_eps*.png, fig_xi_chi4_eps*.png, persistence_eps*.npz, resumen en stdout.

Definiciones (promedios sobre columnas i y tiempos de inicio t):
  p_i(t,tau) = 1 si |h(x_i,t+tau) - h(x_i,t)| < eps
  Pi(t,tau)  = (1/L) sum_i p_i ;  Pi(tau) = <Pi(t,tau)>_t
  C(n,tau)   = <p_i p_{i+n}> - Pi(tau)^2
  chi4(tau)  = L [ <Pi(t,tau)^2>_t - Pi(tau)^2 ]
Identidad exacta (bordes abiertos): chi4 = sum_{|n|<L} (1 - |n|/L) C(n).
Como C(n) incluye la varianza de Pi entre distintos t, tiende a una meseta
B = chi4_glob/L para n grande; se ajusta C(n) = A exp(-n/xi) + B.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit

ap = argparse.ArgumentParser()
ap.add_argument("--h", default="h_xt.npy")
ap.add_argument("--tag", default="", help="sufijo para los archivos de salida")
ap.add_argument("--fps", type=float, default=25.0)
ap.add_argument("--um-per-px", type=float, default=0.17,
                help="escala espacial en µm/px (0 = reportar solo en px)")
ap.add_argument("--eps", type=float, nargs="+", default=None,
                help="umbrales en px (default: 1, 2 y 3*sigma_ruido)")
ap.add_argument("--nmax", type=int, default=300, help="n máximo para C(n)")
ap.add_argument("--min-pers", type=int, default=200,
                help="mínimo de eventos persistentes (sum p) para ajustar xi")
args = ap.parse_args()

h = np.load(args.h)
T, L = h.shape
taus = np.arange(1, T // 2 + 1)

# ruido de segmentación: rugosidad de alta frecuencia a lo largo de x
sigma_noise = (h - gaussian_filter1d(h, 3, axis=1)).std()
# ruido temporal: la pared avanza (h decrece) y no retrocede, así que los
# Delta h(tau=1) > 0 son puro ruido de columnas quietas -> sigma de Delta h
dh1 = h[1:] - h[:-1]
sigma_dh = np.sqrt((dh1[dh1 > 0] ** 2).mean())
eps_list = args.eps if args.eps else [1.0, 2.0, round(3 * sigma_dh, 2)]
print(f"L={L} columnas, T={T} frames, sigma_ruido(alta frec. en x)={sigma_noise:.2f} px, "
      f"sigma_Dh(columnas quietas, retrocesos)={sigma_dh:.2f} px")


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
        # C(n) para todos los n (vía FFT) -> identidad exacta con chi4
        f = np.fft.rfft(p, n=2 * L, axis=1)
        ac = np.fft.irfft(f * np.conj(f), axis=1)[:, :L].mean(0)   # sum_i p_i p_{i+n}, prom. en t
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
        # suma discreta del modelo: A coth(1/2xi) (parte exponencial) + B L (meseta)
        out["chi4_model"][k] = A / np.tanh(1 / (2 * xi)) + B * L
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
    tag = f"{args.tag}eps{eps:g}"
    np.savez(f"persistence_{tag}.npz", eps=eps, **r)
    valid = np.isfinite(r["xi"])
    k_chi = int(np.argmax(r["chi4"]))
    k_xi = int(np.nanargmax(r["xi"])) if valid.any() else None
    tau_star = r["tau"][k_chi]

    # --- C(n,tau) para varios tau ---
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
    axes[0].legend(fontsize=7); axes[0].set_title(rf"$\varepsilon$={eps:g} px  (puntos: datos; líneas: $Ae^{{-n/\xi}}+B$)", fontsize=9)
    axes[1].set_xlim(0, 80); axes[1].set_ylim(1e-3, 1.5)
    axes[1].set_xlabel("n [px]"); axes[1].set_ylabel(r"$(C-B)/A$")
    axes[1].set_title("Escala semilog", fontsize=9)
    fig.tight_layout(); fig.savefig(f"fig_Cn_tau_{tag}.png", dpi=150); plt.close(fig)

    # --- xi(tau) y chi4(tau), lin y log-log ---
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5))
    t = r["tau"]
    for row, logs in enumerate([False, True]):
        ax = axes[row, 0]
        ax.errorbar(t[valid], r["xi"][valid], yerr=r["xi_err"][valid], fmt="o-", ms=3)
        ax.set_ylabel(r"$\xi(\tau)$ [px]")
        ax = axes[row, 1]
        ax.plot(t, r["chi4"], "o-", ms=3, label=r"$L\,\mathrm{Var}_t(\Pi)$")
        ax.plot(t, r["chi4_sumC"], "x", ms=5, label=r"$\sum_n (1-|n|/L)\,C(n)$")
        ax.plot(t[valid], r["chi4_model"][valid], "s", mfc="none", ms=5,
                label=r"modelo: $A\coth(1/2\xi)+BL$")
        ax.axvline(tau_star, color="r", ls="--", lw=0.8)
        ax.set_ylabel(r"$\chi_4(\tau)$")
        ax = axes[row, 2]
        ax.plot(t, r["Pi"], "o-", ms=3)
        ax.set_ylabel(r"$\Pi(\tau)$")
        for ax in axes[row]:
            ax.set_xlabel(r"$\tau$ [frames]")
            if logs:
                ax.set_xscale("log"); ax.set_yscale("log")
        ax = axes[row, 0]
        if logs:   # rango corto: etiquetas en notación simple, no 3x10^0
            for a in (ax.xaxis, ax.yaxis):
                a.set_major_formatter(FormatStrFormatter("%g"))
                a.set_minor_formatter(FormatStrFormatter("%g"))
        if args.um_per_px:
            k = args.um_per_px
            sec = ax.secondary_yaxis("right", functions=(lambda y: y * k, lambda y: y / k))
            sec.set_ylabel(r"$\xi$ [µm]")
            if logs:
                sec.yaxis.set_major_formatter(FormatStrFormatter("%g"))
                sec.yaxis.set_minor_formatter(FormatStrFormatter("%g"))
    axes[0, 1].legend(fontsize=7)
    axes[0, 0].set_title(rf"$\varepsilon$={eps:g} px", fontsize=10)
    axes[0, 1].set_title(rf"$\tau^*$={tau_star} frames", fontsize=10)
    fig.tight_layout(); fig.savefig(f"fig_xi_chi4_{tag}.png", dpi=150); plt.close(fig)

    print(f"\n=== eps = {eps:g} px ===")
    print(" tau   Pi       chi4     sumC     modelo   xi[px]        A        B*L   n_pers")
    for k in range(len(t)):
        if r["npers"][k] == 0 and k > 0 and r["npers"][k - 1] == 0:
            continue
        print(f"{t[k]:4d} {r['Pi'][k]:8.4f} {r['chi4'][k]:8.3f} {r['chi4_sumC'][k]:8.3f} "
              f"{r['chi4_model'][k]:8.3f} {r['xi'][k]:6.2f}±{r['xi_err'][k]:<5.2f} "
              f"{r['A'][k]:8.4f} {r['B'][k] * L:8.3f} {int(r['npers'][k]):6d}")
    s = (f"eps={eps:g} px: tau* (max chi4) = {tau_star} frames = {tau_star / args.fps * 1e3:.0f} ms, "
         f"chi4(tau*)={r['chi4'][k_chi]:.2f}, xi(tau*)={fmt_len(r['xi'][k_chi], r['xi_err'][k_chi])}")
    if k_xi is not None:
        s += (f"; max xi en tau={t[k_xi]} frames ({t[k_xi] / args.fps * 1e3:.0f} ms): "
              f"xi={fmt_len(r['xi'][k_xi])}")
    summary.append(s)

print("\n=== RESUMEN ===")
print(f"sigma_ruido(x) = {sigma_noise:.2f} px; sigma_Dh = {sigma_dh:.2f} px; dt = {1e3 / args.fps:.0f} ms/frame")
v = -np.diff(h.mean(1)).mean()
print(f"velocidad media de la pared = {v:.2f} px/frame"
      + (f" = {v * args.um_per_px * args.fps:.1f} µm/s" if args.um_per_px else ""))
for s in summary:
    print(s)
