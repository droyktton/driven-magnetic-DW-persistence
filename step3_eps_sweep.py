"""Sensibilidad a epsilon: xi(tau) y chi4(tau) para varios eps (h subpíxel).

Uso:  python step2_persistence.py --h h_xt_sub.npy --tag sub_ --eps 0.25 0.5 0.75 1 1.5 2 3
      python step3_eps_sweep.py [--tag sub_] [--sigma-dh 0.25] [--um-per-px 0.17]
Salida: fig_eps_sweep.png
"""
import argparse
import glob
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("--tag", default="sub_")
ap.add_argument("--sigma-dh", type=float, default=0.25, help="ruido de Delta h (px)")
ap.add_argument("--um-per-px", type=float, default=0.17,
                help="escala espacial en µm/px (0 = solo px)")
args = ap.parse_args()

files = sorted(glob.glob(f"persistence_{args.tag}eps*.npz"),
               key=lambda f: float(re.search(r"eps([\d.]+)\.npz", f).group(1)))
runs = [np.load(f) for f in files]
eps = np.array([float(r["eps"]) for r in runs])

fig, axes = plt.subplots(1, 3, figsize=(14, 4))
cmap = plt.cm.viridis(np.linspace(0, 0.9, len(runs)))
for c, r in zip(cmap, runs):
    ok = np.isfinite(r["xi"])
    axes[0].errorbar(r["tau"][ok], r["xi"][ok], yerr=r["xi_err"][ok], fmt="o-", ms=3,
                     color=c, label=rf"$\varepsilon$={float(r['eps']):g}")
    axes[1].semilogy(r["tau"], np.where(r["chi4"] > 0, r["chi4"], np.nan), "o-", ms=3, color=c)
axes[0].set_xlabel(r"$\tau$ [frames]"); axes[0].set_ylabel(r"$\xi(\tau)$ [px]")
axes[0].legend(fontsize=7)
axes[1].set_xlim(0.5, 14); axes[1].set_xlabel(r"$\tau$ [frames]"); axes[1].set_ylabel(r"$\chi_4(\tau)$")

xi1 = np.array([r["xi"][0] for r in runs])
xe1 = np.array([r["xi_err"][0] for r in runs])
chi1 = np.array([r["chi4"][0] for r in runs])
ax = axes[2]
k = args.um_per_px if args.um_per_px else 1.0
unit = "µm" if args.um_per_px else "px"
ax.errorbar(eps, xi1 * k, yerr=xe1 * k, fmt="o-", color="C0")
ax.set_xlabel(r"$\varepsilon$ [px]"); ax.set_ylabel(rf"$\xi(\tau=1)$ [{unit}]", color="C0")
ax.axvspan(0, 3 * args.sigma_dh, color="0.9", zorder=0)
ax.axvline(3 * args.sigma_dh, color="k", ls="--", lw=0.8)
ax.text(3 * args.sigma_dh, ax.get_ylim()[0], r" $3\sigma_{\Delta h}$", va="bottom", fontsize=8)
ax2 = ax.twinx()
ax2.plot(eps, chi1, "s--", color="C3", mfc="none")
ax2.set_ylabel(r"$\chi_4(\tau=1)$", color="C3")
ax.set_title(r"Sensibilidad a $\varepsilon$ ($\tau$=1 frame)", fontsize=10)
if args.um_per_px:
    axes[0].secondary_yaxis("right", functions=(lambda y: y * k, lambda y: y / k)
                            ).set_ylabel(r"$\xi$ [µm]")
fig.tight_layout()
fig.savefig("fig_eps_sweep.png", dpi=150)
for e, x, xe, ch in zip(eps, xi1, xe1, chi1):
    um = f" = {x * args.um_per_px:.2f}±{xe * args.um_per_px:.2f} µm" if args.um_per_px else ""
    print(f"eps={e:5.2f}  xi(1)={x:6.2f}±{xe:.2f} px{um}  chi4(1)={ch:6.2f}")
