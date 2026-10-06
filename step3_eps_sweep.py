"""Sensitivity to epsilon: xi(tau) and chi4(tau) for several eps.

Usage: python step2_persistence.py DIR --eps 0.25 0.5 0.75 1 1.5 2 3
      python step3_eps_sweep.py DIR [--int] [--row-mean] [--sigma-dh S] [--um-per-px X]
sigma_Dh and µm/px are read from the step2 .npz files unless given explicitly.
Output: DIR/fig_eps_sweep_sub.png (_int.png with --int; _sub_rm.png / _rm.png with --row-mean)
"""
import argparse
import glob
import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("dir", help="folder of one movie (output of step1/step2)")
ap.add_argument("--int", action="store_true", help="use the integer-h results")
ap.add_argument("--row-mean", action="store_true", help="use the step2 --row-mean results")
ap.add_argument("--sigma-dh", type=float, default=None, help="noise of Delta h (px)")
ap.add_argument("--um-per-px", type=float, default=None,
                help="spatial scale in µm/px (0 = px only)")
args = ap.parse_args()
tag = ("" if args.int else "sub_") + ("rm_" if args.row_mean else "")

files = sorted(glob.glob(os.path.join(args.dir, f"persistence_{tag}eps*.npz")),
               key=lambda f: float(re.search(r"eps([\d.]+)\.npz", f).group(1)))
if len(files) < 2:
    raise SystemExit(f"at least 2 eps values are needed in {args.dir}; run step2 with --eps")
runs = [np.load(f) for f in files]
eps = np.array([float(r["eps"]) for r in runs])
if args.sigma_dh is None:
    args.sigma_dh = float(runs[0]["sigma_dh"])
if args.um_per_px is None:
    args.um_per_px = float(runs[0]["um_per_px"])

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
if not args.int and np.isfinite(args.sigma_dh):   # with integer h, sigma_Dh cannot resolve the noise
    ax.axvspan(0, 3 * args.sigma_dh, color="0.9", zorder=0)
    ax.axvline(3 * args.sigma_dh, color="k", ls="--", lw=0.8)
    ax.text(3 * args.sigma_dh, ax.get_ylim()[0], r" $3\sigma_{\Delta h}$", va="bottom", fontsize=8)
ax2 = ax.twinx()
ax2.plot(eps, chi1, "s--", color="C3", mfc="none")
ax2.set_ylabel(r"$\chi_4(\tau=1)$", color="C3")
ax.set_title(r"Sensitivity to $\varepsilon$ ($\tau$=1 frame)", fontsize=10)
if args.um_per_px:
    axes[0].secondary_yaxis("right", functions=(lambda y: y * k, lambda y: y / k)
                            ).set_ylabel(r"$\xi$ [µm]")
fig.tight_layout()
fig.savefig(os.path.join(args.dir, f"fig_eps_sweep_{tag or 'int_'}".rstrip("_") + ".png"),
            dpi=150)
for e, x, xe, ch in zip(eps, xi1, xe1, chi1):
    um = f" = {x * args.um_per_px:.2f}±{xe * args.um_per_px:.2f} µm" if args.um_per_px else ""
    print(f"eps={e:5.2f}  xi(1)={x:6.2f}±{xe:.2f} px{um}  chi4(1)={ch:6.2f}")
