"""Sensitivity to epsilon: xi(tau), chi4(tau) and their characteristic values for several eps.

Panels: xi(tau) and chi4(tau) for each eps (stars at tau*, the maximum of chi4); xi(tau*), xi_max
and xi(tau=1) vs eps, with chi4(tau*) on the right axis; tau*, the lag of xi_max and the lag at
which Pi = 1/2 vs eps.

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

fps = float(runs[0]["fps"]) if "fps" in runs[0].files else 0.0
k = args.um_per_px if args.um_per_px else 1.0
unit = "µm" if args.um_per_px else "px"


def at(r, i):
    return r["xi"][i], r["xi_err"][i]


# characteristic values per eps
ks = [int(np.argmax(r["chi4"])) for r in runs]                     # tau*: max chi4
kx = [int(np.nanargmax(r["xi"])) if np.isfinite(r["xi"]).any() else 0 for r in runs]
tau_star = np.array([r["tau"][i] for r, i in zip(runs, ks)])
chi_star = np.array([r["chi4"][i] for r, i in zip(runs, ks)])
xi_star = np.array([at(r, i) for r, i in zip(runs, ks)])          # (n, 2): value, error
xi_max = np.array([r["xi"][i] for r, i in zip(runs, kx)])
tau_xi = np.array([r["tau"][i] for r, i in zip(runs, kx)])
tau_half = np.array([r["tau"][np.argmax(r["Pi"] < 0.5)] if (r["Pi"] < 0.5).any() else np.nan for r in runs])
xi_one = np.array([at(r, 0) for r in runs])

fig, axes = plt.subplots(1, 4, figsize=(19, 4.2))
cmap = plt.cm.viridis(np.linspace(0, 0.9, len(runs)))
for c, r, i in zip(cmap, runs, ks):
    ok = np.isfinite(r["xi"])
    axes[0].errorbar(r["tau"][ok], r["xi"][ok], yerr=r["xi_err"][ok], fmt="o-", ms=3,
                     color=c, label=rf"$\varepsilon$={float(r['eps']):g}")
    if np.isfinite(r["xi"][i]):
        axes[0].plot(r["tau"][i], r["xi"][i], "*", ms=12, color=c, mec="k")
    axes[1].plot(r["tau"], r["chi4"], "o-", ms=3, color=c)
    axes[1].plot(r["tau"][i], r["chi4"][i], "*", ms=12, color=c, mec="k")
axes[0].set_xscale("log"); axes[0].set_xlabel(r"$\tau$ [frames]"); axes[0].set_ylabel(r"$\xi(\tau)$ [px]")
axes[0].legend(fontsize=7); axes[0].set_title(r"$\xi(\tau)$ (stars: $\tau^*$)", fontsize=10)
axes[1].set_xlim(0.5, max(14, 3 * tau_star.max())); axes[1].set_ylim(bottom=0); axes[1].set_xlabel(r"$\tau$ [frames]")
axes[1].set_ylabel(r"$\chi_4(\tau)$"); axes[1].set_title(r"$\chi_4(\tau)$ (stars: $\tau^*$)", fontsize=10)

ax = axes[2]
ax.errorbar(eps, xi_star[:, 0] * k, yerr=xi_star[:, 1] * k, fmt="o-", color="C0", label=r"$\xi(\tau^*)$")
ax.plot(eps, xi_max * k, "D-", color="C2", mfc="none", label=r"$\xi_{max}$")
ax.plot(eps, xi_one[:, 0] * k, ".:", color="0.5", label=r"$\xi(\tau=1)$")
ax.set_xlabel(r"$\varepsilon$ [px]"); ax.set_ylabel(f"length [{unit}]")
ax2 = ax.twinx()
ax2.plot(eps, chi_star, "s--", color="C3", mfc="none")
ax2.set_ylabel(r"$\chi_4(\tau^*)$", color="C3")
ax.legend(fontsize=8, loc="lower right")
ax.set_title(r"Correlation lengths vs $\varepsilon$", fontsize=10)

ax = axes[3]
ax.plot(eps, tau_star, "o-", color="C0", label=r"$\tau^*$ (max $\chi_4$)")
ax.plot(eps, tau_xi, "D-", color="C2", mfc="none", label=r"$\tau$ of $\xi_{max}$")
ax.plot(eps, tau_half, "^:", color="C1", label=r"$\tau$ at $\Pi=1/2$")
ax.set_yscale("log"); ax.set_xlabel(r"$\varepsilon$ [px]"); ax.set_ylabel(r"$\tau$ [frames]")
if fps:
    ax.secondary_yaxis("right", functions=(lambda f: f / fps, lambda s: s * fps)).set_ylabel(r"$\tau$ [s]")
ax.legend(fontsize=8); ax.set_title(r"Time scales vs $\varepsilon$", fontsize=10)

for ax in axes[2:]:
    if not args.int and np.isfinite(args.sigma_dh):   # with integer h, sigma_Dh cannot resolve the noise
        ax.axvspan(0, 3 * args.sigma_dh, color="0.9", zorder=0)
        ax.axvline(3 * args.sigma_dh, color="k", ls="--", lw=0.8)
axes[2].text(3 * args.sigma_dh, axes[2].get_ylim()[0], r" $3\sigma_{\Delta h}$", va="bottom", fontsize=8) \
    if not args.int and np.isfinite(args.sigma_dh) else None
if args.um_per_px:
    axes[0].secondary_yaxis("right", functions=(lambda y: y * k, lambda y: y / k)
                            ).set_ylabel(r"$\xi$ [µm]")
fig.tight_layout()
fig.savefig(os.path.join(args.dir, f"fig_eps_sweep_{tag or 'int_'}".rstrip("_") + ".png"),
            dpi=150)
um = lambda v: f" ({v * args.um_per_px:.2f} µm)" if args.um_per_px else ""
print(" eps   tau*  tau(Pi=1/2)  chi4(tau*)   xi(tau*) [px]        xi_max [px] (tau)")
for e, ts, th, ch, (x, xe), xm, tx in zip(eps, tau_star, tau_half, chi_star, xi_star, xi_max, tau_xi):
    print(f"{e:5.2f}  {ts:4d}  {th:6.0f}       {ch:6.2f}    {x:6.2f}±{xe:.2f}{um(x)}   {xm:6.2f}{um(xm)} ({tx})")
