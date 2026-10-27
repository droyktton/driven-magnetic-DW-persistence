"""Static walls (one snapshot per pulse) -> step7-style inputs for step8 / step10 part A.

Usage: python lps_walls.py SRC_DIR OUT_DIR [--um-per-px 0.171] [--glob '*_u.dat']
SRC_DIR: folder with the walls already extracted as two-column files x [um], u [um] on a uniform x grid
(LPS 2019 GdFeCo B1742 series, extracted by Lucas; single valued by construction).
Outputs in OUT_DIR: local_width.npz, roughness.npz, meta.json, fig_walls.png.
Differences from step7: walls are independent, so each wall is one bootstrap block (instead of 12 time blocks);
there are no folds (curve "unfolded" = "pca"). S(q) from non-overlapping pieces of the shortest wall length
(Hann window, as step7), averaged per wall.
"""
import argparse
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from wall_tools import pca_width, structure_factor, fit_power

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("out")
ap.add_argument("--um-per-px", type=float, default=0.171)
ap.add_argument("--glob", default="*_u.dat")
args = ap.parse_args()
os.makedirs(args.out, exist_ok=True)
kk = args.um_per_px

files = sorted(glob.glob(os.path.join(args.src, args.glob)))
walls = []
for f in files:
    x, u = np.loadtxt(f, unpack=True)
    dx = np.diff(x)
    if not np.allclose(dx, kk, rtol=1e-3):
        raise SystemExit(f"{f}: x step {dx.mean():.4f} um, expected {kk}")
    walls.append(u / kk)                                     # px
lens = np.array([len(u) for u in walls])
Lmax = lens.max()
ells = np.unique(np.round(np.logspace(np.log10(4), np.log10(Lmax), 28)).astype(int)).astype(float)

# local width with local rotation (PCA) and without rotation (line fit along x), sums and counts per wall
nw, nl = len(walls), len(ells)
S_pca, N_pca, S_ctl, N_ctl = (np.zeros((nw, nl)) for _ in range(4))
tilts = [[] for _ in ells]
for i, u in enumerate(walls):
    s = np.arange(len(u), dtype=float)
    w2l, tl, _ = pca_width(s, u, 0.0, ells, step=1.0)
    for j, l in enumerate(ells):
        S_pca[i, j], N_pca[i, j] = w2l[j].sum(), w2l[j].size
        tilts[j].append(tl[j])
        li = int(l)
        if li > len(u):
            continue
        xs = np.arange(li)
        P = np.vstack([xs, np.ones(li)]).T
        P = P @ np.linalg.pinv(P)
        segs = np.array([u[a:a + li] for a in range(0, len(u) - li + 1, max(li // 2, 1))])
        r = segs - segs @ P.T
        S_ctl[i, j], N_ctl[i, j] = (r ** 2).mean(1).sum(), len(segs)
w2_pca = np.where(N_pca.sum(0) > 0, S_pca.sum(0) / np.maximum(N_pca.sum(0), 1), np.nan)
w2_ctl = np.where(N_ctl.sum(0) > 0, S_ctl.sum(0) / np.maximum(N_ctl.sum(0), 1), np.nan)
tilt_sd = np.array([np.std(np.concatenate(t)) if sum(len(a) for a in t) > 1 else np.nan for t in tilts])
rmin, rmax = 16.0, Lmax / 4
z_pca = fit_power(ells, w2_pca, rmin, rmax)[0] / 2
z_ctl = fit_power(ells, w2_ctl, rmin, rmax)[0] / 2
rng = np.random.default_rng(1)
bz = []
for _ in range(300):
    p = rng.integers(0, nw, nw)
    n_ = N_pca[p].sum(0)
    bz.append(fit_power(ells, np.where(n_ > 0, S_pca[p].sum(0) / np.maximum(n_, 1), np.nan), rmin, rmax)[0] / 2)
dz_pca = float(np.nanstd(bz))
blk = {"pca": (S_pca, N_pca), "unfolded": (S_pca, N_pca), "ctrl": (S_ctl, N_ctl)}
np.savez(os.path.join(args.out, "local_width.npz"), ells=ells, w2_pca=w2_pca, w2_ctrl=w2_ctl, w2_unfolded=w2_pca,
         frac_folded=np.zeros(nl), tilt_sd=tilt_sd, n_segments=N_pca.sum(0), zeta_pca=z_pca, dzeta_pca=dz_pca,
         zeta_ctrl=z_ctl, rmin=rmin, rmax=rmax, um_per_px=kk, n_walls=nw, wall_len=lens,
         **{f"blk_sum_{k}": v[0] for k, v in blk.items()}, **{f"blk_cnt_{k}": v[1] for k, v in blk.items()})

# global S(q): pieces of the shortest wall length, line removed, one average per wall
Lc = int(lens.min())
Sblk = []
for u in walls:
    pcs = np.array([u[a:a + Lc] for a in range(0, len(u) - Lc + 1, Lc)])
    xs = np.arange(Lc)
    pcs = pcs - np.array([np.polyval(np.polyfit(xs, p, 1), xs) for p in pcs])
    q, S = structure_factor(pcs)
    Sblk.append(S)
Sblk = np.array(Sblk)
np.savez(os.path.join(args.out, "roughness.npz"), h_q=q, h_S=Sblk.mean(0), h_Sblk=Sblk, um_per_px=kk,
         piece_len=Lc, rmin=rmin, rmax=rmax)
json.dump(dict(source=os.path.abspath(args.src), um_per_px=kk, n_walls=nw, len_px=[int(lens.min()), int(Lmax)],
               static_walls=True), open(os.path.join(args.out, "meta.json"), "w"), indent=1)

fig, ax = plt.subplots(1, 3, figsize=(17, 5))
for u in walls[:: max(1, nw // 8)]:
    ax[0].plot(np.arange(len(u)) * kk, (u - u.mean()) * kk, lw=0.8)
ax[0].set_xlabel("x [µm]"); ax[0].set_ylabel("u − ⟨u⟩ [µm]"); ax[0].set_title(f"{nw} walls, {lens.min() * kk:.0f}–{Lmax * kk:.0f} µm")
ax[1].loglog(ells * kk, w2_pca * kk ** 2, "o-", ms=3, label=f"local rotation (PCA): ζ = {z_pca:.2f} ± {dz_pca:.2f}")
ax[1].loglog(ells * kk, w2_ctl * kk ** 2, "s--", ms=3, mfc="none", label=f"line removed, no rotation: ζ = {z_ctl:.2f}")
ax[1].axvspan(rmin * kk, rmax * kk, color="0.9"); ax[1].set_xlabel("ℓ [µm]"); ax[1].set_ylabel("w² [µm²]"); ax[1].legend(fontsize=8)
ax[2].loglog(q / kk, Sblk.mean(0), ".", ms=3); ax[2].set_xlabel("q [1/µm]"); ax[2].set_ylabel("S(q)")
ax[2].set_title(f"pieces of {Lc * kk:.0f} µm")
fig.suptitle(args.out); fig.tight_layout(); fig.savefig(os.path.join(args.out, "fig_walls.png"), dpi=90)
print(f"{args.out}: {nw} walls, {lens.min() * kk:.0f}-{Lmax * kk:.0f} um; zeta (local width, {rmin * kk:.1f}-{rmax * kk:.1f} um): "
      f"PCA {z_pca:.2f} +- {dz_pca:.2f}, no rotation {z_ctl:.2f}")
