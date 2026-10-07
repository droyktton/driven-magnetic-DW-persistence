"""Step 5: pinning map (TIFF studies): local wall velocity, lagunas and per-row persistence.

Usage: python step5_pinning_map.py DIR [--eps 0.97] [--tau T] [--smin 20] [--sigma 4]
Input: DIR/tiff_extra.npz, DIR/h_xt_sub.npy and DIR/meta.json (step1).
Output: DIR/fig_pinning_map.png and a summary on stdout.

  local velocity  v(x,y) = 1/|grad t_arrival| (px/frame), from the arrival-time map smoothed
                  with a Gaussian of --sigma px. Slow regions are where the wall was held back.
  lagunas         swept pixels that switched in 1-frame patches smaller than --smin px (as in step4).
  per-row values  m_i = <p_i(t,tau)>_t, the mean persistence of row y_i at lag --tau and threshold
                  --eps (the static term removed by step2 --row-mean); v_i, the mean velocity of
                  the row (total advance / duration); the laguna fraction of the row.
--tau defaults to tau* (the lag of the maximum of chi4) for the given eps.
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from scipy.ndimage import binary_dilation, binary_erosion, gaussian_filter, label

ap = argparse.ArgumentParser()
ap.add_argument("dir", help="step1 output folder of a TIFF study")
ap.add_argument("--eps", type=float, default=0.97, help="persistence threshold in px")
ap.add_argument("--tau", type=int, default=None, help="lag in frames (default: tau* for eps)")
ap.add_argument("--smin", type=float, default=20, help="minimum avalanche area in px (lagunas below)")
ap.add_argument("--sigma", type=float, default=4, help="smoothing of the arrival map for the velocity, px")
args = ap.parse_args()

with open(os.path.join(args.dir, "meta.json")) as f:
    meta = json.load(f)
k_um, fps = meta.get("um_per_px") or 0.0, meta["fps"]
e = np.load(os.path.join(args.dir, "tiff_extra.npz"))
ta, defects, xe, y_off = e["arrival"], e["defects"], e["x_eff"], int(e["rows"][0])
h = np.load(os.path.join(args.dir, "h_xt_sub.npy"))
T, L = h.shape
H, W = ta.shape

# --- lagunas (1-frame patches smaller than smin)
fin = np.isfinite(ta) & (ta >= 1)
K = np.where(fin, np.round(ta), -1).astype(int)
size = np.zeros(K.shape, int)
for k in range(1, K.max() + 1):
    mk = K == k
    if mk.any():
        lab, _ = label(mk, structure=np.ones((3, 3), bool))
        sz = np.bincount(lab.ravel())
        sz[0] = 0
        size[mk] = sz[lab[mk]]
ok = fin & ~defects
lag = ok & (size < args.smin)

# --- local velocity
w = gaussian_filter(np.where(fin, ta, 0.0), args.sigma)
n = gaussian_filter(fin.astype(float), args.sigma)
with np.errstate(divide="ignore", invalid="ignore"):
    t = w / n
    gy, gx = np.gradient(t)
    v = 1 / np.hypot(gx, gy)
inner = binary_erosion(fin, iterations=2 * int(args.sigma)) & ~binary_dilation(defects, iterations=2 * int(args.sigma))
v_av, v_lag = np.median(v[inner & ~lag & ok]), np.median(v[inner & lag])

# --- per-row values
if args.tau is None:
    chi4 = [L * (np.abs(h[tau:] - h[:-tau]) < args.eps).mean(1).var() for tau in range(1, min(100, T // 2))]
    args.tau = int(np.argmax(chi4)) + 1
p = np.abs(h[args.tau:] - h[:-args.tau]) < args.eps
m = p.mean(0)
v_row = (xe[-1] - xe[0]) / (T - 1)
lag_row = lag.sum(1) / np.maximum(ok.sum(1), 1)
clean = ~binary_dilation(defects, iterations=3).any(1)        # rows whose path has no defect
r = lambda a, b, s: np.corrcoef(a[s], b[s])[0, 1]
dm = m - m.mean()
ac = np.array([np.mean(dm[:L - k] * dm[k:]) for k in range(L // 4)]) / dm.var()
l_m = int(np.argmax(ac < 1 / np.e))
dl = lag_row - lag_row.mean()
acl = np.array([np.mean(dl[:L - k] * dl[k:]) for k in range(L // 4)]) / dl.var()
l_lag = int(np.argmax(acl < 1 / np.e))

u = f" = {{:.2f}} µm" if k_um else ""
print(f"local velocity (smoothing {args.sigma:g} px): median {v_av:.3f} px/frame in avalanche pixels, "
      f"{v_lag:.3f} in lagunas (ratio {v_lag / v_av:.2f})")
print(f"per row, eps = {args.eps:g} px, tau = {args.tau} frames: m_i = {m.mean():.3f} ± {m.std():.3f}; "
      f"v_i = {v_row.mean():.3f} ± {v_row.std():.3f} px/frame; laguna fraction {lag_row.mean():.2f} ± {lag_row.std():.2f}")
print(f"correlations (all rows / rows without defects, n = {clean.sum()}): "
      f"m-v {r(m, v_row, slice(None)):+.2f} / {r(m, v_row, clean):+.2f}; "
      f"m-lagunas {r(m, lag_row, slice(None)):+.2f} / {r(m, lag_row, clean):+.2f}; "
      f"v-lagunas {r(v_row, lag_row, slice(None)):+.2f} / {r(v_row, lag_row, clean):+.2f}")
print(f"correlation length along y (1/e): m_i {l_m} px" + (u.format(l_m * k_um) if k_um else "")
      + f"; laguna fraction {l_lag} px" + (u.format(l_lag * k_um) if k_um else ""))

# --- figure
y = np.arange(L) + y_off
fig = plt.figure(figsize=(18, 10))
ax = fig.add_subplot(1, 2, 1)
im = ax.imshow(np.where(inner | (fin & ~defects), v, np.nan), cmap="viridis",
               norm=LogNorm(vmin=0.05, vmax=2), extent=(0, W, H + y_off, y_off), interpolation="nearest")
ax.imshow(np.where(defects, 0.0, np.nan), cmap="gray", vmin=0, vmax=1, extent=(0, W, H + y_off, y_off))
plt.colorbar(im, ax=ax, label="local wall velocity [px/frame]", shrink=0.8)
ax.set_title(f"Local velocity 1/|∇t_arrival| (smoothing {args.sigma:g} px); black: defects", fontsize=10)
ax.set_xlabel("x [px]"); ax.set_ylabel("y [px]")
for j, (vals, lab_, c) in enumerate(((m, rf"$m_i$ ($\varepsilon$={args.eps:g} px, $\tau$={args.tau})", "C0"),
                                       (v_row, r"$v_i$ [px/frame]", "C1"),
                                       (lag_row, "laguna fraction", "C2"))):
    a = fig.add_subplot(3, 4, 3 + 4 * j)
    a.plot(y, vals, color=c, lw=0.8)
    a.plot(y[~clean], vals[~clean], ".", ms=1.5, color="k")
    a.set_ylabel(lab_); a.set_xlim(y[0], y[-1])
    if j == 2:
        a.set_xlabel("y [px] (black: rows with a defect in the path)")
a = fig.add_subplot(3, 4, 4)
a.plot(v_row[clean], m[clean], ".", ms=2)
a.set_xlabel(r"$v_i$ [px/frame]"); a.set_ylabel(r"$m_i$"); a.set_title(f"r = {r(m, v_row, clean):+.2f}", fontsize=9)
a = fig.add_subplot(3, 4, 8)
a.plot(lag_row[clean], m[clean], ".", ms=2, color="C2")
a.set_xlabel("laguna fraction"); a.set_ylabel(r"$m_i$"); a.set_title(f"r = {r(m, lag_row, clean):+.2f}", fontsize=9)
a = fig.add_subplot(3, 4, 12)
a.plot(ac, label=r"$m_i$"); a.plot(acl, label="laguna fraction")
a.axhline(1 / np.e, color="k", lw=0.5, ls="--")
a.set_xlabel("n [px] along y"); a.set_ylabel("autocorrelation"); a.legend(fontsize=8)
fig.tight_layout(); fig.savefig(os.path.join(args.dir, "fig_pinning_map.png"), dpi=110); plt.close(fig)
print(f"output: {args.dir}/fig_pinning_map.png")
