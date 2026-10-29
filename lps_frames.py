"""Segment the full LPS 2019 frames (GdFeCo B1742, 220-260 K): one stripe domain bounded by two walls, one frame per pulse.

Usage: python lps_frames.py SRC_DIR --outdir DIR [--sigma 3] [--um-per-px 0.171] [--fov-erode 40] [--min-area 2000]
Frames B1742_<T>K_<nnnnn>.tif, sorted by number (the number is the pulse index; missing numbers are frames Lucas discarded).
Per frame: Gaussian of --sigma px, Otsu threshold inside the field of view (the illuminated disc, from the mean frame),
largest bright component, holes filled. The stripe axis comes from the PCA of the union of all stripes; every frame is
rotated by it so the stripe runs along x. For every column of the rotated frame: top/bottom wall = first/last stripe pixel,
refined to subpixel at the threshold crossing; NaN where the column leaves the field of view or misses the stripe;
n_cross = number of stripe/background transitions in the column (> 2: overhang, hole or extra domain).
Outputs: DIR/frames.npz (h_top, h_bot [um] of shape (n_frames, n_x), n_cross, frame numbers, x [um], theta, thr),
DIR/contours.npz (subpixel contour of the stripe per frame, rotated frame, um), DIR/qc/overlay_<nnnnn>.png, DIR/fig_frames.png.
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
import tifffile
from scipy.ndimage import binary_erosion, binary_fill_holes, gaussian_filter, label, rotate
from skimage.filters import threshold_otsu
from skimage.measure import find_contours

ap = argparse.ArgumentParser()
ap.add_argument("src")
ap.add_argument("--outdir", required=True)
ap.add_argument("--sigma", type=float, default=3.0)
ap.add_argument("--um-per-px", type=float, default=0.171)
ap.add_argument("--fov-erode", type=int, default=40)
ap.add_argument("--min-area", type=int, default=2000)
a = ap.parse_args()
os.makedirs(os.path.join(a.outdir, "qc"), exist_ok=True)

fs = [f for f in glob.glob(os.path.join(a.src, "B1742_*.tif")) if re.search(r"_\d{5}\.tif$", f)]
fs = sorted(fs, key=lambda f: int(re.search(r"_(\d{5})\.tif$", f).group(1)))
nums = np.array([int(re.search(r"_(\d{5})\.tif$", f).group(1)) for f in fs])
ims = np.array([gaussian_filter(tifffile.imread(f).astype(np.float32), a.sigma) for f in fs])
nf, H, W = ims.shape
print(f"{nf} frames {H}x{W}, numbers {nums[0]}-{nums[-1]}", flush=True)

def stripe0(g):   # stripe with a whole-frame Otsu (first pass, only to exclude it from the background)
    return g > threshold_otsu(g[::4, ::4])


# field of view: the illuminated disc of the background (stripe pixels excluded; whole frame if no disc is visible)
bsum = np.zeros((H, W)); bcnt = np.zeros((H, W))
for g in ims:
    m = ~stripe0(g); bsum += g * m; bcnt += m
bg = gaussian_filter(np.where(bcnt > 0, bsum / np.maximum(bcnt, 1), np.nan), 1)
bg = np.where(np.isfinite(bg), bg, np.nanmedian(bg))
bg = gaussian_filter(bg, 30)
fov = bg > threshold_otsu(bg[::4, ::4])
if fov.mean() > 0.9 or fov.mean() < 0.2:
    fov = np.ones((H, W), bool)
fov = binary_fill_holes(fov)
fov = binary_erosion(fov, iterations=a.fov_erode) if a.fov_erode else fov
fov[:a.fov_erode // 4 + 1], fov[-(a.fov_erode // 4 + 1):], fov[:, :a.fov_erode // 4 + 1], fov[:, -(a.fov_erode // 4 + 1):] = 0, 0, 0, 0
print(f"field of view: {fov.mean():.2f} of the frame", flush=True)


def stripe(g):
    thr = threshold_otsu(g[fov])
    lab, n = label((g > thr) & fov)
    if n == 0:
        return np.zeros_like(fov), thr
    area = np.bincount(lab.ravel()); area[0] = 0
    return binary_fill_holes(lab == area.argmax()), thr


masks, thrs = [], []
for g in ims:
    m, t = stripe(g); masks.append(m); thrs.append(t)
masks = np.array(masks); thrs = np.array(thrs)

# stripe axis from the union of all stripes
yy, xx = np.nonzero(masks.any(0))
c = np.cov(np.vstack([xx, yy]))
ev, evec = np.linalg.eigh(c)
vx, vy = evec[:, -1]
theta = np.degrees(np.arctan2(vy, vx))   # angle of the axis in image coordinates (y down)
print(f"stripe axis {theta:.1f} deg", flush=True)

fov_r = rotate(fov.astype(np.float32), theta, reshape=True, order=1) > 0.99
Hr, Wr = fov_r.shape
h_top = np.full((nf, Wr), np.nan); h_bot = np.full((nf, Wr), np.nan); ncross = np.zeros((nf, Wr), int)
contours = {}
for k in range(nf):
    gr = rotate(ims[k], theta, reshape=True, order=1)
    mr = rotate(masks[k].astype(np.float32), theta, reshape=True, order=1) > 0.5
    thr = thrs[k]
    for x in range(Wr):
        col = mr[:, x]; inf = fov_r[:, x]
        if not col.any():
            continue
        r = np.nonzero(col)[0]; r0, r1 = r[0], r[-1]
        ncross[k, x] = int(np.abs(np.diff(col.astype(int))).sum())
        # both walls must lie inside the field of view, with background beyond them
        if r0 < 1 or r1 > Hr - 2 or not (inf[r0 - 1] and inf[r1 + 1]):
            continue
        gc = gr[:, x]
        f0 = (thr - gc[r0 - 1]) / (gc[r0] - gc[r0 - 1]) if gc[r0] != gc[r0 - 1] else 0.5
        f1 = (gc[r1] - thr) / (gc[r1] - gc[r1 + 1]) if gc[r1] != gc[r1 + 1] else 0.5
        h_top[k, x] = (r0 - 1 + np.clip(f0, 0, 1)) * a.um_per_px
        h_bot[k, x] = (r1 + np.clip(f1, 0, 1)) * a.um_per_px
    cs = find_contours(mr.astype(float), 0.5)
    contours[f"c{nums[k]:05d}"] = max(cs, key=len)[:, ::-1] * a.um_per_px if cs else np.zeros((0, 2))
    if k % 10 == 0:
        print(f"frame {k}/{nf}", flush=True)

xum = np.arange(Wr) * a.um_per_px
np.savez_compressed(os.path.join(a.outdir, "frames.npz"), h_top=h_top, h_bot=h_bot, n_cross=ncross, frame=nums,
                    x=xum, theta=theta, thr=thrs, um_per_px=a.um_per_px, sigma=a.sigma)
np.savez_compressed(os.path.join(a.outdir, "contours.npz"), **contours)

# QC overlays: every frame, original orientation, both walls
ct, st = np.cos(np.radians(theta)), np.sin(np.radians(theta))
cxr, cyr = (Wr - 1) / 2, (Hr - 1) / 2; cx, cy = (W - 1) / 2, (H - 1) / 2


def to_orig(xr, yr):   # inverse of scipy rotate(reshape=True) about the centre
    dx, dy = xr - cxr, yr - cyr
    return cx + ct * dx - st * dy, cy + st * dx + ct * dy


for k in range(nf):
    fig, ax = plt.subplots(figsize=(8, 8 * H / W))
    lo, hi = np.percentile(ims[k][fov], [1, 99.5])
    ax.imshow(ims[k], cmap="gray", vmin=lo, vmax=hi)
    for h, col in ((h_top[k], "tab:red"), (h_bot[k], "tab:cyan")):
        ok = np.isfinite(h)
        X, Y = to_orig(np.arange(Wr)[ok], h[ok] / a.um_per_px)
        ax.plot(X, Y, ".", ms=0.6, color=col)
    bad = (ncross[k] > 2) & np.isfinite(h_top[k])
    if bad.any():
        X, Y = to_orig(np.arange(Wr)[bad], h_top[k][bad] / a.um_per_px)
        ax.plot(X, Y, "x", ms=2, color="yellow")
    ax.contour(fov, [0.5], colors="w", linewidths=0.5)
    ax.set_title(f"{os.path.basename(a.src)}  frame {nums[k]}  (yellow: >2 crossings)", fontsize=9); ax.axis("off")
    fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "qc", f"overlay_{nums[k]:05d}.png"), dpi=110); plt.close(fig)

# summary: widths of the stripe, advance per pulse, valid fraction
wid = h_bot - h_top
valid = np.isfinite(wid).mean(1)
adv_t = np.nanmean(np.abs(np.diff(h_top, axis=0)), 1); adv_b = np.nanmean(np.abs(np.diff(h_bot, axis=0)), 1)
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
ax[0].plot(nums, np.nanmean(wid, 1)); ax[0].set_xlabel("frame"); ax[0].set_ylabel("mean stripe width [um]")
ax[1].plot(nums[1:], adv_t, ".-", label="top"); ax[1].plot(nums[1:], adv_b, ".-", label="bottom")
ax[1].set_xlabel("frame"); ax[1].set_ylabel("<|dh|> per frame step [um]"); ax[1].legend()
ax[2].plot(nums, valid); ax[2].set_xlabel("frame"); ax[2].set_ylabel("valid column fraction")
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig_frames.png"), dpi=100)
json.dump(dict(src=a.src, n_frames=int(nf), theta=float(theta), sigma=a.sigma, um_per_px=a.um_per_px,
               adv_top_median=float(np.nanmedian(adv_t)), adv_bot_median=float(np.nanmedian(adv_b)),
               valid_median=float(np.median(valid))), open(os.path.join(a.outdir, "frames_meta.json"), "w"), indent=1)
print(f"advance per step: top {np.nanmedian(adv_t):.2f} um, bottom {np.nanmedian(adv_b):.2f} um; valid {np.median(valid):.2f}")
