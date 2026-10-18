"""Arrival map from a TIFF folder for low-contrast, slow walls (where step1's domain thresholding fails).

Usage: python tiff_arrival_map.py TIFF_DIR --outdir DIR [--drift DIR/drift_raw.npy] [--offset 33000] [--block 10]
       [--x0 0 --x1 W --y0 0 --y1 H] [--sign 1] [--sigma 1.5]
Frames are sorted by the number at the end of the file name, offset-subtracted, registered with the drift of an earlier
step1 run (drift_raw.npy, per frame (dy, dx), smoothed over 3 frames and applied as in step1; optional), cropped, divided by their median (global lamp drift) and
averaged in blocks of --block frames. For every pixel, the block of the largest step of log I (step-function
cross-correlation, as arrival_map.py: s * (mean after - mean before) * sqrt(n1 n2 / N), s = --sign) after a Gaussian
of --sigma px. Pixels whose step amplitude is below the Otsu threshold never switched (-1).
Outputs: DIR/arrival.npz (ta in FRAMES = block start, -1 below threshold; ta_all, the same for every pixel; amp, thr, mean_frame = 255 inside the crop; block, crop),
DIR/fig_arrival.png.
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
from PIL import Image
from scipy.ndimage import gaussian_filter, gaussian_filter1d, shift
from skimage.filters import threshold_otsu

ap = argparse.ArgumentParser()
ap.add_argument("tiffdir")
ap.add_argument("--outdir", required=True)
ap.add_argument("--drift", help="drift_raw.npy of a step1 run on the same folder")
ap.add_argument("--offset", type=float, default=33000)
ap.add_argument("--block", type=int, default=10)
ap.add_argument("--x0", type=int, default=0); ap.add_argument("--x1", type=int, default=None)
ap.add_argument("--y0", type=int, default=0); ap.add_argument("--y1", type=int, default=None)
ap.add_argument("--sign", type=int, default=1, help="+1 if switching makes the pixel brighter")
ap.add_argument("--sigma", type=float, default=1.5)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

fs = sorted(glob.glob(os.path.join(a.tiffdir, "*.tif*")), key=lambda f: int(re.search(r"(\d+)\.tiff?$", f).group(1)))
drift = gaussian_filter1d(np.load(a.drift), 3, axis=0, mode="nearest") if a.drift else np.zeros((len(fs), 2))
H, W = np.asarray(Image.open(fs[0])).shape
x1, y1 = a.x1 or W, a.y1 or H
nb = len(fs) // a.block
stack = np.zeros((nb, y1 - a.y0, x1 - a.x0), np.float32)
for b in range(nb):
    acc = 0
    for i in range(b * a.block, (b + 1) * a.block):
        im = np.asarray(Image.open(fs[i]), np.float32) - a.offset
        if np.any(drift[i]):
            im = shift(im, drift[i], order=1, mode="nearest")   # as step1 (smoothed drift)
        c = im[a.y0:y1, a.x0:x1]
        acc = acc + c / np.median(c)
    stack[b] = gaussian_filter(np.log(np.maximum(acc / a.block, 1e-3)), a.sigma)
    if b % 50 == 0:
        print(f"block {b}/{nb}", flush=True)

# step detection: cut after block k (k = 1..nb-1), cumulative sums over time
cs = np.cumsum(stack, 0, dtype=np.float64); tot = cs[-1]
best = np.full(stack.shape[1:], -np.inf); ta = np.zeros(stack.shape[1:], np.int16)
for k in range(1, nb):
    n1, n2 = k, nb - k
    s = a.sign * ((tot - cs[k - 1]) / n2 - cs[k - 1] / n1) * np.sqrt(n1 * n2 / nb)
    m = s > best; best[m] = s[m]; ta[m] = k
thr = threshold_otsu(best[np.isfinite(best)])
ta_all = (ta * a.block).astype(np.int16)                    # frame of the largest step, every pixel
ta = np.where(best > thr, ta_all, -1).astype(np.int16)
np.savez_compressed(os.path.join(a.outdir, "arrival.npz"), ta=ta, ta_all=ta_all, amp=best.astype(np.float16), thr=thr,
                    mean_frame=np.full(ta.shape, 255, np.uint8), block=a.block, crop=[a.x0, x1, a.y0, y1])
json.dump(dict(tiffdir=a.tiffdir, frames=len(fs), block=a.block, crop=[a.x0, x1, a.y0, y1], sign=a.sign,
               sigma=a.sigma, thr=float(thr), switched=float((ta >= 0).mean())),
          open(os.path.join(a.outdir, "meta.json"), "w"), indent=1)

fig, ax = plt.subplots(1, 3, figsize=(18, 7))
ax[0].imshow(stack[-1] - stack[0], cmap="gray", vmin=-0.03, vmax=0.03); ax[0].set_title("log(last block / first block)")
ax[1].imshow(best, cmap="magma", vmin=0, vmax=np.percentile(best, 99.5)); ax[1].set_title(f"step amplitude (Otsu {thr:.3g})")
im = ax[2].imshow(np.ma.masked_less(ta, 0), cmap="viridis"); ax[2].set_title("arrival frame"); plt.colorbar(im, ax=ax[2])
plt.tight_layout(); plt.savefig(os.path.join(a.outdir, "fig_arrival.png"), dpi=80)
print(f"switched fraction {(ta >= 0).mean():.3f}, Otsu threshold {thr:.3g}")
