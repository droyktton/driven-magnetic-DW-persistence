"""Arrival map of a MOKE movie: for each pixel, the frame at which it switches.

Usage: python arrival_map.py MOVIE --outdir DIR [--sign 1] [--sigma 1] [--detrend] [--frames F0 F1] [--qc K1,K2,...]
MOVIE: any video OpenCV can read (avi, mov, mp4), converted to grey.
For every pixel, the frame of the largest step of its grey level (step-function cross-correlation, as in
Durin et al., PRB 110, L020405 (2024)): maximise s * (mean after - mean before) * sqrt(n1 n2 / N), s = --sign
(+1: switching makes the pixel brighter). Frames smoothed by a Gaussian of --sigma px. Pixels whose step amplitude
is below the Otsu threshold of all amplitudes never switched (-1).
--detrend removes a global brightness change (camera gain, lamp): every frame is shifted by the cumulative median of
the frame-to-frame change over the aperture (few pixels switch between consecutive frames).
Outputs: DIR/arrival.npz (ta: int16 frame of switching or -1; amp: step amplitude; thr; mean_frame, for the aperture), DIR/fig_qc_overlay.png
(arrival map and the switched region {0 <= ta <= k} over frame k for the --qc frames), DIR/meta.json.
"""
import argparse
import json
import os

import cv2
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from skimage.filters import threshold_otsu

ap = argparse.ArgumentParser()
ap.add_argument("movie")
ap.add_argument("--outdir", required=True)
ap.add_argument("--sign", type=int, default=1, help="+1 if switching makes the pixel brighter, -1 darker")
ap.add_argument("--sigma", type=float, default=1.0, help="spatial smoothing of each frame, px")
ap.add_argument("--detrend", action="store_true", help="remove a global brightness change frame by frame")
ap.add_argument("--frames", type=int, nargs=2, help="first and last frame used (default: all)")
ap.add_argument("--qc", default="", help="comma-separated frames for the overlay figure")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)

cap = cv2.VideoCapture(args.movie)
fps = cap.get(cv2.CAP_PROP_FPS)
raw = []
while True:
    ok, im = cap.read()
    if not ok:
        break
    raw.append(cv2.cvtColor(im, cv2.COLOR_BGR2GRAY))
raw = np.array(raw)
f0, f1 = args.frames if args.frames else (0, len(raw) - 1)
fr = raw[f0:f1 + 1].astype(np.float32)
T, H, W = fr.shape
off = np.zeros(T)
if args.detrend:
    aper = fr.mean(0) > 0.3 * fr.mean()
    off = np.concatenate([[0], np.cumsum([np.median((fr[t] - fr[t - 1])[aper]) for t in range(1, T)])])
    fr -= off[:, None, None].astype(np.float32)
if args.sigma > 0:
    fr = np.array([gaussian_filter(f, args.sigma) for f in fr])
ta = np.zeros((H, W), np.int16)
amp = np.zeros((H, W), np.float32)
n1 = np.arange(1, T, dtype=np.float32)[:, None, None]
n2 = T - n1
for r0 in range(0, H, 64):
    cs = np.cumsum(fr[:, r0:r0 + 64], 0)
    d = ((cs[-1] - cs[:-1]) / n2 - cs[:-1] / n1) * args.sign
    k = np.argmax(d * np.sqrt(n1 * n2 / T), 0)
    ta[r0:r0 + 64] = k + 1 + f0                     # first frame after the step
    amp[r0:r0 + 64] = np.take_along_axis(d, k[None], 0)[0]
thr = float(threshold_otsu(amp))
sw = amp > thr
ta = np.where(sw, ta, -1).astype(np.int16)
np.savez_compressed(os.path.join(args.outdir, "arrival.npz"), ta=ta, amp=amp.astype(np.float16), thr=thr,
                    mean_frame=np.round(raw[f0:f1 + 1].mean(0)).astype(np.uint8))
meta = dict(movie=os.path.basename(args.movie), n_frames_video=len(raw), frames=[f0, f1], fps_header=fps, sign=args.sign,
            sigma=args.sigma, detrend=args.detrend, detrend_offset_range=[float(off.min()), float(off.max())],
            otsu_threshold=thr, amp_median_switched=float(np.median(amp[sw])), amp_median_not=float(np.median(amp[~sw])),
            frac_switched=float(sw.mean()), shape=[H, W])
json.dump(meta, open(os.path.join(args.outdir, "meta.json"), "w"), indent=1)
print(f"{meta['movie']}: frames {f0}-{f1}, Otsu threshold {thr:.1f} (median amplitude switched {meta['amp_median_switched']:.1f}, "
      f"not {meta['amp_median_not']:.1f}); switched {sw.mean():.1%}" + (f"; detrend offset {off.min():.1f}..{off.max():.1f}" if args.detrend else ""))

qc = [int(v) for v in args.qc.split(",") if v]
if qc:
    nc = (len(qc) + 2) // 2
    fig, ax = plt.subplots(2, nc, figsize=(6.2 * nc, 9.5))
    ax = ax.ravel()
    im = ax[0].imshow(np.where(ta >= 0, ta, np.nan), cmap="turbo")
    ax[0].set_title("arrival frame (blank = never)")
    plt.colorbar(im, ax=ax[0], fraction=0.04)
    for a, k in zip(ax[1:], qc):
        a.imshow(raw[k], cmap="gray")
        a.contour(((ta >= 0) & (ta <= k)).astype(float), [0.5], colors="r", linewidths=0.7)
        a.set_title(f"frame {k}: switched by the arrival map (red)")
    for a in ax:
        a.axis("off")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "fig_qc_overlay.png"), dpi=50)
