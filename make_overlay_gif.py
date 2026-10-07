"""Animated GIF of a TIFF study: registered frames with the detected wall, for sharing.

Usage: python make_overlay_gif.py DIR TIFF_DIR [--every 5] [--scale 0.5] [--fps-gif 12]
                                  [--out DIR/overlay_movie.gif]
Input: DIR/tiff_extra.npz and DIR/meta.json (step1) and the original TIFF frames.
Each GIF frame is a registered (drift-corrected), cropped camera frame, scaled by --scale, with
the wall front in red and the elapsed time. Every --every-th frame is used.
"""
import argparse
import glob
import json
import os
import re
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import tifffile
from PIL import Image, ImageDraw, ImageFont
from scipy.ndimage import gaussian_filter
from scipy.ndimage import shift as nd_shift

ap = argparse.ArgumentParser()
ap.add_argument("dir", help="step1 output folder of a TIFF study")
ap.add_argument("tiff_dir", help="folder with the original TIFF frames")
ap.add_argument("--every", type=int, default=5, help="use every N-th frame")
ap.add_argument("--scale", type=float, default=0.5, help="image scale factor")
ap.add_argument("--fps-gif", type=float, default=12, help="GIF playback rate (frames per second)")
ap.add_argument("--out", default=None, help="output file (default DIR/overlay_movie.gif)")
args = ap.parse_args()

with open(os.path.join(args.dir, "meta.json")) as f:
    meta = json.load(f)
fps, k_um = meta["fps"], meta.get("um_per_px")
e = np.load(os.path.join(args.dir, "tiff_extra.npz"))
drift, rows, xc, xw = e["drift_smooth"], e["rows"], int(e["crop_right"]), e["x_wall_sub"]
y0, y1 = int(rows[0]), int(rows[-1]) + 1
files = sorted((f for f in glob.glob(os.path.join(args.tiff_dir, "*")) if re.search(r"\.tiff?$", f, re.I)),
               key=lambda f: [int(n) for n in re.findall(r"\d+", os.path.basename(f))])
T = len(drift)
idx = list(range(0, T, args.every)) + ([T - 1] if (T - 1) % args.every else [])


def frame(k):
    im = tifffile.imread(files[k]).astype(np.float32)
    im = nd_shift(gaussian_filter(im, 1.5), drift[k], order=1, mode="nearest")[y0:y1, :xc]
    return im / np.median(im)


lo, hi = np.percentile(frame(0), [1, 99.5])
W, H = int(round(xc * args.scale)), int(round((y1 - y0) * args.scale))
try:
    font = ImageFont.truetype("DejaVuSans-Bold.ttf", max(14, int(32 * args.scale)))
except OSError:
    font = ImageFont.load_default()
bar_px = 20 / k_um if k_um else None           # 20 µm scale bar
# fixed palette: 120 grey levels + the overlay colours, so the thin red line survives quantization
NG = 120
RED, YELLOW, WHITE, BLACK = (255, 30, 30), (255, 230, 0), (255, 255, 255), (0, 0, 0)
pal = [c for i in range(NG) for c in (round(i * 255 / (NG - 1)),) * 3] + [*RED, *YELLOW, *WHITE, *BLACK]
pal_img = Image.new("P", (1, 1))
pal_img.putpalette(pal + [0] * (768 - len(pal)))


def render(k):
    g = np.clip((frame(k) - lo) / (hi - lo) * 255, 0, 255).astype(np.uint8)
    im = Image.fromarray(g).resize((W, H), Image.BILINEAR).convert("RGB")
    d = ImageDraw.Draw(im)
    xs = xw[k]
    ok = np.isfinite(xs)
    pts = list(zip(xs[ok] * args.scale, np.nonzero(ok)[0] * args.scale))
    d.line(pts, fill=RED, width=2)
    t = k / fps
    label = f"t = {int(t // 3600)} h {int(t % 3600 // 60):02d} min"
    box = d.textbbox((8, 6), label, font=font)
    d.rectangle((box[0] - 4, box[1] - 3, box[2] + 4, box[3] + 3), fill=BLACK)
    d.text((8, 6), label, fill=YELLOW, font=font)
    if bar_px:
        L = bar_px * args.scale
        tb = d.textbbox((0, 0), "20 µm", font=font)
        d.rectangle((W - 28 - L, H - 26 - (tb[3] - tb[1]) - 8, W - 12, H - 12), fill=BLACK)
        d.line([(W - 20 - L, H - 20), (W - 20, H - 20)], fill=WHITE, width=4)
        d.text((W - 20 - L, H - 26 - (tb[3] - tb[1])), "20 µm", fill=WHITE, font=font)
    return im.quantize(palette=pal_img, dither=Image.Dither.NONE)


with ThreadPoolExecutor(6) as ex:
    ims = list(ex.map(render, idx))
out = args.out or os.path.join(args.dir, "overlay_movie.gif")
ims[0].save(out, save_all=True, append_images=ims[1:], duration=int(1000 / args.fps_gif), loop=0, optimize=True)
print(f"{len(ims)} frames ({W}x{H}), {os.path.getsize(out) / 1e6:.1f} MB -> {out}")
