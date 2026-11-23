"""Animated GIF of a vincent_walls.py segmentation (CoFeB MAR007 difference movies).

Usage: python cofeb_gif.py DIR SRC_DIR [--dark-npy FILE] [--every 4] [--scale 0.3] [--fps-gif 10]
Input: DIR/meta.json and DIR/tiff_extra.npz (vincent_walls.py) and the TIFF frames SRC_DIR/<seq>_*.tif.
Each GIF frame: the full camera frame minus the dark reference (Gaussian 1.5 px, scaled by --scale), the analysed wall
position (x_eff + subpixel part, as in h_sub; rows used only) in red, the analysed rows marked by yellow ticks at the left border,
time, field and a 20 um scale bar. Output: DIR/segmentation.gif
"""
import argparse
import glob
import json
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi

ap = argparse.ArgumentParser()
ap.add_argument("dir"); ap.add_argument("src")
ap.add_argument("--dark-npy", default=None)
ap.add_argument("--every", type=int, default=4); ap.add_argument("--scale", type=float, default=0.3)
ap.add_argument("--fps-gif", type=float, default=10)
a = ap.parse_args()

m = json.load(open(os.path.join(a.dir, "meta.json")))
e = np.load(os.path.join(a.dir, "tiff_extra.npz"))
# the analysed front: x_eff (switched area per row, h = W-1-x_eff) plus the subpixel part of the first-gap front, i.e. the
# position used for h_sub (a bright speck inside the domain moves the first-gap front, not this one)
xs, rows = e["x_eff"] + (e["x_wall_sub"] - e["x_wall"]), e["rows"]
fs = sorted(glob.glob(os.path.join(a.src, f"{m['seq']}_*.tif")))[:m["n_frames_video"]]
K = len(fs)
flipped = not m["geometry"].startswith("image not flipped")
dark = np.load(a.dark_npy).astype(np.float32) if a.dark_npy else 0.0
ld = lambda k: ndi.gaussian_filter(np.asarray(Image.open(fs[k])).astype(np.float32) - dark, 1.5)
H, W = ld(0).shape
lo, hi = np.percentile(ld(0)[::4, ::4], [1, 99.5])
s = a.scale
Wg, Hg = int(W * s), int(H * s)
try:
    font = ImageFont.truetype("DejaVuSans-Bold.ttf", 15)
except OSError:
    font = ImageFont.load_default()
idx = list(range(0, K, a.every)) + ([K - 1] if (K - 1) % a.every else [])
frames = []
for k in idx:
    g = np.clip((ld(k) - lo) / (hi - lo) * 255, 0, 255).astype(np.uint8)
    im = Image.fromarray(g).resize((Wg, Hg), Image.BILINEAR).convert("RGB")
    d = ImageDraw.Draw(im)
    x = xs[k]
    ok = np.isfinite(x)
    xo = (W - 1 - x) if flipped else x
    # draw contiguous runs only (no lines across gaps)
    br = np.nonzero(np.diff(np.nonzero(ok)[0]) > 1)[0]
    for seg in np.split(np.nonzero(ok)[0], br + 1):
        if seg.size > 1:
            d.line(list(zip(xo[seg] * s, rows[seg] * s)), fill=(255, 40, 40), width=2)
    d.line([(2, rows[0] * s), (12, rows[0] * s)], fill=(255, 230, 0), width=2)
    d.line([(2, rows[-1] * s), (12, rows[-1] * s)], fill=(255, 230, 0), width=2)
    t = k / m["fps"]
    lab = f"{m['seq']}  {m.get('field_uT', 0):.1f} uT  frame {k}/{K - 1}  t = {t:.0f} s"
    d.text((8, 6), lab, fill=(255, 230, 0), font=font, stroke_width=2, stroke_fill=(0, 0, 0))
    L = 20 / m["um_per_px"] * s
    d.line([(Wg - 20 - L, Hg - 20), (Wg - 20, Hg - 20)], fill=(255, 255, 255), width=4)
    d.text((Wg - 20 - L, Hg - 42), "20 um", fill=(255, 255, 255), font=font, stroke_width=2, stroke_fill=(0, 0, 0))
    frames.append(im.quantize(colors=128, dither=Image.Dither.NONE))
out = os.path.join(a.dir, "segmentation.gif")
frames[0].save(out, save_all=True, append_images=frames[1:], duration=int(1000 / a.fps_gif), loop=0, optimize=True)
print(out, len(frames), "frames", f"{os.path.getsize(out) / 1e6:.1f} MB", flush=True)
