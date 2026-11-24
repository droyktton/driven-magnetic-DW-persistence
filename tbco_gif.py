"""Animated GIF of the TbCo pulse segmentation (output of tbco_pulses.py), one GIF frame per image file.

Usage: python tbco_gif.py SRC_DIR PREFIX DIR [--scale 0.25] [--ms 700] [--tb 0.83nm]
Each frame: the registered raw image (contrast stretched), the walls found in it (red) and in the previous image
(cyan, only if both images share a diaphragm state), and the file number, diaphragm state, field and pulse length.
Output: DIR/segmentation.gif
"""
import argparse
import glob
import json
import os

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("prefix"); ap.add_argument("dir")
ap.add_argument("--scale", type=float, default=0.25); ap.add_argument("--ms", type=int, default=700)
ap.add_argument("--tb", default="")
a = ap.parse_args()

J = json.load(open(os.path.join(a.dir, "pulses.json")))
state = {int(k): v for k, v in J["state"].items()}
shift = {int(k): v for k, v in J["shift"].items()}
info = {p["file"]: p for p in J["pulses"]}
Wz = np.load(os.path.join(a.dir, "walls.npz"))
walls = {}
for f in Wz.files:
    k = int(f[1:].split("_")[0]); walls.setdefault(k, []).append(Wz[f])
files = {int(f.rsplit("_", 1)[1][:-4]): f for f in glob.glob(os.path.join(a.src, a.prefix + "_*.tif"))
         if f.rsplit("_", 1)[1][:-4].isdigit()}
ks = [k for k in sorted(files) if k in shift]
try:
    font = ImageFont.truetype("DejaVuSans.ttf", 15)
except OSError:
    font = ImageFont.load_default()
s = a.scale
frames = []
for k in ks:
    im = np.asarray(Image.open(files[k])).astype(np.float32)
    if np.hypot(*shift[k]) > 0.5:
        im = ndi.shift(im, shift[k], order=1, mode="nearest")
    im = ndi.zoom(ndi.gaussian_filter(im, 1 / (2 * s)), s, order=1)
    lo, hi = np.percentile(im, [1, 99.7])
    g = (np.clip((im - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)
    img = Image.fromarray(g).convert("RGB"); d = ImageDraw.Draw(img)
    if k - 1 in walls and state.get(k - 1) == state.get(k):
        for c in walls[k - 1]:
            d.line([(x * s, y * s) for y, x in c[::2]], fill=(0, 220, 255), width=1)
    for c in walls.get(k, []):
        d.line([(x * s, y * s) for y, x in c[::2]], fill=(255, 40, 40), width=2)
    t = f"TbCo {a.tb}  file {k}  diaphragm state {state.get(k, '-')}"
    p = info.get(k)
    t2 = f"pulse: {p['H']:.0f} mT, {p['dur_us']:g} us" if p else "(nucleation / no pulse data)"
    for i, tt in enumerate((t, t2, "red: wall after pulse   cyan: before")):
        d.text((8, 6 + 18 * i), tt, fill=(255, 255, 80), font=font, stroke_width=2, stroke_fill=(0, 0, 0))
    frames.append(img)
out = os.path.join(a.dir, "segmentation.gif")
frames[0].save(out, save_all=True, append_images=frames[1:], duration=a.ms, loop=0, optimize=True)
print(out, len(frames), "frames", f"{os.path.getsize(out) / 1e6:.1f} MB")
