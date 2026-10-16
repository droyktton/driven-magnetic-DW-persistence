"""Roughness and overhangs of critical phi^4 domain walls (vmc-phi4-depinning, --dump-field), straight or bubble.

Usage: python phi4_contour_roughness.py DIR [DIR2 ...] --L L [--Ly LY] [--bubble R0] [--block B] [--outdir OUT]
DIR holds criticafield_h*_seed*.bin (raw float32, Lx*Ly, row-major, index x + y*Lx) written by phi4vmc --dump-field.
The wall is the phi = 0 line (marching squares, subpixel). Straight wall (periodic in x): the field is padded by
64 columns on each side, the longest contour kept and cut to 0 <= x < Lx. Bubble (--bubble R0, the initial
radius): the longest closed contour; points within 3 sites of the initial circle r = R0 (wall that never moved)
are dropped. With --block B the field is first reduced to B x B blocks (sign of the block mean: a binary 'pixel'
image, as a camera sees it) and the contour is taken of that mask smoothed by 1 px, as contour_analysis.py does;
lengths are then in blocks.
Per sample: local width of arc windows rotated to their own axis and P_ov(E) (wall_tools.overhang_segments),
minus the variance of a perfect arc of the mean radius for a bubble. Averaged over samples; zeta from w2(E) over
several ranges, with a bootstrap over samples (16-84 %). Outputs OUT/phi4_roughness.npz and
OUT/fig_phi4_roughness.png; table on stdout.
"""
import argparse
import glob
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter
from skimage.measure import find_contours

from wall_tools import fit_power, overhang_segments

ap = argparse.ArgumentParser()
ap.add_argument("dirs", nargs="+")
ap.add_argument("--L", type=int, required=True, help="Lx")
ap.add_argument("--Ly", type=int, default=0, help="Ly (default: from the file size)")
ap.add_argument("--bubble", type=float, default=0, help="initial radius R0 of the bubble geometry")
ap.add_argument("--block", type=int, default=1, help="reduce to B x B blocks (binary pixels)")
ap.add_argument("--outdir", default=None)
ap.add_argument("--n-boot", type=int, default=300)
args = ap.parse_args()
outdir = args.outdir or args.dirs[0]
os.makedirs(outdir, exist_ok=True)
rng = np.random.default_rng(5)
B = args.block
step = 0.5
TOLS = [2.0, 3.0, 5.0]
RANGES = [(8, 80), (20, 80), (8, 40), (40, 200), (100, 400)]


def resample(c):
    seg = np.hypot(*np.diff(c, axis=0).T)
    arc = np.concatenate([[0], np.cumsum(seg)])
    s = np.arange(0, arc[-1], step)
    return np.interp(s, arc, c[:, 0]), np.interp(s, arc, c[:, 1])


def runs(ok):
    idx = np.nonzero(ok)[0]
    if idx.size == 0:
        return []
    return np.split(idx, np.nonzero(np.diff(idx) > 1)[0] + 1)


def contours(phi):
    Ly, Lx = phi.shape
    if B > 1:
        H, W = (Ly // B) * B, (Lx // B) * B
        m = phi[:H, :W].reshape(H // B, B, W // B, B).mean((1, 3)) > 0
        field = gaussian_filter(m.astype(float), 1.0) - 0.5
    else:
        field = phi
    ny, nx = field.shape
    pieces, R = [], None
    if args.bubble:
        c = max(find_contours(field, 0.0), key=len)
        y, x = resample(c)
        cy, cx = (ny - 1) / 2, (nx - 1) / 2
        r = np.hypot(y - cy, x - cx)
        R = r.mean()
        keep = r > args.bubble / B + 3
        rr = runs(keep)
        if len(rr) > 1 and rr[0][0] == 0 and rr[-1][-1] == keep.size - 1:
            rr[0] = np.concatenate([rr[-1], rr[0]])
            rr = rr[:-1]
        if keep.all():
            pieces = [(y, x)]
        else:
            pieces = [(y[i], x[i]) for i in rr if i.size * step >= 60]
    else:
        pad = 64
        f = np.concatenate([field[:, -pad:], field, field[:, :pad]], axis=1)
        c = max(find_contours(f, 0.0), key=len)
        y, x = resample(c)
        x = x - pad
        ok = (x >= 0) & (x < nx)
        pieces = [(y[i], x[i]) for i in runs(ok) if i.size * step >= 60]
    return pieces, R


files = sorted(f for d in args.dirs for f in glob.glob(os.path.join(d, "criticafield_*.bin")))
if not files:
    raise SystemExit("no criticafield_*.bin files")
per = []
lam_max = None
data = []
for f in files:
    a = np.fromfile(f, dtype=np.float32)
    Ly = args.Ly or a.size // args.L
    phi = a.reshape(Ly, args.L)
    pcs, R = contours(phi)
    if not pcs:
        continue
    data.append((pcs, R))
    longest = max(p[0].size for p in pcs) * step
    lm = min(longest, R) if R else longest / 2
    lam_max = lm if lam_max is None else min(lam_max, lm)
lams = np.unique(np.round(np.logspace(np.log10(4), np.log10(0.95 * lam_max), 28)))
for pcs, R in data:
    NT = np.zeros(lams.size)
    NO = {t: np.zeros(lams.size) for t in TOLS}
    EX = [[] for _ in lams]
    W2 = [[] for _ in lams]
    for y, x in pcs:
        nt, no, ex, w2 = overhang_segments(y, x, lams, step, TOLS)
        NT += nt
        for t in TOLS:
            NO[t] += no[t]
        for j, lam in enumerate(lams):
            circ = 0.0
            if R:
                al = np.linspace(-lam / (2 * R), lam / (2 * R), 200)
                circ = np.var(R * np.cos(al))
            EX[j].extend(ex[j])
            W2[j].extend(np.asarray(w2[j]) - circ)
    per.append(dict(E=np.array([np.median(e) if len(e) else np.nan for e in EX]),
                    w2=np.array([np.mean(w) if len(w) else np.nan for w in W2]),
                    nt=NT, no={t: NO[t] for t in TOLS}, R=R))
n = len(per)
E = np.nanmedian([p["E"] for p in per], axis=0)
W2s = np.array([p["w2"] for p in per])
w2m = np.nanmean(W2s, 0)
P = {t: np.sum([p["no"][t] for p in per], 0) / np.maximum(np.sum([p["nt"] for p in per], 0), 1) for t in TOLS}


def zfit(w2, a, b):
    ok = (E >= a) & (E <= b) & (w2 > 0)
    return fit_power(E[ok], w2[ok], a, b)[0] / 2 if ok.sum() > 3 else np.nan


unit = "blocks" if B > 1 else "sites"
geom = f"bubble R0 = {args.bubble:g}" if args.bubble else "straight"
print(f"{geom}, {n} samples, block {B}: windows up to {lams.max():.0f} {unit}" +
      (f"; mean radius {np.mean([p['R'] for p in per]):.0f} {unit}" if args.bubble else ""))
rows = []
for a, b in RANGES:
    if b > 1.05 * np.nanmax(E):
        continue
    z = zfit(w2m, a, b)
    bs = [zfit(np.nanmean(W2s[rng.integers(0, n, n)], 0), a, b) for _ in range(args.n_boot)]
    lo, hi = np.nanpercentile(bs, [16, 84])
    rows.append((a, b, z, lo, hi))
    print(f"   zeta (local width), E = {a}-{b} {unit}: {z:.2f} [{lo:.2f}, {hi:.2f}]")
print("   w at E = 10/20/50: " + "/".join(f"{np.sqrt(np.interp(np.log(e), np.log(E[w2m > 0]), w2m[w2m > 0])):.2f}" for e in (10, 20, 50)))
print("   P_ov (tol 3) at E = 20/50/100/200: " + "/".join(f"{np.interp(np.log(e), np.log(E), P[3.0]):.3f}" if e <= np.nanmax(E) else "-" for e in (20, 50, 100, 200)))
lx, ly = np.log(E), np.log(np.where(w2m > 0, w2m, np.nan))
zloc = np.gradient(ly, lx) / 2
np.savez(os.path.join(outdir, "phi4_roughness.npz"), E=E, w2=w2m, w2_samples=W2s, zeta_local=zloc,
         lams=lams, ranges=np.array(rows), **{f"P_tol{t:g}": P[t] for t in TOLS}, files=np.array(files))

fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
ax[0].loglog(E, w2m, "o-", ms=3)
ax[0].set_xlabel(f"E [{unit}]")
ax[0].set_ylabel("⟨w²⟩" + (" − circle" if args.bubble else ""))
ax[0].set_title(f"{geom}, {n} samples")
ax[1].semilogx(E, zloc, "o-", ms=3)
for z in (0.5, 0.63, 1.25):
    ax[1].axhline(z, color="0.6", ls=":", lw=0.8)
ax[1].set_ylim(0, 1.8)
ax[1].set_xlabel(f"E [{unit}]")
ax[1].set_ylabel("local ζ = ½ dlog w²/dlog E")
for t, st in zip(TOLS, (":", "-", "--")):
    ax[2].semilogx(E, P[t], st, marker="o", ms=3, label=f"tol {t:g}")
ax[2].axhline(0.5, color="0.5", lw=0.8)
ax[2].set_ylim(-0.02, 1.02)
ax[2].legend()
ax[2].set_xlabel(f"E [{unit}]")
ax[2].set_ylabel("P_ov")
fig.tight_layout()
fig.savefig(os.path.join(outdir, "fig_phi4_roughness.png"), dpi=100)
