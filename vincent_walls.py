"""Central wall h(y,t) from the CoFeB MAR007 difference movies (V. Jeudy's group data, 9 Nov 2026, X5).

Usage: python vincent_walls.py SRC_DIR SEQ --outdir DIR [--last-frame K] [--current A] [--coil-mT-per-A 0.262]
                               [--dark-seq MAR007_seq13] [--um-per-px 0.536] [--sigma 2] [--core-erode 40]

The frames are Hamamatsu difference images (offset 2^15) of a ~0.9 mm polygonal disk. The dark (switched) state of
the whole disk is the mean of the last --dark-n frames of --dark-seq, where the disk has fully reversed; frame minus
dark leaves the domain contrast only (switched ~0, unswitched ~+120-160). No drift (< 0.1 px within and between
sequences), so no registration.

Per frame (no time averaging): Gaussian of --sigma px, threshold = Otsu of frame 0 inside the core (the disk, convex
hull of its border edges, eroded by --core-erode px so the dark rim along the border does not connect domains).
Dark components of >= --min-area px; the central-wall domain A_k = components overlapping A_{k-1} (A_0 = largest dark
component of frame 0). Dark components not connected to A (domains nucleated at the border) are ignored, so the
wall analysed is the central wall only; use --last-frame to stop before a nucleated domain merges with it.

The image is flipped left-right so that A is on the left and the wall moves to +x (step1 convention). Per row:
x_eff = last column of A if A were contiguous (left core edge + area of A in the row), front x_w = last A column before
the first gap, subpixel x_s = threshold crossing between x_w and x_w + 1. h = W-1-x_eff, h_sub = W-1-(x_eff + x_s-x_w)
(median over --median-x rows), as in step1. Rows used: the longest run of rows where the front stays inside the core
(--margin px from both core edges) in all but --max-fail of the frames; the failed frames are interpolated in time
(recorded in tiff_extra.npz: interpolated).

Outputs in DIR: h_xt.npy, h_xt_sub.npy, meta.json, tiff_extra.npz (arrival map, x_wall, x_wall_sub, x_eff, theta_deg,
n_overhang_rows, rows, core, threshold), fig_qc_overlay.png (original orientation), fig_kymograph.png,
fig_tilt_overhang.png. The arrival map uses step1's convention (-1 switched before frame 0, k first frame in A,
inf never), so step7 builds the same domain from it.
"""
import argparse
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
from skimage.morphology import convex_hull_image, remove_small_holes, remove_small_objects

ap = argparse.ArgumentParser()
ap.add_argument("src")
ap.add_argument("seq")
ap.add_argument("--outdir", required=True)
ap.add_argument("--last-frame", type=int, default=None, help="exclusive; default: all frames")
ap.add_argument("--current", type=float, default=None, help="coil current [A], for meta.json")
ap.add_argument("--coil-mT-per-A", type=float, default=0.262)
ap.add_argument("--dark-seq", default="MAR007_seq13")
ap.add_argument("--dark-n", type=int, default=5)
ap.add_argument("--dark-npy", default=None, help="dark (switched-state) image as .npy, overrides --dark-seq")
ap.add_argument("--um-per-px", type=float, default=0.536)
ap.add_argument("--adaptive-thr", action="store_true",
                help="per-frame threshold = midpoint of the medians of the switched domain and the rest of the core in the "
                     "previous frame's segmentation (follows slow illumination/contrast drift)")
ap.add_argument("--stop-frac", type=float, default=None,
                help="stop (exclusive) at the first frame where the switched domain covers more than this fraction of the core")
ap.add_argument("--disk-from", default=None,
                help="raw (non-difference) image of the field, e.g. the saturated reference MAR007_ref*.tif: the disk is "
                     "its Otsu-thresholded largest bright component (convex hull) instead of the edge detection on frame 0")
ap.add_argument("--no-flip", action="store_true", help="switched domain already on the left (do not flip left-right)")
ap.add_argument("--fps", type=float, default=1.0)
ap.add_argument("--sigma", type=float, default=2.0)
ap.add_argument("--core-erode", type=int, default=40)
ap.add_argument("--min-area", type=int, default=1000)
ap.add_argument("--margin", type=int, default=5)
ap.add_argument("--median-x", type=int, default=3)
ap.add_argument("--max-fail", type=float, default=0.02, help="max fraction of frames a row may fail (then interpolated)")
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
out = lambda f: os.path.join(a.outdir, f)

ld = lambda f: np.array(Image.open(f)).astype(np.float32)
fs = sorted(glob.glob(os.path.join(a.src, f"{a.seq}_*.tif")))
K = len(fs) if a.last_frame is None else min(a.last_frame, len(fs))
fs = fs[:K]
dark = np.load(a.dark_npy).astype(np.float32) if a.dark_npy else \
    np.mean([ld(f) for f in sorted(glob.glob(os.path.join(a.src, f"{a.dark_seq}_*.tif")))[-a.dark_n:]], 0)

# disk: convex hull of the long border edge segments (4x4 binned mean of the first 5 frames)
f0 = np.mean([ld(f) for f in fs[:5]], 0)
H, W = f0.shape
r4 = f0.reshape(H // 4, 4, W // 4, 4).mean((1, 3))
g = ndi.gaussian_filter(r4, 1.5)
gm = np.hypot(ndi.sobel(g, 0), ndi.sobel(g, 1))
e = gm > np.percentile(gm, 95)
lab, n = ndi.label(e)
sz = ndi.sum(e, lab, range(1, n + 1))
e = np.isin(lab, 1 + np.where(sz >= 150)[0])
disk = np.kron(convex_hull_image(e), np.ones((4, 4), bool))
if a.disk_from:
    rr = ndi.gaussian_filter(ld(a.disk_from).reshape(H // 4, 4, W // 4, 4).mean((1, 3)), 2)
    bb = rr > threshold_otsu(rr)
    lab, n = ndi.label(bb)
    bb = lab == 1 + int(np.argmax(ndi.sum(bb, lab, range(1, n + 1))))
    disk = np.kron(convex_hull_image(bb), np.ones((4, 4), bool))
core = ndi.binary_erosion(disk, iterations=a.core_erode)

sig = lambda k: ndi.gaussian_filter(ld(fs[k]) - dark, a.sigma)
s0 = sig(0)
thr = float(threshold_otsu(s0[core]))


def darkmask(s, t):
    m = (s < t) & core
    return remove_small_holes(remove_small_objects(m, a.min_area), a.min_area)


# flip so that A is on the left and the wall moves to +x
fl = (lambda x: x) if a.no_flip else (lambda x: x[:, ::-1])
coreF = fl(core)
rows_core = np.where(coreF.any(1))[0]
cl = np.full(H, -1)
cr = np.full(H, -1)
cl[rows_core] = [np.argmax(coreF[y]) for y in rows_core]
cr[rows_core] = [W - 1 - np.argmax(coreF[y][::-1]) for y in rows_core]

ta = np.full((H, W), np.inf)
xw = np.full((K, H), -1, int)
xs = np.full((K, H), np.nan)
xe = np.full((K, H), -1, int)
nA = np.zeros(K)
n_other = np.zeros(K, int)
prev = None
thr_all = []
xx = np.arange(W)
for k in range(K):
    s = sig(k)
    if a.adaptive_thr and prev is not None:
        far = ndi.binary_erosion(prev, iterations=15) & core
        farB = ndi.binary_erosion(core & ~prev, iterations=15)
        if far.sum() > 1000 and farB.sum() > 1000:
            thr_k = 0.5 * (float(np.median(s[far])) + float(np.median(s[farB])))
    else:
        thr_k = thr
    thr_all.append(thr_k)
    m = darkmask(s, thr_k)
    lab, n = ndi.label(m)
    if prev is None:
        sz = ndi.sum(m, lab, range(1, n + 1))
        A = lab == 1 + int(np.argmax(sz))
    else:
        ids = np.unique(lab[prev & m]); ids = ids[ids > 0]
        A = np.isin(lab, ids)
    n_other[k] = len(set(range(1, n + 1)) - set(np.unique(lab[A]).tolist()))
    # small bright specks inside or at the core edge of A (noise where the dark-state contrast is weak) -> A
    A = core & ~remove_small_objects(core & ~A, a.min_area)
    prev = A
    if a.stop_frac is not None and k > 0 and A.sum() > a.stop_frac * core.sum():
        print(f"  {a.seq}: domain covers {A.sum() / core.sum():.0%} of the core at frame {k}; stopping", flush=True)
        K = k
        xw, xs, xe, nA, n_other = xw[:K], xs[:K], xe[:K], nA[:K], n_other[:K]
        ta[ta >= K] = np.inf
        break
    AF, sF = fl(A), fl(s)
    newly = AF & np.isinf(ta)
    ta[newly] = -1 if k == 0 else k
    nA[k] = AF.sum()
    for y in rows_core:
        row = AF[y, cl[y]:cr[y] + 1]
        if not row.any():
            continue
        na = int(row.sum())
        xe[k, y] = cl[y] + na - 1
        gap = np.argmin(row) if not row.all() else row.size
        if gap == 0:                         # A does not touch the left core edge in this row
            continue
        x = cl[y] + gap - 1
        xw[k, y] = x
        if x + 1 <= cr[y]:
            s1, s2 = sF[y, x], sF[y, x + 1]
            xs[k, y] = x + np.clip((thr_k - s1) / (s2 - s1), 0, 1) if s2 != s1 else x + 0.5
    if k % 50 == 0:
        print(f"  {a.seq}: frame {k}/{K}  A = {nA[k] * a.um_per_px ** 2 / 1e3:.0f}e3 um^2, other dark comps {n_other[k]}")

# rows: front inside the core (margin); a row is kept if it fails in <= --max-fail of the frames, and the
# failed frames are interpolated in time (single-frame noise at low contrast)
fail = np.ones((K, H), bool)
for y in rows_core:
    fail[:, y] = (xw[:, y] < cl[y] + a.margin) | (xw[:, y] > cr[y] - a.margin) | ~np.isfinite(xs[:, y])
frac_fail = fail.mean(0)
ok = frac_fail <= a.max_fail
ok[[0, -1]] = False
# a failure run that lasts to the last frame cannot be interpolated (the front has left the core in that row, e.g. the
# wall reached the far diaphragm edge at the top): stop the movie before the earliest such run (if >= half the frames
# remain), otherwise drop those rows
trail = [(y, K - int(np.argmin(fail[::-1, y]))) for y in np.nonzero(ok & fail[-1])[0]]
if trail:
    kcut = min(k0 for _, k0 in trail)
    if kcut >= K // 2:
        print(f"  {a.seq}: front leaves the core at frame {kcut} in {len(trail)} rows; stopping there", flush=True)
        K = kcut
        xw, xs, xe, nA, n_other, fail = xw[:K], xs[:K], xe[:K], nA[:K], n_other[:K], fail[:K]
        ta[ta >= K] = np.inf
        thr_all = thr_all[:K]
        frac_fail = fail.mean(0)
        ok = frac_fail <= a.max_fail
        ok[[0, -1]] = False
    else:
        ok[[y for y, _ in trail]] = False
tt = np.arange(K)
for y in np.nonzero(ok & fail.any(0))[0]:
    g = ~fail[:, y]
    for z in (xw, xe):
        z[~g, y] = np.round(np.interp(tt[~g], tt[g], z[g, y])).astype(int)
    xs[~g, y] = np.interp(tt[~g], tt[g], xs[g, y])
idx = np.nonzero(ok)[0]
runs = np.split(idx, np.nonzero(np.diff(idx) != 1)[0] + 1)
run = max(runs, key=len)
y0, y1 = int(run[0]), int(run[-1]) + 1
print(f"rows used {y0}..{y1 - 1} ({y1 - y0} rows, {(y1 - y0) * a.um_per_px:.0f} um); interpolated (row, frame) pairs in them: "
      f"{fail[:, y0:y1].mean():.2%}")

sl = slice(y0, y1)
Xw, Xs, Xe = xw[:, sl], xs[:, sl], xe[:, sl]
nover = (Xw != Xe).sum(1)
mf = lambda z: ndi.median_filter(z, size=(1, a.median_x), mode="nearest") if a.median_x > 1 else z
h = mf((W - 1 - Xe).astype(float))
h_sub = mf(W - 1 - (Xe + (Xs - Xw)))
yy = np.arange(y0, y1)
theta = np.degrees(np.arctan([np.polyfit(yy, r, 1)[0] for r in Xs]))

# arrival map in the rows used: the A side outside the core is "switched before frame 0"
taC = ta[sl].copy()
for i, y in enumerate(yy):
    taC[i, :cl[y]] = -1
    taC[i, cr[y] + 1:] = np.inf
taC = ndi.median_filter(np.where(np.isinf(taC), 1e9, taC), 5)
taC[taC >= 1e9] = np.inf

np.save(out("h_xt.npy"), h)
np.save(out("h_xt_sub.npy"), h_sub)
np.savez_compressed(out("tiff_extra.npz"), arrival=taC, defects=np.zeros_like(taC, bool), x_wall=Xw, x_wall_sub=Xs,
                    x_eff=Xe, theta_deg=theta, n_overhang_rows=nover, rows=yy, drift=np.zeros((K, 2)),
                    core=coreF[sl], threshold=thr, interpolated=fail[:, sl], frac_fail_all_rows=frac_fail, area_A_px=nA, n_other_dark=n_other)
H_uT = a.current * a.coil_mT_per_A * 1e3 if a.current is not None else None
meta = dict(fps=a.fps, um_per_px=a.um_per_px, n_frames_video=int(K), frame_first=0, frame_last=int(K - 1),
            height=int(W), width=int(y1 - y0), source="tiff", tiff_dir=os.path.abspath(a.src), seq=a.seq,
            sample="CoFeB MAR007, t = 0.7 nm, room temperature", current_A=a.current, field_uT=H_uT,
            coil_mT_per_A=a.coil_mT_per_A, dark_reference=(a.dark_npy if a.dark_npy else f"{a.dark_seq}, last {a.dark_n} frames"),
            threshold=thr, adaptive_threshold=bool(a.adaptive_thr),
            threshold_range=[float(min(thr_all)), float(max(thr_all))], sigma=a.sigma, core_erode=a.core_erode, min_area=a.min_area, margin=a.margin,
            median_x=a.median_x, rows_used=[y0, y1], max_fail=a.max_fail,
            frac_interpolated=float(fail[:, sl].mean()),
            tilt_deg_first_last=[float(theta[0]), float(theta[-1])], overhang_rows_mean=float(nover.mean()),
            frac_rows_front_ne_area=float((Xw != Xe).mean()),
            geometry=("image not flipped (--no-flip)" if a.no_flip else "image flipped left-right") + ": switched domain (central wall) on the left, wall moving to +x; "
                     "h = W-1-x_eff (switched area per row), columns of h = image rows y0..y1-1; "
                     "domains nucleated at the disk border are ignored",
            warnings=[])
hm = h_sub.mean(1)
if (np.diff(hm) > 0).any():
    meta["warnings"].append(f"the mean <h> moves backward in {(np.diff(hm) > 0).sum()} steps")
json.dump(meta, open(out("meta.json"), "w"), indent=2)

# QC in the original orientation
idx = np.unique(np.linspace(0, K - 1, 6).round().astype(int))
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
for ax, k in zip(axes.flat, idx):
    im = ld(fs[k])
    ax.imshow(im, cmap="gray", vmin=np.percentile(im, 1), vmax=np.percentile(im, 99))
    ax.plot((Xe + Xs - Xw)[k] if a.no_flip else W - 1 - (Xe + Xs - Xw)[k], yy, "r-", lw=0.7)
    ax.contour(core, [0.5], colors="y", linewidths=0.5, linestyles=":")
    ax.set_title(f"{a.seq} frame {k} (t = {k / a.fps:.0f} s)", fontsize=9)
    ax.set_xticks([]); ax.set_yticks([])
fig.tight_layout(); fig.savefig(out("fig_qc_overlay.png"), dpi=70); plt.close(fig)

fig, ax = plt.subplots(figsize=(9, 5))
im = ax.imshow((Xs if a.no_flip else W - 1 - Xs).T, aspect="auto", cmap="magma", extent=(-0.5, K - 0.5, y1, y0))
plt.colorbar(im, ax=ax, label="wall x, original image [px]")
ax.set_xlabel("frame"); ax.set_ylabel("y [px]"); ax.set_title(f"{a.seq}: kymograph")
fig.tight_layout(); fig.savefig(out("fig_kymograph.png"), dpi=100); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(15, 3.6))
ax[0].plot(theta); ax[0].set_ylabel("tilt [deg]")
ax[1].plot(nover); ax[1].set_ylabel(f"overhang rows (of {y1 - y0})")
ax[2].plot((hm[0] - hm) * a.um_per_px); ax[2].set_ylabel("mean advance [um]")
for x in ax:
    x.set_xlabel("frame")
fig.tight_layout(); fig.savefig(out("fig_tilt_overhang.png"), dpi=100); plt.close(fig)
v = np.polyfit(np.arange(K) / a.fps, (hm[0] - hm) * a.um_per_px, 1)[0]
print(f"done: {K} frames, rows {y0}..{y1 - 1}, thr {thr:.1f}, mean v = {v:.3f} um/s, "
      f"tilt {theta[0]:.1f} -> {theta[-1]:.1f} deg, overhang rows {nover.mean():.1f}")
