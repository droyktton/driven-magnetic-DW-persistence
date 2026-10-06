"""Step 1: frame extraction and detection of the domain wall position h(x,t).

Usage: python step1_extract_wall.py VIDEO [--outdir DIR] [--um-per-px 0.17] [--fps F]
                                   [--roi X0 X1 Y0 Y1] [--rotate K] [--flipud] [--fliplr]
                                   [--invert] [--sigma-y 2] [--median-x 5] [--edge-margin 2]
       python step1_extract_wall.py TIFF_DIR [--outdir DIR] [--um-per-px X] [--fps F]
                                   [--crop-right auto|N] [--nref 5] [--trim-y 12] [--qc-every 50]

VIDEO: frames are extracted with ffmpeg and h is the first dark row from the top (global Otsu).

TIFF_DIR: a folder of TIFF frames (e.g. Micro-Manager Pos0/), dark domain on the LEFT growing
to the RIGHT (wall roughly vertical, may be tilted). Steps:
  1. rigid drift of every frame vs frame 0 (phase correlation of the high-passed images),
     smoothed in time, and registration of every frame;
  2. crop at a fixed column on the right where the field of view ends (auto or N), and trim
     the rows at the top/bottom that the vertical drift leaves without data;
  3. per-pixel bright (first --nref frames) and dark (last --nref frames) references, after
     normalising every frame by the median of the pixels that never switch. This cancels
     vignetting and static defects;
  4. arrival-time map: the frame at which each pixel switches (number of frames in which it is
     still bright; the advance is assumed monotonic). Static defects inside the swept area,
     which never show contrast, get the arrival time of the nearest valid pixel;
  5. domain at frame t = pixels with arrival <= t connected to the left edge; the drawn wall is
     the right-most domain pixel of each row after filling holes (front of overhangs), subpixel
     by linear interpolation of the 0.5 crossing of the normalised intensity;
  6. position used for the analysis: x_eff(y,t) = switched pixels in row y - 1 (+ the subpixel
     fraction of the front). It equals the front where the wall is single-valued, and changes
     only by the area that really switched where the wall folds around defects.
  Then h = W - 1 - x_eff (equivalent to --rotate 1) so that h decreases with t as in the video
  convention, and the columns of h are the image rows y. fps comes from the Micro-Manager
  metadata.txt timestamps when present.

Convention (after roi/rotate/flip/invert): bright domain on top, dark domain at the
bottom, wall moving upward (h decreases with t). Order of the transformations:
roi (in original video coordinates) -> rotate (K x 90° counter-clockwise) -> flipud ->
fliplr -> invert.

Outputs in DIR (default: the video file name without extension):
  frames/            extracted frames (ffmpeg -vsync 0)
  h_xt.npy           integer h(x,t), (frames used x width), in px
  h_xt_sub.npy       subpixel h(x,t)
  meta.json          video and detection parameters (read by step2 and step3)
  fig_qc_overlay.png, fig_h_mean.png
  TIFF only: tiff_extra.npz (drift, arrival map, tilt, overhangs), fig_drift_arrival.png,
             fig_kymograph.png, fig_tilt_overhang.png, qc/overlay.mp4 and qc/overlay_NNNN.png
"""
import argparse
import glob
import json
import os
import re
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.image import imread
from scipy.ndimage import (binary_dilation, binary_fill_holes, distance_transform_edt,
                           gaussian_filter, gaussian_filter1d, label, median_filter)
from scipy.ndimage import shift as nd_shift
from skimage.filters import threshold_otsu

ap = argparse.ArgumentParser()
ap.add_argument("video", help="video file, or a folder of TIFF frames")
ap.add_argument("--outdir", default=None,
                help="output folder (default: video file name; for TIFF the folder name, or its "
                     "parent if the folder is called PosN)")
ap.add_argument("--um-per-px", type=float, default=0.17,
                help="spatial scale in µm/px, stored in meta.json (0 = unknown)")
ap.add_argument("--fps", type=float, default=None, help="actual frame rate (default: the one stored in the video)")
ap.add_argument("--roi", type=int, nargs=4, metavar=("X0", "X1", "Y0", "Y1"), default=None,
                help="crop in original video px: columns X0:X1, rows Y0:Y1")
ap.add_argument("--rotate", type=int, default=0, help="rotate K x 90° counter-clockwise")
ap.add_argument("--flipud", action="store_true", help="flip top/bottom")
ap.add_argument("--fliplr", action="store_true", help="flip left/right")
ap.add_argument("--invert", action="store_true", help="invert contrast (bright <-> dark)")
ap.add_argument("--sigma-y", type=float, default=2.0, help="vertical smoothing (px) before thresholding")
ap.add_argument("--median-x", type=int, default=5, help="median filter kernel along x")
ap.add_argument("--edge-margin", type=int, default=2,
                help="a frame is valid only if, in every column, the wall is more than this "
                     "many px away from the top and bottom edges")
ap.add_argument("--retreat-px", type=float, default=1.0,
                help="a Delta h(tau=1) larger than this counts as a real retreat")
tg = ap.add_argument_group("TIFF folder input (--roi/--rotate/--flipud/--fliplr/--invert/--sigma-y are ignored)")
tg.add_argument("--crop-right", default="auto",
                help="keep columns x < N (camera coordinates); auto = where the illumination of the "
                     "field of view drops below 90%% of its plateau, minus 20 px")
tg.add_argument("--nref", type=int, default=5, help="frames averaged for the bright/dark references")
tg.add_argument("--trim-y", type=int, default=12,
                help="rows discarded at the top and bottom, on top of the maximum vertical drift")
tg.add_argument("--qc-every", type=int, default=50, help="write qc/overlay_NNNN.png every N frames")
tg.add_argument("--no-movie", action="store_true", help="do not write qc/overlay.mp4")
tg.add_argument("--workers", type=int, default=4, help="threads for reading/processing frames")
args = ap.parse_args()

is_tiff = os.path.isdir(args.video)
if args.outdir:
    outdir = args.outdir
elif is_tiff:
    d = os.path.normpath(os.path.abspath(args.video))
    outdir = os.path.basename(os.path.dirname(d)) if re.fullmatch(r"Pos\d+", os.path.basename(d)) \
        else os.path.basename(d)
else:
    outdir = os.path.splitext(os.path.basename(args.video))[0]
frames_dir = os.path.join(outdir, "frames")


def out(name):
    return os.path.join(outdir, name)


def video_fps(path):
    import imageio_ffmpeg
    gen = imageio_ffmpeg.read_frames(path)
    meta = next(gen)
    gen.close()
    return float(meta["fps"])


def ffmpeg_exe():
    try:
        import imageio_ffmpeg
        return imageio_ffmpeg.get_ffmpeg_exe()
    except ImportError:
        exe = shutil.which("ffmpeg")
        if exe is None:
            raise SystemExit("ERROR: neither imageio-ffmpeg nor an ffmpeg binary is available")
        return exe


def extract_frames():
    os.makedirs(frames_dir, exist_ok=True)
    for f in glob.glob(os.path.join(frames_dir, "f_*.png")):
        os.remove(f)
    subprocess.run([ffmpeg_exe(), "-v", "error", "-i", args.video,
                    "-vsync", "0", os.path.join(frames_dir, "f_%04d.png")], check=True)


def load_stack():
    files = sorted(glob.glob(os.path.join(frames_dir, "f_*.png")))
    frames = []
    for f in files:
        im = imread(f).astype(float)
        if im.ndim == 3:
            im = im[..., :3].mean(-1)
        frames.append(im)
    stack = np.array(frames)          # (T, H, W)
    if args.roi:
        x0, x1, y0, y1 = args.roi
        stack = stack[:, y0:y1, x0:x1]
    if args.rotate % 4:
        stack = np.rot90(stack, args.rotate, axes=(1, 2))
    if args.flipud:
        stack = stack[:, ::-1, :]
    if args.fliplr:
        stack = stack[:, :, ::-1]
    if args.invert:
        stack = stack.max() + stack.min() - stack
    return np.ascontiguousarray(stack)


def detect_wall(stack, thr):
    smooth = gaussian_filter1d(stack, args.sigma_y, axis=1)
    below = smooth < thr
    H = stack.shape[1]
    # first row (from the top) below the threshold; H if the whole column is bright
    r = np.where(below.any(axis=1), below.argmax(axis=1), H)
    h = r.astype(float)
    # subpixel: linear interpolation of the threshold crossing between rows r-1 and r
    ok = (r > 0) & (r < H)
    t_idx, x_idx = np.nonzero(ok)
    rr = r[ok]
    s0 = smooth[t_idx, rr - 1, x_idx]
    s1 = smooth[t_idx, rr, x_idx]
    h_sub = h.copy()
    h_sub[ok] = rr - 1 + (s0 - thr) / (s0 - s1)
    h = median_filter(h, size=(1, args.median_x), mode="nearest")
    h_sub = median_filter(h_sub, size=(1, args.median_x), mode="nearest")
    return r, h, h_sub


def ranges(idx):
    """[0,1,2,5,7,8] -> '0-2, 5, 7-8'"""
    out, start = [], None
    for i, v in enumerate(idx):
        if start is None:
            start = v
        if i + 1 == len(idx) or idx[i + 1] != v + 1:
            out.append(f"{start}" if start == v else f"{start}-{v}")
            start = None
    return ", ".join(out)


def longest_run(mask):
    """(start, end) of the longest contiguous run with mask == True (end exclusive)."""
    best, start = (0, 0), None
    for i, m in enumerate(list(mask) + [False]):
        if m and start is None:
            start = i
        elif not m and start is not None:
            if i - start > best[1] - best[0]:
                best = (start, i)
            start = None
    return best


# ------------------------------------------------------------------ TIFF folder input
def pmap(fn, items):
    """Ordered parallel map in chunks (bounded memory)."""
    items = list(items)
    with ThreadPoolExecutor(args.workers) as ex:
        for c in range(0, len(items), 4 * args.workers):
            yield from ex.map(fn, items[c:c + 4 * args.workers])


def tiff_files(d):
    files = [f for f in glob.glob(os.path.join(d, "*")) if re.search(r"\.tiff?$", f, re.I)]
    num = lambda f: [int(n) for n in re.findall(r"\d+", os.path.basename(f))]
    return sorted(files, key=num)


def mm_fps(d):
    """Frame rate from the Micro-Manager metadata.txt timestamps (None if absent)."""
    try:
        with open(os.path.join(d, "metadata.txt")) as f:
            m = json.load(f)
    except (OSError, ValueError):
        return None
    keys = sorted((k for k in m if k.startswith("FrameKey")), key=lambda k: int(k.split("-")[1]))
    t = np.array([m[k]["ElapsedTime-ms"] for k in keys], float)
    return 1e3 / np.median(np.diff(t)) if len(t) > 1 else None


def read_tiff(f):
    import tifffile
    im = tifffile.imread(f).astype(np.float32)
    return im[..., :3].mean(-1) if im.ndim == 3 else im


def highpass(a):
    return gaussian_filter(a, 1.5) - gaussian_filter(a, 15)


def auto_crop_right(files):
    """First column (right half) where the illumination falls below 90% of its plateau, - 20 px."""
    im = np.mean([read_tiff(files[i]) for i in (0, len(files) // 2, -1)], axis=0)
    H, W = im.shape
    prof = np.median(im[H // 4: 3 * H // 4], axis=0)
    plateau = np.median(prof[W // 4: W // 2])
    low = np.where(prof[W // 2:] < 0.9 * plateau)[0]
    return W if low.size == 0 else int(W // 2 + low[0] - 20)


def wall_from_tiffs():
    from skimage.registration import phase_cross_correlation
    files = tiff_files(args.video)
    if len(files) < 10:
        raise SystemExit(f"ERROR: only {len(files)} TIFF files in {args.video}")
    fps = args.fps or mm_fps(args.video)
    if not fps:
        raise SystemExit("ERROR: no metadata.txt timestamps; pass --fps")
    T = len(files)
    H0, W0 = read_tiff(files[0]).shape
    xc = auto_crop_right(files) if args.crop_right == "auto" else int(args.crop_right)
    print(f"tiff folder={args.video}  frames={T}  size={W0}x{H0}  fps={fps:g} "
          f"(dt={1 / fps:g} s)  crop right at x={xc}")

    # 1. drift vs frame 0, smoothed in time (it is a slow mechanical drift)
    ref = highpass(read_tiff(files[0]))[:, :xc]
    drift = np.array(list(pmap(lambda f: phase_cross_correlation(
        ref, highpass(read_tiff(f))[:, :xc], upsample_factor=10, normalization=None)[0], files)))
    drift_s = gaussian_filter1d(drift, 3, axis=0, mode="nearest")
    print(f"drift: dy {drift[:, 0].min():+.1f}..{drift[:, 0].max():+.1f} px, "
          f"dx {drift[:, 1].min():+.1f}..{drift[:, 1].max():+.1f} px; "
          f"scatter around smoothed curve {np.std(drift - drift_s, 0).round(3)} px")
    # rows valid in every registered frame (shift by +dy moves the content down)
    y0 = int(np.ceil(max(drift_s[:, 0].max(), 0))) + args.trim_y
    y1 = H0 + int(np.floor(min(drift_s[:, 0].min(), 0))) - args.trim_y

    def reg(k):
        g = nd_shift(gaussian_filter(read_tiff(files[k]), 1.5), drift_s[k], order=1, mode="nearest")
        return g[y0:y1, :xc]

    # 2. normalisation region: pixels that never switch (first pass with full-frame medians)
    first = [reg(k) for k in range(args.nref)]
    last = [reg(k) for k in range(T - args.nref, T)]
    B0 = np.mean([g / np.median(g) for g in first], 0)
    D0 = np.mean([g / np.median(g) for g in last], 0)
    c0 = (B0 - D0) / B0
    R = ~binary_dilation(c0 > threshold_otsu(c0), iterations=10) & (B0 > 0.6 * np.median(B0))
    norm = lambda g: g / np.median(g[R])
    B = np.mean([norm(g) for g in first], 0)
    Dk = np.mean([norm(g) for g in last], 0)
    del first, last
    c = (B - Dk) / B
    thr_c = threshold_otsu(c)
    M = c > thr_c                                     # pixels swept during the movie
    den = np.where(M, B - Dk, 1.0)

    # 3. arrival-time map: number of frames in which the pixel is still bright
    ta = np.zeros(M.shape, np.int32)
    for bright in pmap(lambda k: (norm(reg(k)) - Dk) / den > 0.5, range(T)):
        ta += bright & M
    lab, _ = label(~M)
    edge_ids = lambda col: np.setdiff1d(np.unique(col), [0])
    A = np.isin(lab, edge_ids(lab[:, 0]))                       # initial domain (always dark)
    N = np.isin(lab, edge_ids(lab[:, -1])) & ~A                 # never reached
    holes = ~M & ~A & ~N                                        # static defects inside the swept area
    # "switching" specks outside the main swept band are defect halos, not wall motion:
    # they belong to the initial domain if connected to the left edge, otherwise to the never-reached side
    lab, _ = label(M | holes)
    main = lab == np.argmax(np.bincount(lab.ravel())[1:]) + 1
    lab, _ = label(~main)
    left = np.isin(lab, edge_ids(lab[:, 0]))
    holes &= main
    ta = ta.astype(float)
    ta[~main & left], ta[~main & ~left] = -1, np.inf
    if holes.any():
        _, (iy, ix) = distance_transform_edt(holes, return_indices=True)
        ta[holes] = ta[iy[holes], ix[holes]]
    ta = median_filter(ta, 5)

    # 4. wall per frame + QC renders
    Hc, Wc = M.shape
    lo, hi = np.percentile(norm(reg(0)), [1, 99.5])
    os.makedirs(out("qc"), exist_ok=True)
    for f in glob.glob(os.path.join(out("qc"), "overlay_*.png")):
        os.remove(f)

    def wall(k):
        dom = ta <= k
        lab2, _ = label(dom)
        dom = np.isin(lab2, edge_ids(lab2[:, 0]))
        x_eff = dom.sum(1) - 1                  # switched area per row (see docstring)
        dom = binary_fill_holes(dom)
        rows = dom.any(1)
        xr = np.where(rows, Wc - 1 - np.argmax(dom[:, ::-1], axis=1), -1)
        nover = int(((np.diff(dom.astype(np.int8), axis=1) == -1).sum(1) > 1).sum())
        g = norm(reg(k))
        s = (g - Dk) / den
        xs = np.full(Hc, np.nan)
        ok = rows & (xr + 1 < Wc)
        y = np.nonzero(ok)[0]
        s0, s1 = s[y, xr[ok]], s[y, xr[ok] + 1]
        with np.errstate(invalid="ignore", divide="ignore"):
            frac = np.clip((0.5 - s0) / (s1 - s0), 0, 1)
        xs[ok] = xr[ok] + np.nan_to_num(frac, nan=0.5)
        img = np.clip((g - lo) / (hi - lo) * 255, 0, 255).astype(np.uint8)
        return xr, xs, x_eff, nover, img

    movie = None
    if not args.no_movie:
        hp, wp = Hc + Hc % 2, Wc + Wc % 2
        movie = subprocess.Popen([ffmpeg_exe(), "-v", "error", "-y", "-f", "rawvideo", "-pix_fmt", "rgb24",
                                  "-s", f"{wp}x{hp}", "-r", "25", "-i", "-", "-c:v", "libx264",
                                  "-pix_fmt", "yuv420p", "-crf", "20", out("qc/overlay.mp4")],
                                 stdin=subprocess.PIPE)
    from PIL import Image, ImageDraw
    xw = np.empty((T, Hc), int); xw_sub = np.empty((T, Hc)); xe = np.empty((T, Hc), int)
    nover = np.empty(T, int)
    yy = np.arange(Hc)
    for k, (xr, xs, x_eff, no, img) in enumerate(pmap(wall, range(T))):
        xw[k], xw_sub[k], xe[k], nover[k] = xr, xs, x_eff, no
        rgb = np.repeat(img[..., None], 3, axis=2)
        ok = np.isfinite(xs)
        for dx in (0, 1):                                   # 2 px wide red line
            xi = np.clip(np.round(xs[ok]).astype(int) + dx, 0, Wc - 1)
            rgb[yy[ok], xi] = (255, 0, 0)
        im = Image.fromarray(rgb)
        ImageDraw.Draw(im).text((10, 10), f"frame {k}  t = {k / fps / 3600:.2f} h", fill=(255, 255, 0))
        if k % args.qc_every == 0 or k == T - 1:
            im.save(out(f"qc/overlay_{k:04d}.png"))
        if movie:
            fr = np.zeros((hp, wp, 3), np.uint8); fr[:Hc, :Wc] = np.asarray(im)
            movie.stdin.write(fr.tobytes())
        if k % 200 == 0:
            print(f"  wall: frame {k}/{T}")
    if movie:
        movie.stdin.close(); movie.wait()

    theta = np.degrees(np.arctan([np.polyfit(yy[np.isfinite(r)], r[np.isfinite(r)], 1)[0]
                                  if np.isfinite(r).sum() > 2 else np.nan for r in xw_sub]))
    np.savez_compressed(out("tiff_extra.npz"), drift=drift, drift_smooth=drift_s, arrival=ta, swept=M,
             defects=holes, contrast=c, x_wall=xw, x_wall_sub=xw_sub, x_eff=xe, theta_deg=theta,
             n_overhang_rows=nover, rows=np.arange(y0, y1), crop_right=xc)

    # QC figures
    fig, ax = plt.subplots(1, 3, figsize=(18, 5))
    ax[0].plot(drift[:, 0], ".", ms=2, label="dy"); ax[0].plot(drift[:, 1], ".", ms=2, label="dx")
    ax[0].plot(drift_s, "k-", lw=0.8); ax[0].legend()
    ax[0].set_xlabel("frame"); ax[0].set_ylabel("drift vs frame 0 [px]"); ax[0].set_title("Sample drift")
    im = ax[1].imshow(np.where(np.isfinite(ta) & (ta >= 0), ta, np.nan), cmap="viridis",
                      extent=(0, Wc, y1, y0))
    plt.colorbar(im, ax=ax[1], label="arrival frame")
    ax[1].contour(np.arange(Wc), np.arange(y0, y1), holes, [0.5], colors="r", linewidths=0.5)
    ax[1].set_title("Arrival-time map (red: defects, inpainted)")
    ax[2].imshow(c, cmap="gray", vmin=-0.5 * thr_c, vmax=3 * thr_c, extent=(0, Wc, y1, y0))
    ax[2].set_title(f"Switching contrast (B-D)/B, threshold {thr_c:.3f}")
    fig.tight_layout(); fig.savefig(out("fig_drift_arrival.png"), dpi=110); plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    im = ax.imshow(xw_sub.T, aspect="auto", cmap="magma", extent=(-0.5, T - 0.5, y1, y0))
    plt.colorbar(im, ax=ax, label="wall x [px]")
    ax.set_xlabel("frame"); ax.set_ylabel("y [px]"); ax.set_title("Kymograph x(y,t)")
    fig.tight_layout(); fig.savefig(out("fig_kymograph.png"), dpi=110); plt.close(fig)

    fig, ax = plt.subplots(1, 2, figsize=(11, 3.6))
    ax[0].plot(theta); ax[0].set_ylabel(r"tilt $\theta$ [deg]  (fit x = a + y tan$\theta$)")
    ax[1].plot(nover); ax[1].set_ylabel(f"overhang rows (of {Hc})")
    for a in ax:
        a.set_xlabel("frame")
    fig.tight_layout(); fig.savefig(out("fig_tilt_overhang.png"), dpi=110); plt.close(fig)

    # position used for the analysis: the switched area per row, x_eff. It equals the front x where
    # the wall is single-valued; where it folds around defects, the front can leap when a thin
    # unswitched channel is pinched off, while x_eff only changes by the area that really switched.
    # Subpixel: x_eff + the fractional part of the front crossing.
    # step1 convention: h = W - 1 - x decreases as the wall advances; columns of h = rows y
    r = np.where(xw >= 0, Wc - 1 - xe, Wc)               # Wc = wall missing in that row
    h = median_filter(r.astype(float), size=(1, args.median_x), mode="nearest")
    h_sub = median_filter(Wc - 1 - (xe + (xw_sub - xw)), size=(1, args.median_x), mode="nearest")
    extra = dict(source="tiff", tiff_dir=os.path.abspath(args.video), crop_right=xc,
                 rows_used=[y0, y1], drift_max_px=np.abs(drift).max(0).round(2).tolist(),
                 contrast_threshold=float(thr_c), n_defect_px=int(holes.sum()),
                 frac_rows_front_ne_area=float((xw != xe).mean()),
                 tilt_deg_first_last=[float(theta[0]), float(theta[-1])],
                 overhang_rows_mean=float(nover.mean()), nref=args.nref,
                 geometry="dark domain on the left, wall moving to +x; h = W-1-x_eff (switched "
                          "area per row), "
                          "columns of h = image rows y0..y1-1 (registered to frame 0)")

    def qc_overlay(t0, t1, h_sub):
        idx = np.unique(np.linspace(t0, t1 - 1, 6).round().astype(int))
        fig, axes = plt.subplots(2, 3, figsize=(18, 11))
        for ax, i in zip(axes.flat, idx):
            ax.imshow(norm(reg(i)), cmap="gray", vmin=lo, vmax=hi, extent=(0, Wc, y1, y0))
            ax.plot(xw_sub[i], np.arange(y0, y1), "r-", lw=0.8)
            ax.set_title(f"frame {i} (t = {i / fps / 3600:.2f} h), registered", fontsize=9)
        fig.tight_layout(); fig.savefig(out("fig_qc_overlay.png"), dpi=100); plt.close(fig)

    return dict(fps=fps, T=T, H=Wc, W=Hc, r=r, h=h, h_sub=h_sub, extra=extra, qc_overlay=qc_overlay,
                pos_label=r"mean wall position $\langle x \rangle_y$ + 1 [px]")


def wall_from_video():
    fps = args.fps or video_fps(args.video)
    extract_frames()
    stack = load_stack()
    T, H, W = stack.shape
    thr = threshold_otsu(stack)
    r, h, h_sub = detect_wall(stack, thr)
    print(f"video={args.video}  frames={T}  size (after roi/rotation)={W}x{H}  "
          f"fps={fps:g}  Otsu global={thr:.4f}")
    extra = dict(video=os.path.abspath(args.video), otsu_threshold=float(thr), roi=args.roi,
                 rotate=args.rotate, flipud=args.flipud, fliplr=args.fliplr, invert=args.invert,
                 sigma_y=args.sigma_y)

    def qc_overlay(t0, t1, h_sub):
        # QC: h(x,t) overlaid on representative frames of the range used
        idx = np.unique(np.linspace(t0, t1 - 1, 6).round().astype(int))
        fig, axes = plt.subplots(len(idx), 1, figsize=(8, 2.1 * len(idx) * max(H / W, 0.3) / 0.48))
        for ax, i in zip(np.atleast_1d(axes), idx):
            ax.imshow(stack[i], cmap="gray")
            ax.plot(np.arange(W), h_sub[i - t0], "r-", lw=0.8)
            ax.set_title(f"frame {i}", fontsize=9)
            ax.set_xticks([]); ax.set_yticks([])
        fig.tight_layout()
        fig.savefig(out("fig_qc_overlay.png"), dpi=150)
        plt.close(fig)

    return dict(fps=fps, T=T, H=H, W=W, r=r, h=h, h_sub=h_sub, extra=extra, qc_overlay=qc_overlay,
                pos_label=r"$H - \langle h \rangle_x$  [px]")


def main():
    os.makedirs(outdir, exist_ok=True)
    res = wall_from_tiffs() if is_tiff else wall_from_video()
    fps, T, H, W, r, h, h_sub = (res[k] for k in ("fps", "T", "H", "W", "r", "h", "h_sub"))

    # --- edge check: the wall must be inside the frame in every column. If it touches the
    # top edge (or is missing), h gets stuck and those columns would look "still",
    # inflating Pi and chi4. The longest contiguous run of valid frames is kept.
    m = args.edge_margin
    frame_ok = ((r >= m) & (r < H - m)).all(axis=1)
    t0, t1 = longest_run(frame_ok)
    warnings = []
    if t1 - t0 < T:
        bad = np.where(~frame_ok)[0]
        warnings.append(f"{T - (t1 - t0)} frames dropped: the wall touches the edge or is missing "
                        f"in frames {ranges(bad.tolist())}; using frames {t0}..{t1 - 1}")
    if t1 - t0 < 10:
        raise SystemExit(f"ERROR: only {t1 - t0} consecutive valid frames; check the "
                         f"orientation (--rotate/--flipud/--invert) or --roi.")
    h, h_sub = h[t0:t1], h_sub[t0:t1]

    # --- monotonic advance check (assumed by the step2 noise estimate)
    hm = h_sub.mean(axis=1)
    v = -np.diff(hm)
    dh1 = h_sub[1:] - h_sub[:-1]
    frac_retreat = float((dh1 > args.retreat_px).mean())
    if hm[0] - hm[-1] <= 0:
        warnings.append("on average the wall does NOT move upward: check the orientation "
                        "(--flipud/--rotate) or the contrast (--invert)")
    if (v < 0).any():
        warnings.append(f"the mean <h> moves backward in {(v < 0).sum()} steps")
    if frac_retreat > 0.01:
        warnings.append(f"{100 * frac_retreat:.1f}% of Delta h(tau=1) are retreats "
                        f"> {args.retreat_px:g} px: the automatic step2 eps (3 sigma_Dh) "
                        f"is unreliable; pass --eps explicitly")

    np.save(out("h_xt.npy"), h)
    np.save(out("h_xt_sub.npy"), h_sub)
    meta = dict(fps=fps, um_per_px=args.um_per_px or None, n_frames_video=int(T),
                frame_first=int(t0), frame_last=int(t1 - 1), height=int(H), width=int(W),
                **res["extra"], median_x=args.median_x, edge_margin=m,
                frac_retreat=frac_retreat, warnings=warnings)
    with open(out("meta.json"), "w") as f:
        json.dump(meta, f, indent=2, ensure_ascii=False)

    res["qc_overlay"](t0, t1, h_sub)

    # <h>_x vs frame (position measured from the bottom so that it increases as the wall advances)
    frames = np.arange(t0, t1)
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    axes[0].plot(frames, H - hm, "o-", ms=3)
    axes[0].set_xlabel("frame"); axes[0].set_ylabel(res["pos_label"])
    axes[0].set_title("Mean wall position")
    axes[1].plot(frames[1:], v, "o-", ms=3)
    axes[1].axhline(0, color="k", lw=0.5)
    axes[1].set_xlabel("frame"); axes[1].set_ylabel(r"$-\Delta\langle h \rangle$ [px/frame]")
    axes[1].set_title("Mean velocity")
    for ax in axes:
        ax.set_xlim(-0.5, T - 0.5)
        for a, b in ((0, t0), (t1, T)):   # dropped frames
            if b > a:
                ax.axvspan(a - 0.5, b - 0.5, color="0.85", zorder=0)
    fig.tight_layout()
    fig.savefig(out("fig_h_mean.png"), dpi=150)
    plt.close(fig)

    print(f"frames used: {t0}..{t1 - 1} ({t1 - t0} of {T})")
    print(f"total advance = {hm[0] - hm[-1]:.1f} px; mean v = {v.mean():.2f} px/frame; "
          f"retreats > {args.retreat_px:g} px: {100 * frac_retreat:.2f}%")
    for w in warnings:
        print("WARNING:", w)
    print(f"outputs in {outdir}/")


if __name__ == "__main__":
    main()
