"""Step 7: roughness and structure factor of the wall, on its base plane.

Usage: python step7_roughness.py DIR [--every 5] [--rmin 16] [--mask-d 20]
       python step7_roughness.py DIR --from-h [--every 1]       (video studies, simulations)
       python step7_roughness.py --selftest                     (synthetic wall with known zeta)
Outputs in DIR: roughness.npz, fig_roughness.png, contours.npz (TIFF studies), local_width.npz, fig_local_width.png,
local_sq.npz, fig_local_sq.png, height_distribution.npz, fig_height_distribution.png and a summary on stdout.

TIFF studies. The wall is tilted in the image (Study 2: ~ -12 deg), and h(y,t) contains that slope,
which would dominate any roughness measurement. Here every frame is put on the base plane:
  1. theta0 = mean tilt over the run (tiff_extra.npz: theta_deg); one angle for all frames;
  2. the domain of the frame (arrival <= k, connected to the left edge, as in step1) is rotated by
     theta0 as a float image with linear interpolation and thresholded at 0.5 (rotating the binary
     map with nearest neighbours creates staircase edges and spurious overhangs);
  3. along each row s of the rotated frame (now parallel to the motion) three heights are taken:
       u_area  left edge of the valid region + length inside the domain (area-conserving "column
               height", single-valued by construction; equal to the true height without overhangs)
       u_front last crossing of the wall (front of overhangs)
       u_back  first crossing of the wall
     Overhangs are gaps of the domain at least 2 px long between u_back and u_front;
  4. only rows s where the wall is inside the rotated image in every frame used are kept, and the
     residual slope of each frame (theta(t) - theta0) is removed with a straight-line fit.
Observables, averaged over frames (and over thirds of the run):
  S(q)  = |FFT(du * Hann)|^2 / sum(Hann^2) ~ q^-(1+2 zeta)
  B(r)  = <[du(s+r) - du(s)]^2> ~ r^(2 zeta);  zeta_eff(r) = (1/2) d ln B / d ln r
  w2(l) = <variance of du in windows of length l after a linear fit> ~ l^(2 zeta)
  B(r) also with the pairs whose wall point is closer than --mask-d px to a static defect excluded.
zeta is fitted for rmin < r < L/4 and for the corresponding q range, on S(q) averaged in
logarithmic bins; its error comes from a block bootstrap over time (12 blocks of frames, which are
strongly correlated). rmin must be above the optical resolution (~8 px in Study 2) and above the
scale where the three heights start to differ because of overhangs. Self-test (synthetic walls,
--selftest): S(q) is unbiased, B(r) underestimates zeta >= 0.66 (0.58 for 0.66, 0.78 for 1) and
w(l) underestimates zeta ~ 1; S(q) is the reference estimator. Overhangs only distort scales up to their size: where u_area,
u_front and u_back give the same S(q), the exponent is not affected by them.
Local width with local rotation: the wall of each frame is cut into segments of projected length l
along the base plane (windows sliding by l/2); each segment is rotated by its own tilt, i.e.
w2 = smallest eigenvalue of the covariance of its points (orthogonal least-squares line), and
<w2(l)> is averaged over all segments and frames (~ l^(2 zeta)). For TIFF studies the points are
the full subpixel contour of the wall, resampled at uniform arc length, so overhangs are part of the
shape and need no choice of height; control: the same on u_area. Also the spread of the local
tilt of the segments vs l. Error: block bootstrap over frames.
Local structure factor S(q,l): the same segments (4 rmin <= l <= 0.8 L), each rotated by its own tilt; the
points are projected on the segment's own axis and normal, binned at 1 px along the axis (mean
normal displacement per bin, so a fold is averaged), the line through its end points is removed
(no window: for short segments a Hann window biases zeta upwards, end matching does not) and
S(q,l) = |FFT(profile)|^2 / N, averaged in log bins of q over segments and frames.
zeta(l) from S ~ q^-(1+2 zeta) over 2 pi / l < q < 2 pi / rmin (all the scales of the segment
above the resolution) and over the lowest band 2 pi / l < q < --sq-band x that (scales ~ l);
both with all segments and without folded ones (--max-fold). Compared with the slope of the global
S(q) in the same band. Error: block bootstrap over frames.
Height distribution: P(du/sigma) of the deviations du(s,t) = u(s,t) - <u>_s(t) from the base plane of
each frame (rotation by theta0 + per-frame line = rotation by theta(t)), normalised by the width
sigma(t) of each frame, with skewness and excess kurtosis (block-bootstrap errors). Also shown:
the deviation from the mean only (no per-frame line, the residual tilt broadens it), the three
heights (overhangs act on the tails) and local distributions in windows of 5 and 20 um with a
line removed in each window.
--from-h uses DIR/h_xt_sub.npy directly (no rotation, no overhangs), with the per-frame line fit.
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import binary_fill_holes, distance_transform_edt, gaussian_filter, label, rotate

ap = argparse.ArgumentParser()
ap.add_argument("dir", nargs="?", help="study folder")
ap.add_argument("--every", type=int, default=5, help="use every N-th frame")
ap.add_argument("--rmin", type=float, default=16, help="smallest scale used in the fits, px")
ap.add_argument("--rmax-frac", type=float, default=0.25, help="largest scale in the fits, as a fraction of L")
ap.add_argument("--mask-d", type=float, default=20, help="exclusion distance to defects for B(r), px")
ap.add_argument("--from-h", action="store_true", help="use h_xt_sub.npy (no rotation, no overhangs)")
ap.add_argument("--max-tilt", type=float, default=45,
                help="local width: also report the result without segments whose own tilt differs from the "
                     "base plane by more than this (deg); they are folds/tongues of the wall")
ap.add_argument("--max-fold", type=float, default=1.5,
                help="local width: also report the result without folded segments, whose contour length inside "
                     "the window exceeds this factor times the median contour length at the same l")
ap.add_argument("--sq-band", type=float, default=8,
                help="local S(q,l): width (factor in q) of the low-q band used for zeta at scale ~ l")
ap.add_argument("--selftest", action="store_true", help="synthetic tilted wall with known zeta")
args = ap.parse_args()


# ------------------------------------------------------------------ estimators (shared with step8/9)
from wall_tools import (contour_points, detrend, fit_power, heights_rotated, logbin_xy, pca_width, q_edges,
                        rotation_sign, seg_spectra, structure_factor, zeta_all, zeta_band, zeta_S)

# ------------------------------------------------------------------ self test
if args.selftest:
    rng = np.random.default_rng(0)
    H, W, theta = 1006, 1197, -12.0
    Y, X = np.mgrid[0:H, 0:W]
    sgn = rotation_sign(theta, (H, W))
    valid = rotate(np.ones((H, W)), sgn * theta, reshape=True, order=1, cval=0) > 0.999
    for zeta in (0.5, 0.66, 1.0):
        rot, ref = [], []
        for _ in range(12):
            k = np.fft.rfftfreq(H)
            amp = np.zeros_like(k); amp[1:] = k[1:] ** (-(1 + 2 * zeta) / 2)
            u = np.fft.irfft(amp * np.exp(2j * np.pi * rng.random(k.size)), n=H)
            u *= 8 / u.std()
            dom = X < 600 + (Y - H / 2) * np.tan(np.radians(theta)) + u[:, None]
            ua, uf, ub, gaps = heights_rotated(dom, dom, sgn * theta, valid)
            ok = np.nonzero(np.isfinite(ua))[0]
            run = max(np.split(ok, np.nonzero(np.diff(ok) != 1)[0] + 1), key=len)
            r1 = zeta_all(detrend(ua[run][None, :])[0], args.rmin, len(run) * args.rmax_frac)
            r0 = zeta_all(detrend(u[None, :])[0], args.rmin, H * args.rmax_frac)
            rot.append([r1["zS"], r1["zB"], r1["zW"]]); ref.append([r0["zS"], r0["zB"], r0["zW"]])
        rot, ref = np.array(rot), np.array(ref)
        print(f"selftest zeta = {zeta}: rotated wall (tilt {theta}°)  S(q) {rot[:, 0].mean():.2f}±{rot[:, 0].std():.2f}  "
              f"B(r) {rot[:, 1].mean():.2f}±{rot[:, 1].std():.2f}  w(l) {rot[:, 2].mean():.2f}±{rot[:, 2].std():.2f}"
              f"   | unrotated input: {ref[:, 0].mean():.2f} {ref[:, 1].mean():.2f} {ref[:, 2].mean():.2f}")
        zp = []
        ells_t = np.unique(np.round(np.logspace(np.log10(8), np.log10(H / 2), 20)).astype(int))
        for _ in range(6):
            k = np.fft.rfftfreq(H)
            amp = np.zeros_like(k); amp[1:] = k[1:] ** (-(1 + 2 * zeta) / 2)
            u = np.fft.irfft(amp * np.exp(2j * np.pi * rng.random(k.size)), n=H)
            u *= 8 / u.std()
            dom = X < 600 + (Y - H / 2) * np.tan(np.radians(theta)) + u[:, None]
            yy_, xx_ = contour_points(dom)
            w2l, _, _ = pca_width(yy_, xx_, theta, ells_t)
            zp.append(fit_power(ells_t.astype(float), np.array([v.mean() for v in w2l]), args.rmin, H * args.rmax_frac)[0] / 2)
        print(f"   local width with local rotation (PCA on the contour): zeta {np.mean(zp):.2f}±{np.std(zp):.2f}")
        ed = q_edges(H)
        qc_t = np.sqrt(ed[1:] * ed[:-1])
        for l_t in (100, 200, 400):
            zf, zl = [], []
            for _ in range(6):
                k = np.fft.rfftfreq(H)
                amp = np.zeros_like(k); amp[1:] = k[1:] ** (-(1 + 2 * zeta) / 2)
                u = np.fft.irfft(amp * np.exp(2j * np.pi * rng.random(k.size)), n=H)
                u *= 8 / u.std()
                dom = X < 600 + (Y - H / 2) * np.tan(np.radians(theta)) + u[:, None]
                yy_, xx_ = contour_points(dom)
                Ss, Cs, _ = seg_spectra(yy_, xx_, theta, l_t, 0.5, ed)
                Sm = Ss.sum(0) / np.maximum(Cs.sum(0), 1)
                q0 = 2 * np.pi / l_t
                zf.append(zeta_band(qc_t, Sm, q0, 2 * np.pi / args.rmin))
                zl.append(zeta_band(qc_t, Sm, q0, min(args.sq_band * q0, 2 * np.pi / args.rmin)))
            print(f"   local S(q,l), l = {l_t} px: zeta all q above resolution {np.mean(zf):.2f}±{np.std(zf):.2f}, "
                  f"low-q band {np.mean(zl):.2f}±{np.std(zl):.2f}")
    raise SystemExit

# ------------------------------------------------------------------ data
with open(os.path.join(args.dir, "meta.json")) as f:
    meta = json.load(f)
k_um, fps = meta.get("um_per_px") or 0.0, meta["fps"]
unit = "µm" if k_um else "px"
kk = k_um if k_um else 1.0

if args.from_h:
    h = np.load(os.path.join(args.dir, "h_xt_sub.npy"))[::args.every]
    U = {"h": h}
    frames = np.arange(0, h.shape[0] * args.every, args.every)
    wmask = None
    example = None
    theta0 = 0.0
else:
    e = np.load(os.path.join(args.dir, "tiff_extra.npz"))
    ta, defects, theta = e["arrival"], e["defects"], e["theta_deg"]
    theta0 = float(np.mean(theta))
    sgn = rotation_sign(theta0, ta.shape)
    ang = sgn * theta0
    valid = rotate(np.ones(ta.shape), ang, reshape=True, order=1, cval=0) > 0.999
    dist_rot = rotate(distance_transform_edt(~defects), ang, reshape=True, order=1, cval=0)
    T = len(theta)
    frames = np.arange(0, T, args.every)
    ua_all, uf_all, ub_all, gp_all, contours, contours_s1 = [], [], [], [], [], []
    edge = lambda lab_: np.setdiff1d(np.unique(lab_[:, 0]), [0])
    for i, kf in enumerate(frames):
        dom = ta <= kf
        lab_, _ = label(dom)
        dom = np.isin(lab_, edge(lab_))
        domf = binary_fill_holes(dom)
        ua, uf, ub, gp = heights_rotated(dom, domf, ang, valid)
        ua_all.append(ua); uf_all.append(uf); ub_all.append(ub); gp_all.append(gp)
        contours.append(contour_points(domf))
        # for step9 (strictly monotonic local windows): contour of the domain smoothed by 1 px, without
        # the pixel staircase of the binary contour (below the optical resolution, ~3-4 px)
        contours_s1.append(contour_points(gaussian_filter(domf.astype(float), 1.0)))
        if i % 60 == 0:
            print(f"  frame {kf}/{T}")
    ua_all, uf_all, ub_all, gp_all = map(np.array, (ua_all, uf_all, ub_all, gp_all))
    ok = np.isfinite(ua_all).all(0)
    idx = np.nonzero(ok)[0]
    runs = np.split(idx, np.nonzero(np.diff(idx) != 1)[0] + 1)
    run = max(runs, key=len)
    U = {"area": ua_all[:, run], "front": uf_all[:, run], "back": ub_all[:, run]}
    gaps = gp_all[:, run]
    # defect distance of the wall point (u_area, s) for the masked B(r)
    cols = np.clip(np.round(U["area"]).astype(int), 0, dist_rot.shape[1] - 1)
    wmask = (dist_rot[run[None, :], cols] > args.mask_d).astype(float)
    example = (frames[len(frames) // 2], len(frames) // 2, run)
    # ordered subpixel contours of the frames used (for step9, wavelets in local frames)
    off = np.cumsum([0] + [c_[0].size for c_ in contours_s1])
    np.savez_compressed(os.path.join(args.dir, "contours.npz"), y=np.concatenate([c_[0] for c_ in contours_s1]),
                        x=np.concatenate([c_[1] for c_ in contours_s1]), offsets=off, frames=frames, theta0=theta0,
                        step=0.5, smoothing_px=1.0)
    print(f"theta0 = {theta0:.2f} deg (rotation {ang:+.2f}); rows kept on the base plane: {run.size} "
          f"({run.size * kk:.0f} {unit}); frames used: {len(frames)}")

L = next(iter(U.values())).shape[1]
rmax = L * args.rmax_frac
res = {}
for name, u in U.items():
    du, slopes = detrend(u)
    w_ = wmask if (name == "area" and wmask is not None) else None
    res[name] = zeta_all(du, args.rmin, rmax, w_)
    res[name]["tilt_res"] = np.degrees(np.arctan(slopes))
    res[name]["width"] = du.std(1)
    thirds = np.array_split(np.arange(du.shape[0]), 3)
    res[name]["S_thirds"] = [structure_factor(du[t_])[1] for t_ in thirds]
    res[name]["zS_thirds"] = [zeta_S(res[name]["q"], s_, args.rmin, rmax) for s_ in res[name]["S_thirds"]]
    # block bootstrap over time (frames are strongly correlated): 12 blocks resampled with replacement
    blocks = np.array_split(np.arange(du.shape[0]), 12)
    Sblk = [structure_factor(du[b])[1] for b in blocks]
    res[name]["Sblk"] = np.array(Sblk)                    # per time block, for step8
    rng = np.random.default_rng(1)
    zb = [zeta_S(res[name]["q"], np.mean([Sblk[j] for j in rng.integers(0, 12, 12)], axis=0), args.rmin, rmax)
          for _ in range(300)]
    res[name]["dzS"] = float(np.std(zb))

main = "area" if "area" in res else "h"
print(f"residual tilt per frame after rotation: {np.abs(res[main]['tilt_res']).max():.2f} deg max; "
      f"rms width of the detrended wall: {res[main]['width'].mean():.2f} px"
      + (f" = {res[main]['width'].mean() * k_um:.2f} µm" if k_um else ""))
if not args.from_h:
    oh = gaps > 0
    d_fb = U["front"] - U["back"]
    ext = []
    for row in oh:
        dd = np.diff(np.concatenate([[0], row.astype(np.int8), [0]]))
        ext += list(np.nonzero(dd == -1)[0] - np.nonzero(dd == 1)[0])
    ext = np.array(ext) if ext else np.array([0])
    print(f"overhangs on the base plane: {oh.mean():.1%} of (row, frame); extent along the wall median "
          f"{np.median(ext):.0f} px, 90% {np.percentile(ext, 90):.0f} px, max {ext.max()}; depth (front - back) "
          f"median {np.median(d_fb[oh]) if oh.any() else 0:.0f} px, 90% {np.percentile(d_fb[oh], 90) if oh.any() else 0:.0f} px")
print(f"fit range: {args.rmin:g} < r < {rmax:.0f} px" + (f" ({args.rmin * k_um:.1f}-{rmax * k_um:.0f} µm)" if k_um else ""))
for name, r_ in res.items():
    extra = f"  B(r) far from defects {r_['zBm']:.3f}" if r_["Bm"] is not None else ""
    print(f"  {name:6s} zeta: S(q) {r_['zS']:.3f} ± {r_['dzS']:.3f}  B(r) {r_['zB']:.3f}  w(l) {r_['zW']:.3f}{extra}  "
          f"| S(q) by thirds of the run: " + ", ".join(f"{z:.3f}" for z in r_["zS_thirds"]))
np.savez(os.path.join(args.dir, "roughness.npz"), theta0=theta0, frames=frames, every=args.every,
         rmin=args.rmin, rmax=rmax, um_per_px=k_um,
         **{f"{n}_{k}": v for n, r_ in res.items() for k, v in r_.items()
            if k in ("q", "S", "r", "B", "Bm", "l", "W2", "zS", "dzS", "zB", "zBm", "zW", "zS_thirds", "Sblk") and v is not None})

# ------------------------------------------------------------------ local width, local rotation
ells = np.unique(np.round(np.logspace(np.log10(4), np.log10(L), 28)).astype(int)).astype(float)
if args.from_h:
    pts = [(np.arange(L, dtype=float), row) for row in U["h"]]
    th_pts = 0.0
else:
    pts = contours
    th_pts = theta0
per_frame = [pca_width(yy_, xx_, th_pts, ells) for yy_, xx_ in pts]     # [(w2, tilt, arc/l lists)]
W2f = np.array([[v.mean() if v.size else np.nan for v in pf[0]] for pf in per_frame])   # frames x ells
NW = np.array([[v.size for v in pf[0]] for pf in per_frame])
# same, without the segments tilted by more than --max-tilt from the base plane (folds, tongues)
keepT = [[np.abs(t_) <= args.max_tilt for t_ in pf[1]] for pf in per_frame]
W2r = np.array([[v[k_].mean() if k_.any() else np.nan for v, k_ in zip(pf[0], kp)] for pf, kp in zip(per_frame, keepT)])
NR = np.array([[k_.sum() for k_ in kp] for kp in keepT])
w2_restr = np.nansum(np.nan_to_num(W2r) * NR, 0) / np.maximum(NR.sum(0), 1)
frac_excl = 1 - NR.sum(0) / np.maximum(NW.sum(0), 1)
# same, without folded segments: contour inside the window much longer than typical at that l
# (relative to the median, because sub-window wiggles already make arc/l > 1 at every scale)
arc_med = np.array([np.median(np.concatenate([pf[2][j] for pf in per_frame])) if NW[:, j].sum() else np.nan
                    for j in range(len(ells))])
keepF = [[a_ <= args.max_fold * arc_med[j] for j, a_ in enumerate(pf[2])] for pf in per_frame]
W2F = np.array([[v[k_].mean() if k_.any() else np.nan for v, k_ in zip(pf[0], kp)] for pf, kp in zip(per_frame, keepF)])
NF = np.array([[k_.sum() for k_ in kp] for kp in keepF])
w2_fold = np.nansum(np.nan_to_num(W2F) * NF, 0) / np.maximum(NF.sum(0), 1)
frac_fold = 1 - NF.sum(0) / np.maximum(NW.sum(0), 1)
w2_pca = np.nansum(W2f * NW, 0) / np.maximum(NW.sum(0), 1)
tilt_sd = np.array([np.std(np.concatenate([pf[1][j] for pf in per_frame])) if NW[:, j].sum() else np.nan
                    for j in range(len(ells))])
blocks = np.array_split(np.arange(W2f.shape[0]), 12)
rng_w = np.random.default_rng(3)
boot = []
for _ in range(300):
    pick = np.concatenate([blocks[j] for j in rng_w.integers(0, 12, 12)])
    boot.append(np.nansum(W2f[pick] * NW[pick], 0) / np.maximum(NW[pick].sum(0), 1))
boot = np.array(boot)
dw2 = boot.std(0)
zfit = lambda w2: fit_power(ells, w2, args.rmin, rmax)[0] / 2
z_pca = zfit(w2_pca)
dz_pca = np.std([zfit(b_) for b_ in boot])
ref_name = "area" if "area" in U else "h"
# control: same PCA on the single-valued height of the base plane (rows s, height u)
ctrl = [pca_width(np.arange(L, dtype=float), row, 0.0, ells) for row in U[ref_name]]
w2_ctrl = np.array([np.mean(np.concatenate([c[0][j] for c in ctrl])) if any(c[0][j].size for c in ctrl) else np.nan
                    for j in range(len(ells))])
z_ctrl = zfit(w2_ctrl)
# per time block sums and counts of w2 (12 blocks of frames), for the bootstrap of step8
csum = np.array([[c[0][j].sum() for j in range(len(ells))] for c in ctrl])
ccnt = np.array([[c[0][j].size for j in range(len(ells))] for c in ctrl])
blk_sum = {k_: np.array([np.nansum(np.nan_to_num(W_[b]) * N_[b], 0) for b in blocks])
           for k_, W_, N_ in (("pca", W2f, NW), ("unfolded", W2F, NF))}
blk_sum["ctrl"] = np.array([csum[b].sum(0) for b in blocks])
blk_cnt = {"pca": np.array([NW[b].sum(0) for b in blocks]), "unfolded": np.array([NF[b].sum(0) for b in blocks]),
           "ctrl": np.array([ccnt[b].sum(0) for b in blocks])}
z_shear = res[ref_name]["zW"]
z_restr = zfit(w2_restr)
z_fold = zfit(w2_fold)
print(f"local width w2(l), each segment rotated by its own tilt (PCA), {'contour' if not args.from_h else 'h'}: "
      f"zeta = {z_pca:.3f} ± {dz_pca:.3f}; same on u_{ref_name}: {z_ctrl:.3f}; "
      f"line removed without rotation (w(l) above): {z_shear:.3f}")
print(f"   without segments tilted > {args.max_tilt:g} deg from the base plane: zeta = {z_restr:.3f} "
      f"(segments excluded: {frac_excl[(ells >= args.rmin) & (ells <= rmax)].mean():.0%} in the fit range)")
print(f"   without folded segments (contour > {args.max_fold:g} x median at that l): zeta = {z_fold:.3f} "
      f"(segments excluded: {frac_fold[(ells >= args.rmin) & (ells <= rmax)].mean():.0%} in the fit range)")
for j in np.unique(np.searchsorted(ells, [8, 25, 85, 250, 850])).clip(0, len(ells) - 1):
    print(f"   l = {ells[j]:5.0f} px ({ells[j] * kk:6.1f} {unit}): w = {np.sqrt(w2_pca[j]) * kk:6.3f} {unit}, "
          f"local tilt sd {tilt_sd[j]:5.1f} deg, {int(NW[:, j].sum())} segments, "
          f"{frac_excl[j]:.0%} tilted > {args.max_tilt:g} deg, median contour/l {arc_med[j]:.2f}, "
          f"{frac_fold[j]:.0%} folded")
np.savez(os.path.join(args.dir, "local_width.npz"), ells=ells, w2_pca=w2_pca, dw2=dw2, w2_ctrl=w2_ctrl,
         w2_restricted=w2_restr, frac_excluded=frac_excl, max_tilt=args.max_tilt, zeta_restricted=z_restr,
         w2_unfolded=w2_fold, frac_folded=frac_fold, arc_median=arc_med, max_fold=args.max_fold, zeta_unfolded=z_fold,
         tilt_sd=tilt_sd, n_segments=NW.sum(0), zeta_pca=z_pca, dzeta_pca=dz_pca, zeta_ctrl=z_ctrl,
         theta0=th_pts, every=args.every, rmin=args.rmin, rmax=rmax, um_per_px=k_um,
         **{f"blk_sum_{k_}": v for k_, v in blk_sum.items()}, **{f"blk_cnt_{k_}": v for k_, v in blk_cnt.items()})

fig, ax = plt.subplots(1, 3, figsize=(18, 5))
yy_, xx_ = pts[len(pts) // 2]
jshow = np.argmin(np.abs(ells * kk - (20 if k_um else 160)))
lshow = ells[jshow]
t0r = np.radians(th_pts)
sx = yy_ * np.cos(t0r) + xx_ * np.sin(t0r)
ax[0].plot(xx_, yy_, ".", ms=0.6, color="0.6")
for j, s0 in enumerate(np.arange(sx.min(), sx.max() - lshow, lshow)):
    m = (sx >= s0) & (sx < s0 + lshow)
    if m.sum() < 4:
        continue
    py, px = yy_[m], xx_[m]
    cy, cx = py.mean(), px.mean()
    cov = np.cov(np.vstack([py, px]))
    ev, evec = np.linalg.eigh(cov)
    d = evec[:, 1]
    tt = np.array([-lshow / 2, lshow / 2]) / max(abs(d[0] * np.cos(t0r) + d[1] * np.sin(t0r)), 0.2)
    c = plt.cm.tab10(j % 10)
    folded = m.sum() * 0.5 / lshow > args.max_fold * arc_med[jshow]
    ax[0].plot(px, py, ".", ms=1.2, color=c)
    ax[0].plot(cx + tt * d[1], cy + tt * d[0], ":" if folded else "-", color="r" if folded else "k", lw=1.2)
ax[0].invert_yaxis(); ax[0].set_aspect("equal")
ax[0].set_title(f"Segments of {lshow * kk:.0f} {unit} and their own best lines (frame {frames[len(pts) // 2]})\n"
                f"red dotted: folded (contour > {args.max_fold:g} x median), excluded", fontsize=10)
ax[0].set_xlabel("x [px]" if not args.from_h else "h [px]"); ax[0].set_ylabel("y [px]" if not args.from_h else "column")
ax[1].loglog(ells * kk, w2_pca * kk ** 2, "o-", ms=3, color="C0", label=f"local rotation, {'contour' if not args.from_h else 'h'}: ζ = {z_pca:.2f} ± {dz_pca:.2f}")
ax[1].fill_between(ells * kk, (w2_pca - 2 * dw2) * kk ** 2, (w2_pca + 2 * dw2) * kk ** 2, color="C0", alpha=0.2)
ax[1].loglog(ells * kk, w2_restr * kk ** 2, "D-", ms=3, mfc="none", color="C4",
             label=f"local rotation, segments tilted ≤ {args.max_tilt:g}°: ζ = {z_restr:.2f}")
ax[1].loglog(ells * kk, w2_fold * kk ** 2, "v-", ms=3, mfc="none", color="C3",
             label=f"local rotation, folded segments excluded: ζ = {z_fold:.2f}")
ax[1].loglog(ells * kk, w2_ctrl * kk ** 2, "s--", ms=3, mfc="none", color="C1", label=f"local rotation, u_{ref_name}: ζ = {z_ctrl:.2f}")
ax[1].loglog(res[ref_name]["l"] * kk, res[ref_name]["W2"] * kk ** 2, "^:", ms=3, color="C2", label=f"line removed, no rotation: ζ = {z_shear:.2f}")
ax[1].axvspan(args.rmin * kk, rmax * kk, color="0.92", zorder=0)
ax[1].set_xlabel(f"segment length ℓ [{unit}]"); ax[1].set_ylabel(f"⟨w²(ℓ)⟩ [{unit}²]"); ax[1].legend(fontsize=8)
ax[1].set_title(r"Local width, $\langle w^2\rangle \sim \ell^{2\zeta}$ (grey: fit range; band: ±2σ)", fontsize=10)
ax[2].semilogx(ells * kk, tilt_sd, "o-", ms=3, label="sd of the local tilt")
ax2b = ax[2].twinx()
ax2b.semilogx(ells * kk, 100 * frac_excl, "s--", ms=3, color="C4", mfc="none", label=f"% tilted > {args.max_tilt:g}°")
ax2b.semilogx(ells * kk, 100 * frac_fold, "v--", ms=3, color="C3", mfc="none", label=f"% folded (contour > {args.max_fold:g} x median)")
ax2b.set_ylabel("% of segments excluded"); ax2b.legend(fontsize=8, loc="upper right")
ax[2].set_xlabel(f"segment length ℓ [{unit}]"); ax[2].set_ylabel("sd of the local tilt [deg]")
ax[2].set_title("Spread of the local tilt of the segments", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(args.dir, "fig_local_width.png"), dpi=120); plt.close(fig)

# ------------------------------------------------------------------ local structure factor S(q, l)
step_pts = 0.5                     # same convention as pca_width, so contour/l matches arc_med
edges = q_edges(L)
qc = np.sqrt(edges[1:] * edges[:-1])
qres = 2 * np.pi / args.rmin
jl = np.nonzero((ells >= 4 * args.rmin) & (ells <= 0.8 * L))[0]
nblk = 12
blk_of = np.repeat(np.arange(nblk), [len(b_) for b_ in np.array_split(np.arange(len(pts)), nblk)])
SQ = {}                      # l -> dict(all=(Ssum_blocks, C_blocks), unf=(...), nseg, frac_fold)
for j in jl:
    l = ells[j]
    Sb = {k_: np.zeros((nblk, len(qc))) for k_ in ("all", "unf")}
    Cb = {k_: np.zeros((nblk, len(qc))) for k_ in ("all", "unf")}
    nseg = nfold = 0
    for f_, (yy_, xx_) in enumerate(pts):
        Ss, Cs, arc = seg_spectra(yy_, xx_, th_pts, l, step_pts, edges)
        if not arc.size:
            continue
        unf = arc <= args.max_fold * arc_med[j]
        b_ = blk_of[f_]
        Sb["all"][b_] += Ss.sum(0); Cb["all"][b_] += Cs.sum(0)
        Sb["unf"][b_] += Ss[unf].sum(0); Cb["unf"][b_] += Cs[unf].sum(0)
        nseg += arc.size; nfold += (~unf).sum()
    SQ[l] = dict(Sb=Sb, Cb=Cb, nseg=nseg, frac_fold=nfold / max(nseg, 1))


def sq_zetas(Sb, Cb, l, pick=None):
    pick = np.arange(nblk) if pick is None else pick
    Sm = Sb[pick].sum(0) / np.maximum(Cb[pick].sum(0), 1)
    Sm = np.where(Cb[pick].sum(0) > 0, Sm, np.nan)
    q0 = 2 * np.pi / l
    return Sm, zeta_band(qc, Sm, q0, qres), zeta_band(qc, Sm, q0, min(args.sq_band * q0, qres))


qg, Sg = logbin_xy(res[ref_name]["q"], res[ref_name]["S"])
rng_q = np.random.default_rng(5)
ls_sq = np.array(sorted(SQ))
Z = {k_: np.full((len(ls_sq), 2), np.nan) for k_ in ("all", "unf")}
dZ = {k_: np.full((len(ls_sq), 2), np.nan) for k_ in ("all", "unf")}
Sq_mean = {k_: [] for k_ in ("all", "unf")}
zg_band = np.full(len(ls_sq), np.nan)
for i, l in enumerate(ls_sq):
    for k_ in ("all", "unf"):
        Sm, zf, zl = sq_zetas(SQ[l]["Sb"][k_], SQ[l]["Cb"][k_], l)
        Sq_mean[k_].append(Sm)
        Z[k_][i] = zf, zl
        bt = np.array([sq_zetas(SQ[l]["Sb"][k_], SQ[l]["Cb"][k_], l, rng_q.integers(0, nblk, nblk))[1:]
                       for _ in range(200)])
        dZ[k_][i] = np.nanstd(bt, 0)
    q0 = 2 * np.pi / l
    zg_band[i] = zeta_band(qg, Sg, q0, min(args.sq_band * q0, qres))
# local slope of the local width (folds excluded), over a factor 3 in l, for comparison
zw_loc = np.array([fit_power(ells, w2_fold, l / np.sqrt(3), l * np.sqrt(3))[0] / 2 for l in ells])
print(f"local structure factor S(q,l), segments rotated by their own tilt (fit: 2pi/l < q < 2pi/{args.rmin:g} px; "
      f"low band: factor {args.sq_band:g} above 2pi/l):")
for i, l in enumerate(ls_sq):
    print(f"   l = {l:5.0f} px ({l * kk:6.1f} {unit}): {SQ[l]['nseg']:6d} segments, {SQ[l]['frac_fold']:.0%} folded | "
          f"zeta all q: {Z['all'][i, 0]:.2f}±{dZ['all'][i, 0]:.2f} (unfolded {Z['unf'][i, 0]:.2f}±{dZ['unf'][i, 0]:.2f}) | "
          f"low band: {Z['all'][i, 1]:.2f}±{dZ['all'][i, 1]:.2f} (unfolded {Z['unf'][i, 1]:.2f}±{dZ['unf'][i, 1]:.2f}); "
          f"global S(q) in the same band {zg_band[i]:.2f}")
np.savez(os.path.join(args.dir, "local_sq.npz"), q=qc, ells=ls_sq, S_all=np.array(Sq_mean["all"]),
         S_unfolded=np.array(Sq_mean["unf"]), zeta_all=Z["all"], dzeta_all=dZ["all"], zeta_unfolded=Z["unf"],
         dzeta_unfolded=dZ["unf"], zeta_global_band=zg_band, n_segments=[SQ[l]["nseg"] for l in ls_sq],
         frac_folded=[SQ[l]["frac_fold"] for l in ls_sq], q_res=qres, sq_band=args.sq_band, um_per_px=k_um,
         note="zeta columns: [fit over 2pi/l < q < q_res, fit over the low band 2pi/l < q < sq_band*2pi/l]")

fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
show = ls_sq[np.unique(np.round(np.linspace(0, len(ls_sq) - 1, min(7, len(ls_sq)))).astype(int))]
for i, l in enumerate(ls_sq):
    if l not in show:
        continue
    c = plt.cm.viridis(np.log(l / ls_sq[0]) / max(np.log(ls_sq[-1] / ls_sq[0]), 1e-9) * 0.9)
    Sm = Sq_mean["unf"][i]
    m = np.isfinite(Sm) & (qc >= 0.9 * 2 * np.pi / l)
    ax[0].loglog(qc[m] / kk, Sm[m] * kk ** 3, "o-", ms=3, color=c,
                 label=f"ℓ = {l * kk:.0f} {unit}: ζ = {Z['unf'][i, 0]:.2f}")
    q0 = 2 * np.pi / l
    ax[0].axvline(q0 / kk, color=c, lw=0.6, ls=":")
ax[0].loglog(qg / kk, Sg * kk ** 3, "k-", lw=1.2, alpha=0.7, label=f"whole wall (u_{ref_name}), global S(q)")
ax[0].axvline(qres / kk, color="0.4", ls="--", lw=0.8)
ax[0].text(qres / kk, ax[0].get_ylim()[1], " resolution", va="top", fontsize=8, color="0.4")
ax[0].set_xlabel(f"q [{unit}$^{{-1}}$]"); ax[0].set_ylabel(f"S(q, ℓ) [{unit}$^3$]")
ax[0].set_title("Local structure factor, segments rotated by their own tilt (folds excluded)\n"
                "dotted: q = 2π/ℓ, lower end of the fits", fontsize=10)
ax[0].legend(fontsize=7)
for k_, mk, lab in (("all", "o", "all segments"), ("unf", "v", "folds excluded")):
    ax[1].errorbar(ls_sq * kk, Z[k_][:, 0], dZ[k_][:, 0], fmt=mk + "-", ms=4, capsize=2, mfc="none" if k_ == "all" else None,
                   label=f"S(q,ℓ), all q above resolution, {lab}")
    ax[1].errorbar(ls_sq * kk, Z[k_][:, 1], dZ[k_][:, 1], fmt=mk + "--", ms=4, capsize=2, mfc="none" if k_ == "all" else None,
                   label=f"S(q,ℓ), low-q band (scales ~ ℓ), {lab}")
ax[1].plot(ls_sq * kk, zg_band, "ks:", ms=4, mfc="none", label="global S(q), same low-q band")
okw = (ells >= 2 * args.rmin) & (ells <= L / 2)
ax[1].plot(ells[okw] * kk, zw_loc[okw], "-", color="0.6", lw=1.5, label="local width (folds excluded), local slope")
for zz, lab in ((0.5, "1/2"), (2 / 3, "2/3"), (1.25, "1.25")):
    ax[1].axhline(zz, color="0.8", lw=0.8, zorder=0)
    ax[1].text(ax[1].get_xlim()[0] if False else ls_sq[0] * kk * 0.8, zz, lab, fontsize=8, color="0.5", va="bottom")
ax[1].set_xscale("log"); ax[1].set_xlabel(f"segment length ℓ [{unit}]"); ax[1].set_ylabel("ζ(ℓ)")
ax[1].set_title(r"Exponent from $S(q,\ell) \sim q^{-(1+2\zeta(\ell))}$ (error: block bootstrap)", fontsize=10)
ax[1].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(args.dir, "fig_local_sq.png"), dpi=120); plt.close(fig)

# ------------------------------------------------------------------ height distribution
def moments(x):
    x = x - x.mean()
    sd = x.std()
    return np.mean(x ** 3) / sd ** 3, np.mean(x ** 4) / sd ** 4 - 3


def local_residuals(du, l):
    """Residuals of du after a line fit in non-overlapping windows of length l, per frame."""
    l = int(l)
    out = []
    x = np.arange(l)
    A = np.vstack([x, np.ones(l)]).T
    P = A @ np.linalg.pinv(A)
    for s0 in range(0, du.shape[1] - l + 1, l):
        seg = du[:, s0:s0 + l]
        out.append(seg - seg @ P.T)
    return np.concatenate(out, axis=1)


ZB = np.linspace(-5, 5, 61)
hdist = {}
rng_h = np.random.default_rng(2)
for name, u in U.items():
    du, _ = detrend(u)                                   # deviation from the base plane of each frame
    variants = {"base plane": du, "mean only": u - u.mean(1, keepdims=True)}
    for l_um in (5, 20):
        l_px = l_um / kk if k_um else l_um * 8
        if l_px * 3 < du.shape[1]:
            variants[f"local {l_um:g} {unit}"] = local_residuals(du, l_px)
    for vname, d in variants.items():
        z = d / d.std(1, keepdims=True)                  # normalised by the width of each frame
        hist, _ = np.histogram(z.ravel(), ZB, density=True)
        sk, ku = moments(z.ravel())
        blocks = np.array_split(np.arange(z.shape[0]), 12)
        bs = []
        for _ in range(200):
            pick = np.concatenate([blocks[j] for j in rng_h.integers(0, 12, 12)])
            bs.append(moments(z[pick].ravel()))
        bs = np.array(bs)
        hdist[(name, vname)] = dict(hist=hist, skew=sk, kurt=ku, dskew=bs[:, 0].std(), dkurt=bs[:, 1].std(),
                                    width=d.std(1).mean())
print("height distribution (deviation from the mean of each frame, normalised by its width):")
for (name, vname), hd_ in hdist.items():
    print(f"  {name:6s} {vname:14s} width {hd_['width'] * kk:6.2f} {unit}  skewness {hd_['skew']:+.2f} ± {hd_['dskew']:.2f}  "
          f"excess kurtosis {hd_['kurt']:+.2f} ± {hd_['dkurt']:.2f}")
np.savez(os.path.join(args.dir, "height_distribution.npz"), bins=ZB,
         **{f"{n}|{v}|{k}": val for (n, v), d_ in hdist.items() for k, val in d_.items()})

fig, ax = plt.subplots(1, 3, figsize=(17, 4.6))
zc = 0.5 * (ZB[1:] + ZB[:-1])
gauss = np.exp(-zc ** 2 / 2) / np.sqrt(2 * np.pi)
mainh = "area" if "area" in U else "h"
for a_ in ax:
    a_.semilogy(zc, gauss, "k--", lw=1, label="Gaussian")
for j, (vname, c) in enumerate((("base plane", "C0"), ("mean only", "C7"))):
    d_ = hdist[(mainh, vname)]
    ax[0].semilogy(zc, np.where(d_["hist"] > 0, d_["hist"], np.nan), "o-", ms=3, color=c,
                   label=f"{vname}: S={d_['skew']:+.2f}, K={d_['kurt']:+.2f}")
ax[0].set_title(f"P(δu/σ), {mainh}: base plane vs mean only", fontsize=10)
for name, c in (("area", "C0"), ("front", "C3"), ("back", "C2"), ("h", "C0")):
    if (name, "base plane") in hdist:
        d_ = hdist[(name, "base plane")]
        ax[1].semilogy(zc, np.where(d_["hist"] > 0, d_["hist"], np.nan), "o-", ms=3, color=c,
                       label=f"{name}: S={d_['skew']:+.2f}, K={d_['kurt']:+.2f}")
ax[1].set_title("Overhang sensitivity (area / front / back)", fontsize=10)
for vname, c in ((f"local 5 {unit}", "C1"), (f"local 20 {unit}", "C4"), ("base plane", "C0")):
    if (mainh, vname) in hdist:
        d_ = hdist[(mainh, vname)]
        ax[2].semilogy(zc, np.where(d_["hist"] > 0, d_["hist"], np.nan), "o-", ms=3, color=c,
                       label=f"{vname}: S={d_['skew']:+.2f}, K={d_['kurt']:+.2f}")
ax[2].set_title("Local (windows, line removed) vs whole wall", fontsize=10)
for a_ in ax:
    a_.set_xlabel(r"$\delta u/\sigma$"); a_.set_ylabel("probability density"); a_.set_ylim(1e-5, 1); a_.legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(args.dir, "fig_height_distribution.png"), dpi=120); plt.close(fig)

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(1, 4, figsize=(22, 5))
cols = {"area": "C0", "front": "C3", "back": "C2", "h": "C0"}
if example is not None:
    kf, i_, run_ = example
    dom = ta <= kf
    lab_, _ = label(dom)
    dom = np.isin(lab_, np.setdiff1d(np.unique(lab_[:, 0]), [0]))
    R = rotate(dom.astype(float), ang, reshape=True, order=1, cval=np.nan)
    ax[0].imshow(R, cmap="gray", vmin=-0.2, vmax=1.2, aspect="auto")
    for name in ("area", "front", "back"):
        ax[0].plot(U[name][i_], run_, "-", lw=0.8, color=cols[name], label=f"u_{name}")
    ax[0].set_title(f"frame {kf} rotated by {ang:+.1f}° (base plane vertical)", fontsize=10)
    ax[0].legend(fontsize=8); ax[0].set_xlabel("u [px]"); ax[0].set_ylabel("s [px]")
else:
    ax[0].plot(res["h"]["width"]); ax[0].set_xlabel("frame index"); ax[0].set_ylabel("rms width [px]")
    ax[0].set_title("Width of the detrended wall", fontsize=10)
for name, r_ in res.items():
    ax[1].loglog(r_["q"] / kk, r_["S"] * kk ** 3, "-", lw=0.6, alpha=0.5, color=cols[name])
    qb, Sb = logbin_xy(r_["q"], r_["S"])
    ax[1].loglog(qb / kk, Sb * kk ** 3, "o", ms=3, color=cols[name], label=f"{name}: ζ = {r_['zS']:.2f} ± {r_['dzS']:.2f}")
qa, qb = 2 * np.pi / rmax, 2 * np.pi / args.rmin
ax[1].axvspan(qa / kk, qb / kk, color="0.92", zorder=0)
ax[1].set_xlabel(f"q [1/{unit}]"); ax[1].set_ylabel(f"S(q) [{unit}³]"); ax[1].legend(fontsize=8)
ax[1].set_title(r"Structure factor, $S \sim q^{-(1+2\zeta)}$ (grey: fit range)", fontsize=10)
r_ = res[main]
ax[2].loglog(r_["r"] * kk, r_["B"] * kk ** 2, "o-", ms=3, label=f"B(r): ζ = {r_['zB']:.2f}")
if r_["Bm"] is not None:
    ax[2].loglog(r_["r"] * kk, r_["Bm"] * kk ** 2, "s--", ms=3, mfc="none", label=f"far from defects: ζ = {r_['zBm']:.2f}")
ax[2].loglog(r_["l"] * kk, r_["W2"] * kk ** 2, "^:", ms=3, label=f"w²(ℓ): ζ = {r_['zW']:.2f}")
ax[2].axvspan(args.rmin * kk, rmax * kk, color="0.92", zorder=0)
ax[2].set_xlabel(f"r, ℓ [{unit}]"); ax[2].set_ylabel(f"[{unit}²]"); ax[2].legend(fontsize=8)
ax[2].set_title(f"Height differences and local width ({main})", fontsize=10)
axz = ax[2].inset_axes([0.58, 0.08, 0.38, 0.32])
lr = np.log(r_["r"]); zeff = 0.5 * np.gradient(np.log(r_["B"]), lr)
axz.semilogx(r_["r"] * kk, zeff, "-", color="k", lw=1)
for z_, ls_ in ((2 / 3, ":"), (1.25, "--")):
    axz.axhline(z_, color="0.5", ls=ls_, lw=0.8)
p2 = os.path.join(args.dir, "persistence_sub_rm_eps0.97.npz")
if os.path.exists(p2):
    pr = np.load(p2)
    for xv, c in ((pr["xi"][int(np.argmax(pr["chi4"]))], "C1"), (np.nanmax(pr["xi"]), "C3")):
        axz.axvline(xv * kk, color=c, lw=0.8)
axz.set_ylim(0, 1.6); axz.set_title(r"$\zeta_{eff}(r)$ (2/3, 1.25; ξ(τ*), ξ_max)", fontsize=7); axz.tick_params(labelsize=7)
for j, s_ in enumerate(res[main]["S_thirds"]):
    ax[3].loglog(r_["q"] / kk, s_ * kk ** 3, "-", lw=1, color=plt.cm.viridis(j / 2.5),
                 label=f"third {j + 1}: ζ = {res[main]['zS_thirds'][j]:.2f}")
ax[3].axvspan(qa / kk, qb / kk, color="0.92", zorder=0)
ax[3].set_xlabel(f"q [1/{unit}]"); ax[3].set_ylabel(f"S(q) [{unit}³]"); ax[3].legend(fontsize=8)
ax[3].set_title("S(q) in the three thirds of the run", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(args.dir, "fig_roughness.png"), dpi=120); plt.close(fig)
print(f"outputs: {args.dir}/roughness.npz, fig_roughness.png")
