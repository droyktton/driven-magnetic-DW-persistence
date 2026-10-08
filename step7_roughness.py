"""Step 7: roughness and structure factor of the wall, on its base plane.

Usage: python step7_roughness.py DIR [--every 5] [--rmin 16] [--mask-d 20]
       python step7_roughness.py DIR --from-h [--every 1]       (video studies, simulations)
       python step7_roughness.py --selftest                     (synthetic wall with known zeta)
Outputs in DIR: roughness.npz, fig_roughness.png, local_width.npz, fig_local_width.png,
height_distribution.npz, fig_height_distribution.png and a summary on stdout.

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
from scipy.ndimage import binary_fill_holes, distance_transform_edt, label, rotate

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
ap.add_argument("--selftest", action="store_true", help="synthetic tilted wall with known zeta")
args = ap.parse_args()


# ------------------------------------------------------------------ estimators
def detrend(u):
    """Remove a straight line from each row of u (frames x L); returns residuals and slopes."""
    s = np.arange(u.shape[1])
    A = np.vstack([s, np.ones_like(s)]).T
    coef, *_ = np.linalg.lstsq(A, u.T, rcond=None)
    return u - (A @ coef).T, coef[0]


def structure_factor(du):
    L = du.shape[1]
    win = np.hanning(L)
    U = np.fft.rfft(du * win, axis=1)
    S = (np.abs(U) ** 2).mean(0) / np.sum(win ** 2)
    q = 2 * np.pi * np.arange(U.shape[1]) / L
    return q[1:], S[1:]


def height_diff(du, rs, w=None):
    """B(r) = <[du(s+r)-du(s)]^2>, optionally with pair weights w (frames x L, 0/1)."""
    B = np.empty(len(rs))
    for j, r in enumerate(rs):
        d2 = (du[:, r:] - du[:, :-r]) ** 2
        if w is None:
            B[j] = d2.mean()
        else:
            ww = w[:, r:] * w[:, :-r]
            B[j] = (d2 * ww).sum() / max(ww.sum(), 1)
    return B


def local_width(du, ls):
    out = []
    for l in ls:
        l = int(l)
        starts = np.arange(0, du.shape[1] - l + 1, max(l // 2, 1))
        x = np.arange(l)
        A = np.vstack([x, np.ones(l)]).T
        P = A @ np.linalg.pinv(A)                         # projector on lines
        vals = []
        for s0 in starts:
            seg = du[:, s0:s0 + l]
            res = seg - seg @ P.T
            vals.append((res ** 2).mean())
        out.append(np.mean(vals))
    return np.array(out)


def logbin_xy(x, y, per_decade=8):
    """Average y in logarithmic bins of x (geometric bin centres)."""
    e = np.logspace(np.log10(x.min()), np.log10(x.max() * 1.0001), int(np.log10(x.max() / x.min()) * per_decade) + 2)
    i = np.digitize(x, e)
    xs = np.array([np.exp(np.log(x[i == j]).mean()) for j in np.unique(i)])
    ys = np.array([y[i == j].mean() for j in np.unique(i)])
    return xs, ys


def zeta_S(q, S, rmin, rmax):
    qb, Sb = logbin_xy(q, S)
    return -(fit_power(qb, Sb, 2 * np.pi / rmax, 2 * np.pi / rmin)[0] + 1) / 2


def fit_power(x, y, lo, hi):
    sel = (x >= lo) & (x <= hi) & (y > 0)
    if sel.sum() < 3:
        return np.nan, np.nan
    p = np.polyfit(np.log(x[sel]), np.log(y[sel]), 1)
    return p[0], p[1]


def zeta_all(du, rmin, rmax, w=None):
    L = du.shape[1]
    q, S = structure_factor(du)
    rs = np.unique(np.round(np.logspace(0, np.log10(L // 2), 40)).astype(int))
    B = height_diff(du, rs)
    Bm = height_diff(du, rs, w) if w is not None else None
    ls = np.unique(np.round(np.logspace(np.log10(4), np.log10(L // 2), 25)).astype(int))
    W2 = local_width(du, ls)
    zS = zeta_S(q, S, rmin, rmax)
    zB = fit_power(rs, B, rmin, rmax)[0] / 2
    zBm = fit_power(rs, Bm, rmin, rmax)[0] / 2 if Bm is not None else np.nan
    zW = fit_power(ls, W2, rmin, rmax)[0] / 2
    return dict(q=q, S=S, r=rs, B=B, Bm=Bm, l=ls, W2=W2, zS=zS, zB=zB, zBm=zBm, zW=zW)


# ------------------------------------------------------------------ rotation onto the base plane
def rotation_sign(theta0, shape):
    """Sign of the scipy rotation angle that makes a wall x = a + y tan(theta0) vertical."""
    H, W = shape
    Y, X = np.mgrid[0:H, 0:W]
    m = (X < W / 2 + (Y - H / 2) * np.tan(np.radians(theta0))).astype(float)
    best = None
    for sgn in (1, -1):
        R = rotate(m, sgn * theta0, reshape=True, order=1, cval=-1)
        V = rotate(np.ones(shape), sgn * theta0, reshape=True, order=1, cval=0) > 0.999
        rows = np.nonzero(V.sum(1) > W / 2)[0]
        xs = [np.nonzero((R[r] >= 0.5) & V[r])[0].max() for r in rows if ((R[r] >= 0.5) & V[r]).any()]
        slope = abs(np.polyfit(np.arange(len(xs)), xs, 1)[0]) if len(xs) > 10 else np.inf
        if best is None or slope < best[0]:
            best = (slope, sgn)
    return best[1]


def heights_rotated(dom_area, dom_front, ang, valid):
    """Three heights per row of the rotated frame. valid: rotated mask of the image support."""
    Ra = rotate(dom_area.astype(float), ang, reshape=True, order=1, cval=0.0)
    Rf = rotate(dom_front.astype(float), ang, reshape=True, order=1, cval=0.0)
    n = Ra.shape[0]
    ua, uf, ub = (np.full(n, np.nan) for _ in range(3))
    gaps = np.zeros(n, int)
    for r in range(n):
        v = np.nonzero(valid[r])[0]
        if v.size < 10:
            continue
        c0, c1 = v[0], v[-1]
        a = np.clip(Ra[r, c0:c1 + 1], 0, 1)
        f = Rf[r, c0:c1 + 1]
        inside = f >= 0.5
        if not inside[0] or inside[-1] or not inside.any():
            continue                                   # wall not inside the valid part of the row
        last = np.nonzero(inside)[0][-1]
        first_out = np.nonzero(~inside)[0][0]
        # subpixel crossings by linear interpolation of f around 0.5
        fr = last + (f[last] - 0.5) / max(f[last] - f[last + 1], 1e-9)
        bk = first_out - 1 + (f[first_out - 1] - 0.5) / max(f[first_out - 1] - f[first_out], 1e-9)
        ua[r], uf[r], ub[r] = c0 + a.sum(), c0 + fr, c0 + bk
        seg = ~inside[first_out:last + 1]
        if seg.any():                                  # gaps (runs of outside) of length >= 2
            d = np.diff(np.concatenate([[0], seg.astype(np.int8), [0]]))
            runs = np.nonzero(d == -1)[0] - np.nonzero(d == 1)[0]
            gaps[r] = np.sum(runs >= 2)
    return ua, uf, ub, gaps


# ------------------------------------------------------------------ local width with local rotation (PCA)
def contour_points(domf, step=0.5, border=2):
    """Subpixel wall contour of a filled domain, resampled at uniform arc length; (y, x) arrays."""
    from skimage.measure import find_contours
    cs = find_contours(domf.astype(float), 0.5)
    if not cs:
        return np.empty(0), np.empty(0)
    c = max(cs, key=len)
    seg = np.hypot(*np.diff(c, axis=0).T)
    arc = np.concatenate([[0], np.cumsum(seg)])
    a = np.arange(0, arc[-1], step)
    y, x = np.interp(a, arc, c[:, 0]), np.interp(a, arc, c[:, 1])
    H, W = domf.shape
    keep = (y > border) & (y < H - 1 - border) & (x > border) & (x < W - 1 - border)
    return y[keep], x[keep]


def pca_width(y, x, theta0_deg, ells, step=0.5):
    """For windows of projected length l along the base direction (tilt theta0), w2 = smallest
    eigenvalue of the covariance of the contour points (= variance normal to the segment's own
    best line, i.e. after rotating the segment by its local tilt). Returns per-l lists of w2 and
    local tilt (deg, relative to theta0)."""
    t0 = np.radians(theta0_deg)
    s = y * np.cos(t0) + x * np.sin(t0)
    o = np.argsort(s)
    s, y, x = s[o], y[o], x[o]
    C = [np.concatenate([[0], np.cumsum(v)]) for v in (np.ones_like(y), y, x, y * y, x * x, x * y)]
    out_w2, out_tilt = [], []
    for l in ells:
        s0 = np.arange(s[0], s[-1] - l, l / 2)
        a = np.searchsorted(s, s0)
        b = np.searchsorted(s, s0 + l)
        n = b - a
        ok = n >= max(4, 0.5 * l / step)
        a, b, n = a[ok], b[ok], n[ok].astype(float)
        if n.size == 0:
            out_w2.append(np.empty(0)); out_tilt.append(np.empty(0)); continue
        my, mx = (C[1][b] - C[1][a]) / n, (C[2][b] - C[2][a]) / n
        syy = (C[3][b] - C[3][a]) / n - my ** 2
        sxx = (C[4][b] - C[4][a]) / n - mx ** 2
        sxy = (C[5][b] - C[5][a]) / n - mx * my
        lam = (syy + sxx) / 2 - np.sqrt(((syy - sxx) / 2) ** 2 + sxy ** 2)
        alpha = 0.5 * np.arctan2(2 * sxy, syy - sxx)          # major axis, from y towards x
        out_w2.append(np.clip(lam, 0, None)); out_tilt.append(np.degrees(alpha - t0))
    return out_w2, out_tilt


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
            w2l, _ = pca_width(yy_, xx_, theta, ells_t)
            zp.append(fit_power(ells_t.astype(float), np.array([v.mean() for v in w2l]), args.rmin, H * args.rmax_frac)[0] / 2)
        print(f"   local width with local rotation (PCA on the contour): zeta {np.mean(zp):.2f}±{np.std(zp):.2f}")
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
    ua_all, uf_all, ub_all, gp_all, contours = [], [], [], [], []
    edge = lambda lab_: np.setdiff1d(np.unique(lab_[:, 0]), [0])
    for i, kf in enumerate(frames):
        dom = ta <= kf
        lab_, _ = label(dom)
        dom = np.isin(lab_, edge(lab_))
        domf = binary_fill_holes(dom)
        ua, uf, ub, gp = heights_rotated(dom, domf, ang, valid)
        ua_all.append(ua); uf_all.append(uf); ub_all.append(ub); gp_all.append(gp)
        contours.append(contour_points(domf))
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
            if k in ("q", "S", "r", "B", "Bm", "l", "W2", "zS", "dzS", "zB", "zBm", "zW", "zS_thirds") and v is not None})

# ------------------------------------------------------------------ local width, local rotation
ells = np.unique(np.round(np.logspace(np.log10(4), np.log10(L), 28)).astype(int)).astype(float)
if args.from_h:
    pts = [(np.arange(L, dtype=float), row) for row in U["h"]]
    th_pts = 0.0
else:
    pts = contours
    th_pts = theta0
per_frame = [pca_width(yy_, xx_, th_pts, ells) for yy_, xx_ in pts]     # [(w2 lists, tilt lists)]
W2f = np.array([[v.mean() if v.size else np.nan for v in pf[0]] for pf in per_frame])   # frames x ells
NW = np.array([[v.size for v in pf[0]] for pf in per_frame])
# same, without the segments tilted by more than --max-tilt from the base plane (folds, tongues)
keepT = [[np.abs(t_) <= args.max_tilt for t_ in pf[1]] for pf in per_frame]
W2r = np.array([[v[k_].mean() if k_.any() else np.nan for v, k_ in zip(pf[0], kp)] for pf, kp in zip(per_frame, keepT)])
NR = np.array([[k_.sum() for k_ in kp] for kp in keepT])
w2_restr = np.nansum(np.nan_to_num(W2r) * NR, 0) / np.maximum(NR.sum(0), 1)
frac_excl = 1 - NR.sum(0) / np.maximum(NW.sum(0), 1)
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
z_shear = res[ref_name]["zW"]
z_restr = zfit(w2_restr)
print(f"local width w2(l), each segment rotated by its own tilt (PCA), {'contour' if not args.from_h else 'h'}: "
      f"zeta = {z_pca:.3f} ± {dz_pca:.3f}; same on u_{ref_name}: {z_ctrl:.3f}; "
      f"line removed without rotation (w(l) above): {z_shear:.3f}")
print(f"   without segments tilted > {args.max_tilt:g} deg from the base plane: zeta = {z_restr:.3f} "
      f"(segments excluded: {frac_excl[(ells >= args.rmin) & (ells <= rmax)].mean():.0%} in the fit range)")
for j in np.unique(np.searchsorted(ells, [8, 25, 85, 250, 850])).clip(0, len(ells) - 1):
    print(f"   l = {ells[j]:5.0f} px ({ells[j] * kk:6.1f} {unit}): w = {np.sqrt(w2_pca[j]) * kk:6.3f} {unit}, "
          f"local tilt sd {tilt_sd[j]:5.1f} deg, {int(NW[:, j].sum())} segments, "
          f"{frac_excl[j]:.0%} tilted > {args.max_tilt:g} deg")
np.savez(os.path.join(args.dir, "local_width.npz"), ells=ells, w2_pca=w2_pca, dw2=dw2, w2_ctrl=w2_ctrl,
         w2_restricted=w2_restr, frac_excluded=frac_excl, max_tilt=args.max_tilt, zeta_restricted=z_restr,
         tilt_sd=tilt_sd, n_segments=NW.sum(0), zeta_pca=z_pca, dzeta_pca=dz_pca, zeta_ctrl=z_ctrl,
         theta0=th_pts, every=args.every, rmin=args.rmin, rmax=rmax, um_per_px=k_um)

fig, ax = plt.subplots(1, 3, figsize=(18, 5))
yy_, xx_ = pts[len(pts) // 2]
lshow = ells[np.argmin(np.abs(ells * kk - (20 if k_um else 160)))]
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
    ax[0].plot(px, py, ".", ms=1.2, color=c)
    ax[0].plot(cx + tt * d[1], cy + tt * d[0], "-", color="k", lw=1)
ax[0].invert_yaxis(); ax[0].set_aspect("equal")
ax[0].set_title(f"Segments of {lshow * kk:.0f} {unit} and their own best lines (frame {frames[len(pts) // 2]})", fontsize=10)
ax[0].set_xlabel("x [px]" if not args.from_h else "h [px]"); ax[0].set_ylabel("y [px]" if not args.from_h else "column")
ax[1].loglog(ells * kk, w2_pca * kk ** 2, "o-", ms=3, color="C0", label=f"local rotation, {'contour' if not args.from_h else 'h'}: ζ = {z_pca:.2f} ± {dz_pca:.2f}")
ax[1].fill_between(ells * kk, (w2_pca - 2 * dw2) * kk ** 2, (w2_pca + 2 * dw2) * kk ** 2, color="C0", alpha=0.2)
ax[1].loglog(ells * kk, w2_restr * kk ** 2, "D-", ms=3, mfc="none", color="C4",
             label=f"local rotation, segments tilted ≤ {args.max_tilt:g}°: ζ = {z_restr:.2f}")
ax[1].loglog(ells * kk, w2_ctrl * kk ** 2, "s--", ms=3, mfc="none", color="C1", label=f"local rotation, u_{ref_name}: ζ = {z_ctrl:.2f}")
ax[1].loglog(res[ref_name]["l"] * kk, res[ref_name]["W2"] * kk ** 2, "^:", ms=3, color="C2", label=f"line removed, no rotation: ζ = {z_shear:.2f}")
ax[1].axvspan(args.rmin * kk, rmax * kk, color="0.92", zorder=0)
ax[1].set_xlabel(f"segment length ℓ [{unit}]"); ax[1].set_ylabel(f"⟨w²(ℓ)⟩ [{unit}²]"); ax[1].legend(fontsize=8)
ax[1].set_title(r"Local width, $\langle w^2\rangle \sim \ell^{2\zeta}$ (grey: fit range; band: ±2σ)", fontsize=10)
ax[2].semilogx(ells * kk, tilt_sd, "o-", ms=3, label="sd of the local tilt")
ax2b = ax[2].twinx(); ax2b.semilogx(ells * kk, 100 * frac_excl, "s--", ms=3, color="C3", mfc="none")
ax2b.set_ylabel(f"% of segments tilted > {args.max_tilt:g}°", color="C3")
ax[2].set_xlabel(f"segment length ℓ [{unit}]"); ax[2].set_ylabel("sd of the local tilt [deg]")
ax[2].set_title("Spread of the local tilt of the segments", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(args.dir, "fig_local_width.png"), dpi=120); plt.close(fig)

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
