"""Step 9: roughness exponent from wavelets (DOG-3, cross-checked with DOG-4) in local frames.

Usage: python step9_wavelet.py DIR [--every 1] [--min-kept 0.98] [--n-boot 300]
       python step9_wavelet.py --selftest              (synthetic walls with known zeta)
Input: DIR/contours.npz (TIFF studies, written by step7) or DIR/h_xt_sub.npy (video studies and
simulations: the curve (column, h) of each frame), and DIR/local_width.npz etc. for the comparison.
Outputs: DIR/wavelet.npz, DIR/fig_wavelet.png and a table on stdout.

Method (the estimator is the validated one of wall_exponents.py, imported from there):
  every curve is resampled at uniform arc length (1 px); for each wavelet scale a, windows of
  n = 2 H 2.2 ~ 26 a points (H = 6 a) are put in their own PCA frame (t along the window's
  principal axis, u normal to it); windows where t is not monotonic, or that do not cover t_c +- H,
  are discarded; u(t) is interpolated on a unit grid and W_m(a) = sum psi_m((t - t_c)/a) u(t) / a,
  psi_m = m-th derivative of a Gaussian (m vanishing moments, L1 normalisation), so that
  F(a) = sqrt(<W^2>) ~ a^zeta. The local frame and the vanishing moments make W insensitive to the
  local tilt and to the curvature of the wall up to order m - 1.
Scale selection: a scale is used only if (i) the kept fraction of windows is >= --min-kept, (ii) the
window arc length 26 a <= 0.6 R_c(a), with R_c(a) the 10th percentile of the radius of curvature of
the curve smoothed over the window half-length (13 a), the largest ratio of the validated "curved"
geometry of wall_exponents.py, and (iii) l ~ 4 a >= rmin of step7 (resolution). Fit range: the
longest run of usable scales whose local slopes d ln F / d ln a lie within +-0.15 of their median.
Errors: block bootstrap over contiguous segments of the wall (6 along the arc x 12 time blocks of
frames = 72 blocks); also over time blocks only, as a conservative check.
"""
import argparse
import json
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from scipy.ndimage import gaussian_filter, gaussian_filter1d

from wall_exponents import dog_kernel, lin_width, local_frame_windows, resample_arclength, wavelet_rms

ap = argparse.ArgumentParser()
ap.add_argument("dir", nargs="?", help="study folder")
ap.add_argument("--every", type=int, default=1, help="use every N-th stored curve")
ap.add_argument("--min-kept", type=float, default=0.98, help="minimum kept fraction of windows per scale")
ap.add_argument("--curv-ratio", type=float, default=0.6, help="max window arc length / R_c")
ap.add_argument("--slope-tol", type=float, default=0.15, help="tolerance of the local slopes in the fit range")
ap.add_argument("--n-boot", type=int, default=300, help="bootstrap resamples")
ap.add_argument("--selftest", action="store_true", help="synthetic walls with known zeta")
args = ap.parse_args()
MS = (3, 4)
SCALES = np.round(np.geomspace(2, 90, 23), 2)       # factor ~1.19 between scales
NSEG, NTB = 6, 12


# ------------------------------------------------------------------ wavelet coefficients per window
def window_geometry(a):
    H = int(np.ceil(6 * a))
    return H, int(np.ceil(2 * H * 2.2)), max(1, int(a))


def wavelet_windows(X, Y, a, ms=MS):
    """Same computation as wall_exponents.wavelet_rms, but returning every kept window's coefficient
    and the arc-length index of its centre. Returns {m: W}, centres, n_kept, n_windows."""
    H, n, stride = window_geometry(a)
    if len(X) < n + 1:
        return {m: np.empty(0) for m in ms}, np.empty(0, int), 0, 0
    T, U, vf = local_frame_windows(X, Y, n, stride)
    # the same monotonicity mask as local_frame_windows, to know which windows it kept
    Xw = sliding_window_view(X, n + 1)[::stride]
    Yw = sliding_window_view(Y, n + 1)[::stride]
    xc, yc = Xw - Xw.mean(1, keepdims=True), Yw - Yw.mean(1, keepdims=True)
    phi = 0.5 * np.arctan2(2 * (xc * yc).mean(1), (xc * xc).mean(1) - (yc * yc).mean(1))
    Tm = xc * np.cos(phi)[:, None] + yc * np.sin(phi)[:, None]
    Tm = Tm * np.where(Tm[:, -1] - Tm[:, 0] < 0, -1.0, 1.0)[:, None]
    valid = np.all(np.diff(Tm, axis=1) > 0, axis=1)
    assert valid.sum() == T.shape[0]
    starts = (np.arange(Xw.shape[0]) * stride)[valid]
    c = n // 2
    tc = T[:, c]
    ok = (T[:, 0] <= tc - H) & (T[:, -1] >= tc + H)
    T, U, tc, starts = T[ok], U[ok], tc[ok], starts[ok]
    M = T.shape[0]
    if M == 0:
        return {m: np.empty(0) for m in ms}, np.empty(0, int), 0, Xw.shape[0]
    off = (np.arange(M) * (4.0 * n + 10))[:, None]
    Tf, Uf = (T + off).ravel(), U.ravel()
    Gf = (tc[:, None] + np.arange(-H, H + 1)[None, :] + off).ravel()
    idx = np.clip(np.searchsorted(Tf, Gf), 1, Tf.size - 1)
    x0, x1 = Tf[idx - 1], Tf[idx]
    w = (Gf - x0) / (x1 - x0)
    val = (Uf[idx - 1] * (1 - w) + Uf[idx] * w).reshape(M, -1)
    return {m: val @ dog_kernel(a, m, H) for m in ms}, starts + c, M, Xw.shape[0]


def radius_p10(X, Y, a):
    """10th percentile of the radius of curvature of the curve smoothed over 13 a points."""
    sig = 13 * a
    if len(X) < 6 * sig:
        return np.nan
    xs, ys = gaussian_filter1d(X, sig, mode="nearest"), gaussian_filter1d(Y, sig, mode="nearest")
    dx, dy = np.gradient(xs), np.gradient(ys)
    ddx, ddy = np.gradient(dx), np.gradient(dy)
    k = np.abs(dx * ddy - dy * ddx) / np.maximum((dx * dx + dy * dy) ** 1.5, 1e-12)
    m = int(2 * sig)
    k = k[m:-m]
    return np.percentile(1 / np.maximum(k, 1e-9), 10) if k.size else np.nan


# ------------------------------------------------------------------ analysis of a set of curves
def analyse(curves, scales):
    """curves: list of (X, Y) at unit arc-length spacing, in time order. Returns per-scale block sums
    of W^2 (m), counts, kept fractions and R_c."""
    nf = len(curves)
    S2 = {m: np.zeros((len(scales), NTB, NSEG)) for m in MS}
    CNT = np.zeros((len(scales), NTB, NSEG))
    kept_n, tot_n = np.zeros(len(scales)), np.zeros(len(scales))
    Rc = np.full(len(scales), np.nan)
    rc_samples = [[] for _ in scales]
    for f, (X, Y) in enumerate(curves):
        tb = f * NTB // nf
        for j, a in enumerate(scales):
            Wm, cen, M, Mall = wavelet_windows(X, Y, a)
            kept_n[j] += M
            tot_n[j] += Mall
            if M:
                seg = np.minimum(cen * NSEG // len(X), NSEG - 1)
                for m in MS:
                    S2[m][j, tb] += np.bincount(seg, Wm[m] ** 2, minlength=NSEG)
                CNT[j, tb] += np.bincount(seg, minlength=NSEG)
            if f % max(1, nf // 24) == 0:                    # curvature on a subset of frames
                rc_samples[j].append(radius_p10(X, Y, a))
    for j in range(len(scales)):
        v = np.array(rc_samples[j], float)
        Rc[j] = np.nanmedian(v) if np.isfinite(v).any() else np.nan
    kept = np.where(tot_n > 0, kept_n / np.maximum(tot_n, 1), 0.0)
    # R_c cannot be measured when the curve is shorter than ~6 x 13 a; smoothing over more points only
    # straightens the curve, so the last measured R_c is carried forward as a lower bound
    Rc_ext = np.isnan(Rc) & (np.cumsum(np.isfinite(Rc)) > 0)
    Rc = np.fmax.accumulate(np.where(np.isfinite(Rc), Rc, -np.inf))
    Rc = np.where(np.isfinite(Rc), Rc, np.nan)
    return dict(S2=S2, CNT=CNT, kept=kept, Rc=Rc, Rc_ext=Rc_ext, n_windows=kept_n)


def F_of(R, m, pick=None):
    """F(a) from the block sums; pick: list of (time block, segment) pairs (bootstrap)."""
    if pick is None:
        s, c = R["S2"][m].sum((1, 2)), R["CNT"].sum((1, 2))
    else:
        s = R["S2"][m][:, pick[:, 0], pick[:, 1]].sum(1)
        c = R["CNT"][:, pick[:, 0], pick[:, 1]].sum(1)
    return np.where(c > 0, np.sqrt(s / np.maximum(c, 1)), np.nan)


def choose_range(a, F, usable, tol):
    """Longest run of consecutive usable scales whose local slopes are within tol of their median."""
    idx = np.nonzero(usable & np.isfinite(F) & (F > 0))[0]
    runs = np.split(idx, np.nonzero(np.diff(idx) != 1)[0] + 1) if idx.size else []
    best = None
    for r in runs:
        for i in range(len(r)):
            for j in range(i + 2, len(r)):
                sub = r[i:j + 1]
                sl = np.diff(np.log(F[sub])) / np.diff(np.log(a[sub]))
                if np.max(np.abs(sl - np.median(sl))) <= tol:
                    key = (len(sub), -np.std(sl))
                    if best is None or key > best[0]:
                        best = (key, sub)
    return None if best is None else best[1]


def zfit(a, F, sel):
    return np.polyfit(np.log(a[sel]), np.log(F[sel]), 1)[0]


def fit_all(R, scales, rmin, n_boot, rng):
    a = np.asarray(scales, float)
    usable = (R["kept"] >= args.min_kept) & (26 * a <= args.curv_ratio * R["Rc"]) & (4 * a >= rmin)
    out = dict(usable=usable)
    for m in MS:
        F = F_of(R, m)
        rng_m = choose_range(a, F, usable, args.slope_tol)
        out[m] = dict(F=F, range=rng_m)
        if rng_m is None:
            out[m].update(z=np.nan, dz=np.nan, dz_t=np.nan)
            continue
        z = zfit(a, F, rng_m)
        blocks = np.array([(t, s) for t in range(NTB) for s in range(NSEG)])
        bz = [zfit(a, F_of(R, m, blocks[rng.integers(0, len(blocks), len(blocks))]), rng_m) for _ in range(n_boot)]
        bt = []
        for _ in range(n_boot):
            tt = rng.integers(0, NTB, NTB)
            bt.append(zfit(a, F_of(R, m, np.array([(t, s) for t in tt for s in range(NSEG)])), rng_m))
        # a second plateau (disjoint from the first) if there is one
        rest = usable.copy()
        rest[max(rng_m[0] - 1, 0):rng_m[-1] + 2] = False
        rng2 = choose_range(a, F, rest, args.slope_tol)
        out[m].update(z=z, dz=float(np.nanstd(bz)), dz_t=float(np.nanstd(bt)), range2=rng2,
                      z2=zfit(a, F, rng2) if rng2 is not None else np.nan)
    return out


# ------------------------------------------------------------------ self test
if args.selftest:
    from wall_tools import contour_points
    rng = np.random.default_rng(3)
    H_, W_, theta, nw = 1006, 1197, -12.0, 12
    Y_, X_ = np.mgrid[0:H_, 0:W_]
    print(f"selftest: {nw} synthetic walls per case, {H_} rows, tilted {theta} deg, self-affine below 256 px, slope 6 % "
          f"at 64 px; contour -> unit arc length -> DOG wavelets in local frames; automatic range (rmin = 16 px)")
    cases = [(z, None, b) for z in (0.5, 1.0, 1.25) for b in (0, 3)] + [(None, 150, 0), (None, 150, 3)]
    for zeta, l0, blur in cases:
        curves = []
        for _ in range(nw):
            k = np.fft.rfftfreq(H_)
            amp = np.zeros_like(k)
            if l0 is None:                                  # self-affine below xi = 256 px, flat above
                amp[1:] = np.maximum(k[1:], 1 / 256) ** (-(1 + 2 * zeta) / 2)
            else:
                amp[1:] = np.sqrt(k[1:] ** (-(1 + 2 * 0.5)) * (1 + (k[1:] * l0) ** 2) ** (-(1.25 - 0.5)))
            u = np.fft.irfft(amp * np.exp(2j * np.pi * rng.random(k.size)), n=H_)
            u *= 0.06 * 64 / lin_width(u, 64)          # slope 6 % at l = 64 px, as w(100 px) ~ 9 px in Study 2
            dom = (X_ < 600 + (Y_ - H_ / 2) * np.tan(np.radians(theta)) + u[:, None]).astype(float)
            # contour of the domain smoothed by sigma = 1 px (no pixel staircase), as step7 saves it;
            # "blur" adds an optical PSF on top
            yy, xx = contour_points(gaussian_filter(dom, np.hypot(1, blur)))
            curves.append(resample_arclength(xx, yy))
        if zeta == 0.5 and blur == 0:                       # the wrapper reproduces wavelet_rms
            X, Y = curves[0]
            for a in (3, 8):
                Wm, _, M, _ = wavelet_windows(X, Y, a)
                if M == 0:                                  # wavelet_rms cannot handle zero kept windows
                    print(f"   check a = {a}: no kept windows")
                    continue
                r, _ = wavelet_rms(X, Y, a, MS)
                print(f"   check a = {a}: wavelet_rms DOG-3 {r[3]:.6g} vs wrapper {np.sqrt(np.mean(Wm[3] ** 2)):.6g}; "
                      f"DOG-4 {r[4]:.6g} vs {np.sqrt(np.mean(Wm[4] ** 2)):.6g}")
        sc = SCALES[np.array([window_geometry(a_)[1] for a_ in SCALES]) < 0.9 * np.median([len(c[0]) for c in curves])]
        R = analyse(curves, sc)
        res = fit_all(R, sc, 16, 100, rng)
        txt = "  ".join(f"DOG-{m}: {res[m]['z']:.2f}±{res[m]['dz']:.2f} (a {sc[res[m]['range'][0]]:.0f}-{sc[res[m]['range'][-1]]:.0f})"
                        if res[m]["range"] is not None else f"DOG-{m}: no range" for m in MS)
        lab = f"zeta = {zeta}" if l0 is None else f"two regimes, l0 = {l0} px (1.25 -> 0.5)"
        fixed = "  | fixed a 4-8: " + ", ".join(f"{zfit(sc, res[m]['F'], (sc >= 4) & (sc <= 8)):.2f}" for m in MS) + \
                "; a 8-19: " + ", ".join(f"{zfit(sc, res[m]['F'], (sc >= 7.9) & (sc <= 19)):.2f}" for m in MS) + \
                f"; kept at a = 8, 19: {R['kept'][np.argmin(abs(sc - 8))]:.3f}, {R['kept'][np.argmin(abs(sc - 19))]:.3f}"
        print(f"   {lab}, blur {blur} px: {txt}{fixed}")
    raise SystemExit


# ------------------------------------------------------------------ data
d = args.dir
with open(os.path.join(d, "meta.json")) as f:
    meta = json.load(f)
k_um = meta.get("um_per_px") or 0.0
kk, unit = (k_um, "µm") if k_um else (1.0, "px")
lw = np.load(os.path.join(d, "local_width.npz"))
rmin = float(lw["rmin"])
curves, kind = [], ""
if os.path.exists(os.path.join(d, "contours.npz")):
    c = np.load(os.path.join(d, "contours.npz"))
    off = c["offsets"]
    for i in range(0, len(off) - 1, args.every):
        x, y = c["x"][off[i]:off[i + 1]], c["y"][off[i]:off[i + 1]]
        curves.append(resample_arclength(x, y))
    kind = f"subpixel contours of {len(curves)} frames (step7)"
else:
    h = np.load(os.path.join(d, "h_xt_sub.npy"))[::args.every]
    for row in h:
        ok = np.isfinite(row)
        curves.append(resample_arclength(np.nonzero(ok)[0].astype(float), row[ok]))
    kind = f"curves (column, h) of {len(curves)} frames (h_xt_sub.npy)"
lens = np.array([len(c_[0]) for c_ in curves])
scales = SCALES[np.array([window_geometry(a)[1] for a in SCALES]) < 0.9 * np.median(lens)]
print(f"{d}: {kind}, arc length median {np.median(lens) * kk:.0f} {unit}; scales a = {scales[0]:g}-{scales[-1]:g} px")
R = analyse(curves, scales)
rng = np.random.default_rng(9)
res = fit_all(R, scales, rmin, args.n_boot, rng)
a = scales.astype(float)

# comparison with the earlier estimators over the same l ~ 4 a range
cmp = {}
rr = res[3]["range"]
if rr is not None:
    lo, hi = 4 * a[rr[0]], 4 * a[rr[-1]]
    ells = lw["ells"]
    sel = (ells >= lo) & (ells <= hi)
    for key, lab in (("w2_pca", "local width, contour"), ("w2_unfolded", "local width, folds excluded"),
                     ("w2_ctrl", "local width, u_area / h")):
        if key in lw.files and sel.sum() >= 3:
            v = lw[key]
            ok = sel & np.isfinite(v) & (v > 0)
            cmp[lab] = np.polyfit(np.log(ells[ok]), np.log(v[ok]), 1)[0] / 2
    rg = np.load(os.path.join(d, "roughness.npz"))
    pre = "area" if "area_q" in rg.files else "h"
    q, S = rg[f"{pre}_q"], rg[f"{pre}_S"]
    ok = (q >= 2 * np.pi / hi) & (q <= 2 * np.pi / lo)
    if ok.sum() >= 3:
        cmp["global S(q), same scales"] = -(np.polyfit(np.log(q[ok]), np.log(S[ok]), 1)[0] + 1) / 2
    cmp["global S(q), step7 fit range"] = float(rg[f"{pre}_zS"])
else:
    # no usable wavelet range: show the earlier estimators over their own fit ranges
    rg = np.load(os.path.join(d, "roughness.npz"))
    pre = "area" if "area_q" in rg.files else "h"
    for key, lab in (("zeta_pca", "local width, contour (own range)"), ("zeta_unfolded", "local width, folds excluded (own range)"),
                     ("zeta_ctrl", "local width, u_area / h (own range)")):
        if key in lw.files:
            cmp[lab] = float(lw[key])
    cmp["global S(q), step7 fit range"] = float(rg[f"{pre}_zS"])
if os.path.exists(os.path.join(d, "crossover.npz")):
    cr = np.load(os.path.join(d, "crossover.npz"), allow_pickle=True)
    cmp["l0 (step8, a / b / c) [µm]"] = tuple(float(cr[f"l0_{e}"]) for e in "abc" if f"l0_{e}" in cr.files)

print(f"   {'a[px]':>6s} {'l~4a':>8s} {'kept':>6s} {'26a/Rc':>7s} {'usable':>6s} {'F DOG-3':>10s} {'F DOG-4':>10s} "
      f"{'slope3':>7s} {'slope4':>7s}")
sl = {m: np.r_[np.nan, np.diff(np.log(res[m]["F"])) / np.diff(np.log(a))] for m in MS}
for j in range(len(a)):
    print(f"   {a[j]:6.2f} {4 * a[j] * kk:7.1f}{unit[0] if unit == 'px' else ''} {R['kept'][j]:6.3f} "
          f"{26 * a[j] / R['Rc'][j]:6.2f}{'*' if R['Rc_ext'][j] else ' '} {'yes' if res['usable'][j] else 'no':>6s} {res[3]['F'][j] * kk:10.4g} "
          f"{res[4]['F'][j] * kk:10.4g} {sl[3][j]:7.2f} {sl[4][j]:7.2f}")
for m in MS:
    r_ = res[m]
    if r_["range"] is None:
        print(f"   DOG-{m}: no usable range of scales")
        continue
    a0, a1 = a[r_["range"][0]], a[r_["range"][-1]]
    s2 = ""
    if r_.get("range2") is not None:
        b0, b1 = a[r_["range2"][0]], a[r_["range2"][-1]]
        s2 = f"; second plateau a = {b0:g}-{b1:g} px (l ~ {4 * b0 * kk:.1f}-{4 * b1 * kk:.1f} {unit}): zeta = {r_['z2']:.2f}"
    print(f"   DOG-{m}: zeta = {r_['z']:.3f} ± {r_['dz']:.3f} (segments x time blocks; time blocks only ± {r_['dz_t']:.3f}) "
          f"over a = {a0:g}-{a1:g} px (l ~ 4a = {4 * a0 * kk:.1f}-{4 * a1 * kk:.1f} {unit}){s2}")
if R["Rc_ext"].any():
    print("   * R_c not measurable at this scale (curve too short): lower bound from the largest measured scale")
for k_, v in cmp.items():
    print(f"   compare: {k_}: " + (", ".join(f"{x:.1f}" for x in v) if isinstance(v, tuple) else f"{v:.3f}"))
np.savez(os.path.join(d, "wavelet.npz"), a=a, ell=4 * a * kk, kept=R["kept"], Rc=R["Rc"], usable=res["usable"],
         n_windows=R["n_windows"], **{f"F_dog{m}": res[m]["F"] * kk for m in MS},
         **{f"zeta_dog{m}": res[m]["z"] for m in MS}, **{f"dzeta_dog{m}": res[m]["dz"] for m in MS},
         **{f"dzeta_time_dog{m}": res[m]["dz_t"] for m in MS},
         **{f"range_dog{m}": (a[res[m]["range"][0]], a[res[m]["range"][-1]]) if res[m]["range"] is not None else (np.nan, np.nan)
            for m in MS}, compare=json.dumps({k_: v for k_, v in cmp.items()}), um_per_px=k_um, rmin=rmin)

# ------------------------------------------------------------------ figure
fig, ax = plt.subplots(1, 4, figsize=(22, 5))
for m, col in zip(MS, ("C0", "C3")):
    r_ = res[m]
    ax[0].loglog(a, r_["F"] * kk, "o-", ms=3, color=col, label=f"DOG-{m}")
    if r_["range"] is not None:
        s_ = r_["range"]
        p = np.polyfit(np.log(a[s_]), np.log(r_["F"][s_] * kk), 1)
        ax[0].loglog(a[s_], np.exp(np.polyval(p, np.log(a[s_]))) * 1.6, "-", color=col, lw=2,
                     label=f"ζ = {r_['z']:.2f} ± {r_['dz']:.2f}")
ax[0].set_xlabel("wavelet scale a [px]"); ax[0].set_ylabel(f"F(a) = √⟨W²⟩ [{unit}]")
ax[0].set_title("Wavelet fluctuation function in local frames\n(fit line shifted up ×1.6)", fontsize=10)
ax[0].legend(fontsize=8)
sx = ax[0].secondary_xaxis("top", functions=(lambda x: 4 * x * kk, lambda x: x / (4 * kk)))
sx.set_xlabel(f"ℓ ≈ 4a [{unit}]")
for m, col in zip(MS, ("C0", "C3")):
    ax[1].semilogx(a, sl[m], "o-", ms=3, color=col, label=f"DOG-{m}")
    if res[m]["range"] is not None:
        ax[1].axvspan(a[res[m]["range"][0]], a[res[m]["range"][-1]], color=col, alpha=0.1)
for zz in (0.5, 2 / 3, 1.25):
    ax[1].axhline(zz, color="0.8", lw=0.8, zorder=0)
ax[1].plot(a[~res["usable"]], np.full((~res["usable"]).sum(), 0.3), "x", color="0.5", label="scale not usable")
ax[1].set_xlabel("a [px]"); ax[1].set_ylabel("d ln F / d ln a")
ax[1].set_title("Local slopes (shaded: fit range)", fontsize=10); ax[1].legend(fontsize=8)
ax[2].semilogx(a, R["kept"], "o-", ms=3, color="k", label="kept fraction of windows")
ax[2].axhline(args.min_kept, color="k", ls=":", lw=0.8)
ax2 = ax[2].twinx()
ax2.semilogx(a, 26 * a / R["Rc"], "s--", ms=3, color="C2", label="26a / R_c")
ax2.axhline(args.curv_ratio, color="C2", ls=":", lw=0.8)
ax2.set_ylabel("window length / radius of curvature", color="C2")
ax[2].axvline(rmin / 4, color="0.5", ls="--", lw=0.8)
ax[2].set_xlabel("a [px]"); ax[2].set_ylabel("kept fraction")
ax[2].set_title(f"Scale selection: kept ≥ {args.min_kept:g}, 26a ≤ {args.curv_ratio:g} R_c, 4a ≥ rmin (dashed)", fontsize=10)
ax[2].legend(loc="lower left", fontsize=8); ax2.legend(loc="lower right", fontsize=8)
labs = [f"DOG-{m}" + ("" if np.isfinite(res[m]["z"]) else " (none)") for m in MS] + [k_ for k_, v in cmp.items() if not isinstance(v, tuple)]
vals = [res[m]["z"] for m in MS] + [v for v in cmp.values() if not isinstance(v, tuple)]
errs = [res[m]["dz"] for m in MS] + [0] * (len(vals) - len(MS))
ax[3].errorbar(vals, range(len(vals)), xerr=errs, fmt="o", color="k", capsize=3)
ax[3].set_yticks(range(len(vals))); ax[3].set_yticklabels(labs, fontsize=8)
for zz in (0.5, 2 / 3, 1.25):
    ax[3].axvline(zz, color="0.8", lw=0.8, zorder=0)
if rr is None:
    ax[3].text(0.5, 0.97, "no usable wavelet scale\n(kept < " + f"{args.min_kept:g}" + " or 26a > "
               + f"{args.curv_ratio:g}" + " R_c at every a)\n→ no wavelet ζ; earlier estimators\nover their own fit ranges",
               transform=ax[3].transAxes, ha="center", va="top", fontsize=8, color="C3")
    ax[3].set_ylim(-0.5, len(vals) + 1.6)
ax[3].set_xlabel("ζ")
ax[3].set_title("ζ: wavelets vs earlier estimators" + (" (same ℓ ≈ 4a range)" if rr is not None else ""), fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(d, "fig_wavelet.png"), dpi=110); plt.close(fig)
print(f"   outputs: {d}/wavelet.npz, fig_wavelet.png")
