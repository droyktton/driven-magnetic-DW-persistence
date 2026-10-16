"""Contour-based roughness, overhangs and KPZ-term test from an arrival map (frame of switching per pixel, -1 never).
Works for a bubble (--centre, optional excluded sector) or for walls crossing the field (ROI mask).
Usage: python contour_analysis.py ARRIVAL.npz OUTDIR --frames K0 K1 [--um-per-px s] [--centre y x] [--exclude-deg a b]
       [--roi-erode px] [--ddef 8] [--every n] [--n-boot 300]
ARRIVAL.npz from arrival_map.py. ROI = bright aperture of its mean frame, eroded by --roi-erode px (0: whole image).
Per frame k: subpixel contour of {0 <= ta <= k} (1 px smoothing), points outside the ROI, near defects (never-switched
islands) or inside the excluded sector dropped; contiguous pieces >= 60 px kept.
 - local width w2(E) of arc windows rotated to their own axis (wall_tools.overhang_segments), minus the variance of a
   perfect arc of radius R for a bubble; zeta from E in [8, 80] px; P_ov(E) at tol 2/3/5 px.
 - KPZ term (as step10 B): segments of arc length L; pixels switching in (k, k+dt] assigned to the nearest contour
   point; v_n = count / (L dt); dphi = segment normal - normal of the 8L window around it; v/<v>_k = a[1 + b dphi +
   kappa/2 dphi^2], Lambda = 1 + kappa (0: elastic qEW, 1: isotropic, >> 1: qKPZ); 12 time blocks bootstrap.
 - bubble only: macroscopic v_n(phi) (radial advance between frames, Fourier fit n <= 4 on the allowed angles),
   Lambda(phi) = (v + v'')/v; skewness of r(theta) - low modes on the allowed arc."""
import argparse, os, sys, json
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.ndimage import binary_fill_holes, label, gaussian_filter, distance_transform_edt, binary_erosion
from scipy.spatial import cKDTree
from scipy.stats import skew

from wall_tools import overhang_segments, fit_power
from skimage.measure import find_contours

ap = argparse.ArgumentParser(); ap.add_argument("arrival"); ap.add_argument("outdir")
ap.add_argument("--frames", type=int, nargs=2, required=True); ap.add_argument("--um-per-px", type=float, default=np.nan)
ap.add_argument("--centre", type=float, nargs=2); ap.add_argument("--exclude-deg", type=float, nargs=2)
ap.add_argument("--roi-erode", type=int, default=0); ap.add_argument("--ddef", type=float, default=8)
ap.add_argument("--skip-kpz", action="store_true"); ap.add_argument("--n-boot", type=int, default=300); ap.add_argument("--every", type=int, default=1)
a = ap.parse_args(); os.makedirs(a.outdir, exist_ok=True); rng = np.random.default_rng(3)
A = np.load(a.arrival); ta = A["ta"].astype(int); H, W = ta.shape
k0, k1 = a.frames; bubble = a.centre is not None
# ROI: bright aperture of the movie, eroded
roi = np.ones((H, W), bool)
if a.roi_erode:
    mean = A["mean_frame"]
    roi = mean > 0.4 * np.median(mean[H // 4:3 * H // 4, W // 4:3 * W // 4])
    roi = binary_fill_holes(roi); l, _ = label(roi); s = np.bincount(l.ravel()); s[0] = 0; roi = l == s.argmax()
    roi = binary_erosion(roi, iterations=a.roi_erode)
roi[:3] = roi[-3:] = False; roi[:, :3] = roi[:, -3:] = False
# defects: never-switched islands (components of ta < 0 not touching the ROI border / image border, < 20000 px)
never = ta < 0; l, n = label(never); sz = np.bincount(l.ravel())
touch = np.unique(np.concatenate([l[~roi]]))
defect = np.isin(l, [i for i in range(1, n + 1) if sz[i] < 20000 and i not in touch])
ddef = distance_transform_edt(~defect) if defect.any() else np.full((H, W), np.inf)
print(f"ROI {roi.mean():.1%} of the image; defect islands {defect.sum()} px")
if bubble:
    y0, x0 = a.centre

def allowed(y, x):
    yi, xi = np.clip(np.round(y).astype(int), 0, H - 1), np.clip(np.round(x).astype(int), 0, W - 1)
    ok = roi[yi, xi] & (ddef[yi, xi] > a.ddef)
    if bubble and a.exclude_deg:
        th = np.degrees(np.arctan2(y - y0, x - x0)); lo, hi = a.exclude_deg
        ok &= ~((th >= lo) & (th <= hi))
    return ok

def pieces(k, step=0.5, minlen=60):
    dom = (ta >= 0) & (ta <= k)
    if bubble:
        dom = binary_fill_holes(dom); l, _ = label(dom); dom = l == l[int(y0), int(x0)]
    out = []
    for c in find_contours(gaussian_filter(dom.astype(float), 1.0), 0.5):
        seg = np.hypot(*np.diff(c, axis=0).T); arc = np.concatenate([[0], np.cumsum(seg)])
        if arc[-1] < minlen: continue
        s = np.arange(0, arc[-1], step); y, x = np.interp(s, arc, c[:, 0]), np.interp(s, arc, c[:, 1])
        ok = allowed(y, x)
        # split into runs of allowed points (closed contours: join the wrap-around run)
        idx = np.nonzero(ok)[0]
        if idx.size == 0: continue
        br = np.nonzero(np.diff(idx) > 1)[0]
        runs = np.split(idx, br + 1)
        if len(runs) > 1 and runs[0][0] == 0 and runs[-1][-1] == ok.size - 1 and np.hypot(y[0] - y[-1], x[0] - x[-1]) < 2:
            runs[0] = np.concatenate([runs[-1], runs[0]]); runs = runs[:-1]
        for r in runs:
            if r.size * step >= minlen: out.append((y[r], x[r]))
    return out

frames = list(range(k0, k1 + 1, a.every)); k1 = frames[-1]
C = {k: pieces(k) for k in frames}
Rk = {}
if bubble:
    for k in frames:
        ys = np.concatenate([p[0] for p in C[k]]); xs = np.concatenate([p[1] for p in C[k]])
        Rk[k] = np.hypot(ys - y0, xs - x0).mean()
print(f"frames {k0}-{k1}: contour length per frame {np.mean([sum(p[0].size for p in C[k]) * 0.5 for k in frames]):.0f} px "
      f"in {np.mean([len(C[k]) for k in frames]):.1f} pieces" + (f"; mean radius {Rk[k0]:.0f} -> {Rk[k1]:.0f} px" if bubble else ""))

# ---------------------------------------------------------------- local width and overhangs
TOLS = [2.0, 3.0, 5.0]
lmax = 0.95 * min(max(p[0].size * 0.5 for p in C[k]) for k in frames)
if bubble: lmax = min(lmax, min(Rk.values()))
lams = np.unique(np.round(np.logspace(np.log10(4), np.log10(lmax), 24)))
NT = np.zeros(lams.size); NO = {t: np.zeros(lams.size) for t in TOLS}; EX = [[] for _ in lams]; W2 = [[] for _ in lams]
for k in frames:
    for y, x in C[k]:
        nt, no, ex, w2 = overhang_segments(y, x, lams, 0.5, TOLS)
        NT += nt
        for t in TOLS: NO[t] += no[t]
        for j, lam in enumerate(lams):
            circ = 0.0
            if bubble:
                al = np.linspace(-lam / (2 * Rk[k]), lam / (2 * Rk[k]), 200); circ = np.var(Rk[k] * np.cos(al))
            EX[j].append(ex[j]); W2[j].append(np.asarray(w2[j]) - circ)
E = np.array([np.median(np.concatenate(e)) if sum(len(v) for v in e) else np.nan for e in EX])
w2m = np.array([np.mean(np.concatenate(w)) if sum(len(v) for v in w) else np.nan for w in W2])
P = {t: NO[t] / np.maximum(NT, 1) for t in TOLS}
ok = (E >= 8) & (E <= 80) & (w2m > 0)
zW = fit_power(E[ok], w2m[ok], 8, 80)[0] / 2
ok2 = (E >= 20) & (E <= min(80, E[np.isfinite(E)].max())) & (w2m > 0)
zW2 = fit_power(E[ok2], w2m[ok2], 20, 80)[0] / 2 if ok2.sum() > 2 else np.nan
def half(xv, pv):
    i = np.nonzero(pv >= 0.5)[0]
    if i.size == 0: return np.nan
    j = i[0]
    if j == 0: return xv[0]
    f = (0.5 - pv[j - 1]) / (pv[j] - pv[j - 1]); return float(np.exp(np.log(xv[j - 1]) + f * (np.log(xv[j]) - np.log(xv[j - 1]))))
l0g = {t: half(E, P[t]) for t in TOLS}
s = a.um_per_px
print(f"local width: zeta = {zW:.2f} (8-80 px), {zW2:.2f} (20-80 px); w at E = 10/20/50 px: " +
      "/".join(f"{np.sqrt(np.interp(np.log(e), np.log(E[w2m > 0]), w2m[w2m > 0])):.2f}" for e in (10, 20, 50)) + " px")
print("P_ov (tol 3) at E = 10/20/50/100 px: " + "/".join(f"{np.interp(np.log(e), np.log(E), P[3.0]):.3f}" if e <= np.nanmax(E) else "-" for e in (10, 20, 50, 100)) +
      "; l0_geom: " + ", ".join(f"tol {t:g}: " + (f"{v:.0f} px ({v * s:.1f} um)" if np.isfinite(v) else f"> {np.nanmax(E):.0f} px") for t, v in l0g.items()))

# ---------------------------------------------------------------- KPZ term (step10 B on the contour)
def kpz(L, dt, step=0.5):
    rows = []
    n = int(round(L / step)); nb = 8 * n
    for k in frames:
        if k + dt > k1 + 0: continue
        sw = (ta > k) & (ta <= k + dt)
        for y, x in C[k]:
            if y.size < nb + n: continue
            tree = cKDTree(np.c_[y, x]); py, px = np.nonzero(sw)
            if py.size == 0: continue
            dist, near = tree.query(np.c_[py, px], distance_upper_bound=40); near = near[np.isfinite(dist)]
            cnt = np.bincount(near, minlength=y.size)
            for c0 in range(nb // 2, y.size - nb // 2 - n, n):
                sl = slice(c0, c0 + n); bl = slice(c0 + n // 2 - nb // 2, c0 + n // 2 + nb // 2)
                Y, X = y[sl], x[sl]; my, mx = Y.mean(), X.mean()
                cyy, cxx, cxy = ((Y - my) ** 2).mean(), ((X - mx) ** 2).mean(), ((Y - my) * (X - mx)).mean()
                ang = 0.5 * np.arctan2(2 * cxy, cxx - cyy)                  # tangent angle (x towards y)
                Yb, Xb = y[bl], x[bl]; mby, mbx = Yb.mean(), Xb.mean()
                byy, bxx, bxy = ((Yb - mby) ** 2).mean(), ((Xb - mbx) ** 2).mean(), ((Yb - mby) * (Xb - mbx)).mean()
                angb = 0.5 * np.arctan2(2 * bxy, bxx - byy)
                d = np.angle(np.exp(2j * (ang - angb))) / 2                # tilt relative to the larger-scale wall, mod pi
                rows.append((k, d, cnt[sl].sum() / (n * step) / dt))
    r = np.array(rows)
    if r.size == 0 or len(r) < 200: return None
    # normalise by the mean velocity of the frame
    for k in np.unique(r[:, 0]):
        m = r[:, 0] == k; r[m, 2] /= max(r[m, 2].mean(), 1e-12)
    keep = np.abs(r[:, 1]) < np.radians(60); r = r[keep]
    blocks = np.array_split(np.argsort(r[:, 0], kind="stable"), 12)
    def one(idx):
        d, v = r[idx, 1], r[idx, 2] / r[idx, 2].mean(); A = np.vstack([np.ones_like(d), d, d * d]).T
        c = np.linalg.lstsq(A, v, rcond=None)[0]; return 2 * c[2] / c[0]
    k0_ = one(np.concatenate(blocks))
    bs = [one(np.concatenate([blocks[j] for j in rng.integers(0, 12, 12)])) for _ in range(a.n_boot)]
    m2 = np.mean(np.tan(r[:, 1]) ** 2)
    return dict(L=L, dt=dt, n=len(r), Lambda=1 + k0_, ci=1 + np.percentile(bs, [16, 84]), R=abs(1 + k0_) * m2 / 2, m2=m2)
# dt so that the mean advance is ~2 and ~5 px
if bubble:
    vpx = (Rk[k1] - Rk[k0]) / (k1 - k0)
else:
    sw_all = ((ta > k0) & (ta <= k1) & roi).sum(); Ltot = np.mean([sum(p[0].size for p in C[k]) * 0.5 for k in frames])
    vpx = sw_all / Ltot / (k1 - k0)
dts = sorted({max(1, int(round(2 / vpx))), max(1, int(round(5 / vpx)))})
print(f"mean advance {vpx:.2f} px/frame ({vpx * s:.3f} um/frame); dt = {dts} frames")
print("KPZ term: Lambda = 1 + kappa (0 elastic qEW, 1 isotropic, >> 1 qKPZ); R = |Lambda| <m^2>/2")
KP = []
for L in (() if a.skip_kpz else (8, 16, 32)):
    for dt in dts:
        f = kpz(L, dt)
        if f is None: continue
        KP.append(f)
        print(f"   L = {L:3d} px ({L * s:4.1f} um), dt = {dt:2d}: n = {f['n']:6d}, Lambda = {f['Lambda']:+.2f} [{f['ci'][0]:+.2f}, {f['ci'][1]:+.2f}], R = {f['R']:.3f}")

out = dict(E=E, w2=w2m, lams=lams, zeta_W=zW, zeta_W20=zW2, n_win=NT, **{f"P_tol{t:g}": P[t] for t in TOLS},
           **{f"l0g_tol{t:g}": l0g[t] for t in TOLS}, kpz=json.dumps([{k: (v.tolist() if hasattr(v, 'tolist') else v) for k, v in f.items()} for f in KP]))

# ---------------------------------------------------------------- bubble: macroscopic anisotropy and heights
if bubble:
    NTH = 720; g = -np.pi + (np.arange(NTH) + 0.5) * 2 * np.pi / NTH
    okth = np.ones(NTH, bool)
    if a.exclude_deg: okth = ~((np.degrees(g) >= a.exclude_deg[0]) & (np.degrees(g) <= a.exclude_deg[1]))
    Rth = []
    for k in frames:
        ys = np.concatenate([p[0] for p in C[k]]); xs = np.concatenate([p[1] for p in C[k]])
        th = np.arctan2(ys - y0, xs - x0); rr = np.hypot(ys - y0, xs - x0); o = np.argsort(th)
        Rth.append(np.interp(g, th[o], rr[o], period=2 * np.pi))
    Rth = np.array(Rth)
    V = (Rth[1:] - Rth[:-1]); V = V / V[:, okth].mean(1, keepdims=True)
    vbar = V.mean(0)
    def fourier_fit(f, nmax):
        A = [np.ones(okth.sum())]
        for n_ in range(1, nmax + 1): A += [np.cos(n_ * g[okth]), np.sin(n_ * g[okth])]
        c = np.linalg.lstsq(np.array(A).T, f[okth], rcond=None)[0]; vs = c[0] * np.ones(NTH); ll = c[0] * np.ones(NTH)
        amps = []
        for n_ in range(1, nmax + 1):
            cs, sn = c[2 * n_ - 1], c[2 * n_]; mode = cs * np.cos(n_ * g) + sn * np.sin(n_ * g)
            vs += mode; ll += (1 - n_ ** 2) * mode; amps.append(np.hypot(cs, sn) / c[0])
        return vs, ll, np.array(amps)
    blocks = np.array_split(np.arange(len(V)), 6)
    for nmax in (2, 4):
        vs, ll, amps = fourier_fit(vbar, nmax); Lam = (ll / vs)[okth]
        nulls = []
        for _ in range(200):
            VR = np.array([np.roll(v, rng.integers(NTH)) for v in V]); nulls.append(fourier_fit(VR.mean(0), nmax)[2])
        print(f"bubble anisotropy (n <= {nmax}, allowed angles): |dv_n|/v = " + " ".join(f"{x:.3f}" for x in amps) +
              f"  (noise 95%: " + " ".join(f"{x:.3f}" for x in np.percentile(nulls, 95, 0)) + f"); Lambda(phi) {Lam.min():.2f} .. {Lam.max():.2f}")
    # heights: r(theta) minus modes n <= 4 (fitted on the allowed arc), skewness per frame
    sk = []
    for r_ in Rth:
        vs, _, _ = fourier_fit(r_, 4); u = (r_ - vs)[okth]; sk.append(skew(u))
    sk = np.array(sk)
    print(f"heights: skewness of r(theta) - modes n <= 4 (allowed arc): by thirds " +
          ", ".join(f"{np.mean(sk[p]):+.2f}" for p in np.array_split(np.arange(len(sk)), 3)) + f" (frame-to-frame std {sk.std():.2f})")
    out.update(Rth=Rth, theta=g, okth=okth, vbar=vbar, skew=sk, Rk=np.array([Rk[k] for k in frames]))
np.savez(os.path.join(a.outdir, "contour_analysis.npz"), frames=np.array(frames), **out)

fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
ax[0].loglog(E, w2m, "o-", ms=3, label=f"ζ = {zW:.2f} (8–80 px)"); ax[0].set_xlabel("E [px]"); ax[0].set_ylabel("⟨w²⟩" + (" − circle" if bubble else "") + " [px²]"); ax[0].legend()
for t, st in zip(TOLS, (":", "-", "--")): ax[1].semilogx(E, P[t], st, marker="o", ms=3, label=f"tol {t:g} px")
ax[1].axhline(0.5, color="0.5", lw=0.8); ax[1].set_ylim(-0.02, 1.02); ax[1].legend(); ax[1].set_xlabel("E [px]"); ax[1].set_ylabel("P_ov")
for dt in dts:
    fs = [f for f in KP if f["dt"] == dt]
    ax[2].errorbar([f["L"] for f in fs], [f["Lambda"] for f in fs], [[f["Lambda"] - f["ci"][0] for f in fs], [f["ci"][1] - f["Lambda"] for f in fs]], marker="o", capsize=3, label=f"dt = {dt}")
ax[2].axhline(1, color="0.5", ls=":"); ax[2].axhline(0, color="0.5", ls=":"); ax[2].set_xscale("log"); ax[2].set_xlabel("L [px]"); ax[2].set_ylabel("Λ"); ax[2].legend()
fig.tight_layout(); fig.savefig(os.path.join(a.outdir, "fig_contour_analysis.png"), dpi=100)
