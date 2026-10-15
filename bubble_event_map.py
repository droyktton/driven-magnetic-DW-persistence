"""Expanding bubble from an event map (one avalanche index per pixel): roughness, overhangs, avalanches.

Usage: python bubble_event_map.py EVENTS.csv --outdir OUT [--snapshots 300] [--rmin 8] [--tol 3]

Input: a 2D map (CSV, no header) where each pixel holds the index of the event (avalanche) in which it
switched, in time order, or -1 if no switch was detected (format of the mokas / Baiesi event files).
There are no frame times: the event index is used as the clock. No pixel size: lengths in px.

1. Arrival map. The -1 pixels are split into the initial bubble (the component that contains the
   centroid of the switched pixels), the outside (components touching the image border) and holes
   (the rest; filled with the index of the nearest switched pixel, as the defects of step1).
   Snapshot k: domain = initial bubble + pixels with index <= k, holes filled; --snapshots values of k
   equally spaced in the number of switched pixels (i.e. in area) between 10 % and 100 %.
2. Roughness of the closed wall, two estimators:
   a. Polar: centre = centroid of the domain; r(theta) of the subpixel contour (domain smoothed by 1 px,
      as in step7), interpolated on a uniform grid of --sector px of arc (where the wall folds back in theta,
      the points are sorted by theta: P_ov is ~0 in these data, see 3). Counting pixels in angular sectors
      (area-conserving radius) is not usable here: the wall is rough by ~1 px only, and the pixel-counting
      noise of thin wedges makes S rise at high q; u(s) = r - <r>, s = <r> theta. The wall is closed, so u is periodic: S(q) needs no
      window. S(q) ~ q^-(1+2 zeta), averaged over snapshots in log bins of q; fit over
      2 pi / (L/4) < q < 2 pi / --rmin-sq, L = 2 pi <r> of the smallest snapshot used (the lowest modes are
      the shape of the bubble: ellipticity, off-centre growth).
   b. Local width on the contour (subpixel, domain smoothed by 1 px as in step7), arc windows of length
      lambda each rotated to its own axis (wall_tools.overhang_segments); the variance of a perfect arc of
      the same length on a circle of radius <r> is subtracted (curvature). lambda <= <r>.
3. Overhangs: P_ov(E) of step11 on the same windows (back step > --tol px along the own axis), with
   lambda <= <r>, below which the curvature of the circle cannot fold a window back.
4. Avalanches: P(S >= s); exponent tau of P(s) ~ s^-tau by maximum likelihood truncated to
   [5 px, --smax-tau] (P(S) has a cutoff); extent of
   each event along the wall (tangent to the circle at its barycentre) l_par and across it l_perp;
   S ~ l_par^(1 + zeta) gives zeta from the avalanche shapes (fit over S >= 10 px).
Outputs in OUT: bubble.npz, fig_bubble_map.png, fig_bubble_roughness.png, fig_bubble_avalanches.png,
summary on stdout.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import binary_fill_holes, distance_transform_edt, gaussian_filter, label

from wall_tools import contour_points, fit_power, logbin_xy, overhang_segments

ap = argparse.ArgumentParser()
ap.add_argument("events", help="CSV event map")
ap.add_argument("--outdir", required=True)
ap.add_argument("--snapshots", type=int, default=300)
ap.add_argument("--rmin", type=float, default=8, help="smallest scale of the local-width fit, px (optical resolution)")
ap.add_argument("--rmin-sq", type=float, default=20, help="smallest scale 2 pi / q of the polar S(q) fit, px")
ap.add_argument("--sector", type=float, default=1, help="arc length of the grid of r(theta), px")
ap.add_argument("--smax-tau", type=float, default=30, help="upper end of the power-law fit of P(S), px")
ap.add_argument("--tol", type=float, default=3.0, help="overhang tolerance, px (as step11)")
args = ap.parse_args()
os.makedirs(args.outdir, exist_ok=True)
out = lambda n: os.path.join(args.outdir, n)

# ------------------------------------------------------------------ 1. arrival map
ev = np.loadtxt(args.events, delimiter=",", dtype=np.int64)
H, W = ev.shape
sw = ev >= 0
cy, cx = np.argwhere(sw).mean(0)
lab, _ = label(~sw)
init_id = lab[int(round(cy)), int(round(cx))]
if init_id == 0:
    raise SystemExit("ERROR: the centroid of the switched pixels is switched, no initial bubble found")
border = np.setdiff1d(np.unique(np.concatenate([lab[0], lab[-1], lab[:, 0], lab[:, -1]])), [0])
init = lab == init_id
outside = np.isin(lab, border)
holes = ~sw & ~init & ~outside
ta = np.where(sw, ev, 0).astype(float)
ta[init] = -1
ta[outside] = np.inf
if holes.any():
    _, (iy, ix) = distance_transform_edt(~sw, return_indices=True)
    ta[holes] = ev[iy[holes], ix[holes]]
ids, sizes = np.unique(ev[sw], return_counts=True)
print(f"event map {H}x{W}: {ids.size} events (index {ids.min()}..{ids.max()}), {sw.sum()} switched px, "
      f"initial bubble {init.sum()} px (r_eq = {np.sqrt(init.sum() / np.pi):.0f} px), {holes.sum()} hole px filled")

# snapshots equally spaced in switched area
order = np.sort(ta[np.isfinite(ta) & (ta >= 0)])
ks = order[np.linspace(0.1 * (order.size - 1), order.size - 1, args.snapshots).astype(int)]

conts, polar, Rbar = [], [], []
yy, xx = np.mgrid[0:H, 0:W]
for k in ks:
    dom = binary_fill_holes(ta <= k)
    lab2, _ = label(dom)
    dom = lab2 == lab2[int(round(cy)), int(round(cx))]
    y0, x0 = yy[dom].mean(), xx[dom].mean()
    R = np.sqrt(dom.sum() / np.pi)
    c_ = contour_points(gaussian_filter(dom.astype(float), 1.0))
    conts.append(c_)
    th = np.arctan2(c_[0] - y0, c_[1] - x0); rr = np.hypot(c_[0] - y0, c_[1] - x0)
    o_ = np.argsort(th); th, rr = th[o_], rr[o_]
    nb = int(round(2 * np.pi * R / args.sector))
    grid = -np.pi + (np.arange(nb) + 0.5) * 2 * np.pi / nb
    r = np.interp(grid, th, rr, period=2 * np.pi)
    polar.append(r); Rbar.append(r.mean())
Rbar = np.array(Rbar)
print(f"snapshots: {len(ks)}, mean radius {Rbar[0]:.0f} -> {Rbar[-1]:.0f} px")

# ------------------------------------------------------------------ 2a. polar S(q)
qs, Ss = [], []
for r, R in zip(polar, Rbar):
    u = r - r.mean(); n = u.size; L = 2 * np.pi * R; ds = L / n
    S = np.abs(np.fft.rfft(u) * ds) ** 2 / L
    q = 2 * np.pi * np.arange(S.size) / L
    qs.append(q[1:]); Ss.append(S[1:])
qa, Sa = np.concatenate(qs), np.concatenate(Ss)
qb, Sb = logbin_xy(qa, Sa, 10)
Lmin = 2 * np.pi * Rbar.min()
qlo, qhi = 2 * np.pi / (Lmin / 4), 2 * np.pi / args.rmin_sq
ok = (qb >= qlo) & (qb <= qhi)
slope = np.polyfit(np.log(qb[ok]), np.log(Sb[ok]), 1)[0]
zS = (-slope - 1) / 2
# by thirds of the growth
zS_th = []
for part in np.array_split(np.arange(len(polar)), 3):
    qa_, Sa_ = np.concatenate([qs[i] for i in part]), np.concatenate([Ss[i] for i in part])
    qb_, Sb_ = logbin_xy(qa_, Sa_, 10)
    m_ = (qb_ >= 2 * np.pi / (2 * np.pi * Rbar[part].min() / 4)) & (qb_ <= qhi)
    zS_th.append((-np.polyfit(np.log(qb_[m_]), np.log(Sb_[m_]), 1)[0] - 1) / 2)
print(f"polar S(q), r(theta) of the contour: zeta = {zS:.2f} over {args.rmin_sq:g} < 2pi/q < {Lmin / 4:.0f} px "
      f"(thirds of the growth: " + ", ".join(f"{z:.2f}" for z in zS_th) + ")")

# ------------------------------------------------------------------ 2b, 3. local width and overhangs on the contour
TOLS = sorted({2.0, args.tol, 5.0})
lam_max = Rbar.min()
lams = np.unique(np.round(np.logspace(np.log10(4), np.log10(lam_max), 24)))
NT, NO, EX, W2 = [], {t: [] for t in TOLS}, [[] for _ in lams], [[] for _ in lams]
for (y, x), R in zip(conts, Rbar):
    nt, no, ex, w2 = overhang_segments(y, x, lams, 0.5, TOLS)
    NT.append(nt)
    for t in TOLS:
        NO[t].append(no[t])
    for j, lam in enumerate(lams):
        a = np.linspace(-lam / (2 * R), lam / (2 * R), 200)          # perfect arc of the same length
        circ = np.var(R * np.cos(a))
        EX[j].append(ex[j]); W2[j].append(w2[j] - circ)
NT = np.array(NT)
E = np.array([np.median(np.concatenate(e)) for e in EX])
w2m = np.array([np.mean(np.concatenate(w)) for w in W2])
P = {t: np.array(NO[t]).sum(0) / np.maximum(NT.sum(0), 1) for t in TOLS}
okw = (E >= args.rmin) & (E <= lam_max / 2) & (w2m > 0)
zW = fit_power(E[okw], w2m[okw], args.rmin, lam_max / 2)[0] / 2 if okw.sum() > 2 else np.nan
print(f"local width on the contour (arc windows, own axis, circle subtracted): zeta = {zW:.2f} "
      f"over {args.rmin:g} < E < {lam_max / 2:.0f} px")
for j in np.unique(np.searchsorted(E, [5, 10, 20, 50, 100, 200])).clip(0, len(E) - 1):
    print(f"   E = {E[j]:6.1f} px (arc {lams[j]:5.0f}): w = {np.sqrt(max(w2m[j], 0)):5.2f} px, "
          f"P_ov(tol {args.tol:g}) = {P[args.tol][j]:.3f}  ({int(NT[:, j].sum())} windows)")


def half(xv, pv):
    a = np.nonzero(pv >= 0.5)[0]
    if a.size == 0:
        return np.nan
    j = a[0]
    if j == 0:
        return xv[0]
    f = (0.5 - pv[j - 1]) / (pv[j] - pv[j - 1])
    return float(np.exp(np.log(xv[j - 1]) + f * (np.log(xv[j]) - np.log(xv[j - 1]))))


l0g = {t: half(E, P[t]) for t in TOLS}
print("   l0_geom (P_ov = 1/2): " + ", ".join(f"tol {t:g}: {'> %.0f' % E.max() if np.isnan(v) else '%.0f' % v} px"
                                           for t, v in l0g.items()))

# ------------------------------------------------------------------ 4. avalanches
flat = ev.ravel(); ii = np.nonzero(flat >= 0)[0]
o = np.argsort(flat[ii], kind="stable"); ii = ii[o]
cut = np.nonzero(np.diff(flat[ii]))[0] + 1
groups = np.split(ii, cut)
S = np.array([g.size for g in groups])
lpar, lperp = np.zeros(S.size), np.zeros(S.size)
for n, g in enumerate(groups):
    gy, gx = g // W, g % W
    by, bx = gy.mean(), gx.mean()
    rad = np.array([by - cy, bx - cx]); rad /= max(np.hypot(*rad), 1e-9)
    tan = np.array([-rad[1], rad[0]])
    pt = gy * tan[0] + gx * tan[1]; pr = gy * rad[0] + gx * rad[1]
    lpar[n] = pt.max() - pt.min() + 1; lperp[n] = pr.max() - pr.min() + 1
smin, smax = 5, args.smax_tau
sv = np.arange(smin, int(smax) + 1)
cnt = np.array([(S == v).sum() for v in sv])
taus = np.linspace(0.5, 3.5, 3001)
ll = np.array([-(t * np.log(sv) * cnt).sum() - cnt.sum() * np.log((sv ** -t).sum()) for t in taus])
tau = taus[np.argmax(ll)]
curv = (ll[np.argmax(ll) + 1] - 2 * ll.max() + ll[np.argmax(ll) - 1]) / (taus[1] - taus[0]) ** 2
dtau = 1 / np.sqrt(-curv)
big = S >= 10
lb, Sbn = logbin_xy(lpar[big], S[big].astype(float), 8)
zA = np.polyfit(np.log(lpar[big]), np.log(S[big]), 1)[0] - 1
lpb, _ = logbin_xy(lpar[big], lperp[big], 8)
print(f"avalanches: {S.size} events, S median {np.median(S):.0f} px, max {S.max()}; "
      f"tau (ML truncated to {smin}-{smax:g} px) = {tau:.2f} ± {dtau:.2f}; P(S >= 50) = {(S >= 50).mean():.3f}")
print(f"   shape: S ~ l_par^(1+zeta) over S >= 10 px ({big.sum()} events): zeta = {zA:.2f}; "
      f"median l_par / l_perp = {np.median(lpar[big] / lperp[big]):.1f}")

np.savez(out("bubble.npz"), ks=ks, Rbar=Rbar, qb=qb, Sb=Sb, zeta_S=zS, zeta_S_thirds=zS_th, lams=lams, E=E,
         w2=w2m, zeta_W=zW, n_win=NT.sum(0), **{f"P_tol{t:g}": P[t] for t in TOLS},
         **{f"l0_geom_tol{t:g}": l0g[t] for t in TOLS}, S=S, l_par=lpar, l_perp=lperp, tau=tau, dtau=dtau,
         zeta_aval=zA, center=(cy, cx))

# ------------------------------------------------------------------ figures
fig, ax = plt.subplots(1, 2, figsize=(14, 7))
show = np.where(np.isfinite(ta) & (ta >= 0), ta, np.nan)
ax[0].imshow(show, cmap="viridis"); ax[0].set_title("arrival map = event index (initial bubble and outside blank)")
for j in np.linspace(0, len(conts) - 1, 8).astype(int):
    ax[1].plot(conts[j][1], conts[j][0], lw=0.8, color=plt.cm.viridis(j / (len(conts) - 1)))
ax[1].set_aspect("equal"); ax[1].invert_yaxis(); ax[1].set_title("wall at 8 snapshots")
fig.tight_layout(); fig.savefig(out("fig_bubble_map.png"), dpi=90); plt.close(fig)

fig, ax = plt.subplots(1, 3, figsize=(18, 5))
ax[0].loglog(qb, Sb, "o", ms=3)
ax[0].loglog(qb[ok], np.exp(np.polyval(np.polyfit(np.log(qb[ok]), np.log(Sb[ok]), 1), np.log(qb[ok]))), "r-",
             label=f"ζ = {zS:.2f}")
for z_, c_ in ((1.25, "C2"), (0.63, "C1")):
    ax[0].loglog(qb[ok], Sb[ok][0] * (qb[ok] / qb[ok][0]) ** (-(1 + 2 * z_)), ":", color=c_, label=f"ζ = {z_}")
ax[0].set_xlabel("q [1/px]"); ax[0].set_ylabel("S(q)"); ax[0].legend(); ax[0].set_title("polar S(q), r(θ) of the subpixel contour")
ax[1].loglog(E, w2m, "o-", ms=3, label=f"ζ = {zW:.2f}"); ax[1].set_xlabel("E [px]"); ax[1].set_ylabel("⟨w²⟩ − circle [px²]")
ax[1].axvspan(args.rmin, lam_max / 2, color="0.92", zorder=0); ax[1].legend(); ax[1].set_title("local width, arc windows")
for t, st in zip(TOLS, (":", "-", "--")):
    ax[2].semilogx(E, P[t], st, marker="o", ms=3, label=f"tol {t:g} px")
ax[2].axhline(0.5, color="0.5", lw=0.8); ax[2].set_ylim(-0.02, 1.02); ax[2].legend()
ax[2].set_xlabel("E [px]"); ax[2].set_ylabel("P_ov"); ax[2].set_title("overhangs (step11 definition)")
fig.tight_layout(); fig.savefig(out("fig_bubble_roughness.png"), dpi=100); plt.close(fig)

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
s_ = np.sort(S)[::-1]
ax[0].loglog(s_, np.arange(1, s_.size + 1) / s_.size, ".", ms=2)
sx = np.logspace(np.log10(smin), np.log10(smax), 20)
ax[0].loglog(sx, (S >= smin).mean() * (sx / smin) ** (1 - tau), "r-", label=f"τ = {tau:.2f} ({smin}-{smax:g} px)")
ax[0].set_xlabel("S [px]"); ax[0].set_ylabel("P(S ≥ s)"); ax[0].legend()
ax[1].loglog(lpar[big], S[big], ".", ms=1, alpha=0.3, color="0.6")
ax[1].loglog(lb, Sbn, "o", ms=4)
ax[1].loglog(lb, np.exp(np.polyval(np.polyfit(np.log(lpar[big]), np.log(S[big]), 1), np.log(lb))), "r-",
             label=f"S ~ ℓ^(1+ζ), ζ = {zA:.2f}")
ax[1].set_xlabel("extent along the wall ℓ_par [px]"); ax[1].set_ylabel("S [px]"); ax[1].legend()
fig.tight_layout(); fig.savefig(out("fig_bubble_avalanches.png"), dpi=100); plt.close(fig)
print(f"outputs in {args.outdir}")
