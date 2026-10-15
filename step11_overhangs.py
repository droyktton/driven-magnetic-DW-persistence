"""Step 11: onset of overhangs measured directly on the wall contour, l0_geom.

Usage: python step11_overhangs.py DIR [DIR ...] [--tol 3] [--n-boot 300]
       python step11_overhangs.py --selftest
Needs DIR/contours.npz (step7, TIFF studies) and, for the comparison, DIR/crossover.npz (step8).
Outputs in DIR: overhangs.npz, fig_overhangs.png; summary on stdout (one table for all DIRs).

Why. step8 finds the crossover l0 as a change of slope of w2(l) or S(q), from zeta ~ 1-1.25 to
zeta_eff ~ 0.5. But zeta = 1/2 is expected on both sides of two different crossovers: above l0
(effective, the wall is multivalued; Kolton, Ferrero & Rosso) and above the avalanche scale xi_av of
creep (genuine, thermal, the wall is still single-valued). The break measured by the exponent is
min(xi_av, l0). Whether the wall is single-valued at a given scale does not depend on zeta, so it
separates the two: here l0 is measured as the scale where segments stop being single-valued.

Estimator. The contour of each frame (ordered along the arc, 0.5 px steps, domain smoothed by 1 px
as saved by step7) is cut into windows of arc length lambda, sliding by lambda/2. Each window is
rotated to its own principal axis (orthogonal least-squares line, as the local width of step7) and
its points are projected on that axis, p(s), oriented so that p increases from the first to the
last point. The segment is single-valued in its own frame if p(s) never goes back; it has an
overhang if max_s [max_{s'<=s} p(s') - p(s)] > --tol px (pixel staircase and the 1 px smoothing
make back steps of < 1-2 px; the optical resolution is ~3-4 px, so the default is 3 px).
  P_ov(lambda)  fraction of segments with an overhang (increases with lambda, -> 1 above l0);
  E(lambda)     median extent of the segments along their own axis (= the l of step7/8, < lambda
                because of the wiggles);
  l0_geom       E at the lambda where P_ov = 1/2 (log interpolation); 16-84 % block bootstrap over
                time (12 blocks of frames); repeated with tol = 2 and 5 px.
Control: overhangs on the base plane (fixed direction theta0, no local rotation): number of back
steps > tol of the projection on the base direction, per unit projected length; 1/density is the
mean distance between overhangs along the wall (it includes the macroscopic tilt of the wall).
Comparison with step8: P_ov at the step8 crossovers (l0 of the local width, a, and of the global
S(q), b). P_ov(l_x) << 1/2 means the wall is still single-valued at the break, i.e. the break is
not the overhang crossover (candidate: xi_av); P_ov(l_x) ~ 1/2 means it is.
Self-test (--selftest): synthetic walls tilted by -12 deg, self-affine with zeta = 1 (rms 8 px):
without overhangs (false-positive rate) and with mushroom-shaped tongues every D px (D = 100, 200,
400): l0_geom must scale with D.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from wall_tools import contour_points

ap = argparse.ArgumentParser()
ap.add_argument("dirs", nargs="*", help="study folders")
ap.add_argument("--tol", type=float, default=3.0, help="smallest back step along the own axis counted as overhang, px")
ap.add_argument("--n-boot", type=int, default=300, help="bootstrap resamples")
ap.add_argument("--selftest", action="store_true", help="synthetic walls with and without overhangs")
args = ap.parse_args()
TOLS = sorted({2.0, args.tol, 5.0})


def segments(y, x, lams, step, tols):
    """Per arc length lambda: number of windows, number with overhang (per tol), extents along own axis."""
    n_tot, n_ov, ext = [], {t: [] for t in tols}, []
    for lam in lams:
        n = int(round(lam / step)) + 1
        if y.size < n:
            n_tot.append(0); ext.append(np.empty(0)); [n_ov[t].append(0) for t in tols]; continue
        hop = max(n // 2, 1)
        Y = sliding_window_view(y, n)[::hop]; X = sliding_window_view(x, n)[::hop]
        my, mx = Y.mean(1, keepdims=True), X.mean(1, keepdims=True)
        dy, dx = Y - my, X - mx
        cyy, cxx, cxy = (dy * dy).mean(1), (dx * dx).mean(1), (dy * dx).mean(1)
        al = 0.5 * np.arctan2(2 * cxy, cyy - cxx)                  # major axis, from y towards x
        p = dy * np.cos(al)[:, None] + dx * np.sin(al)[:, None]
        p *= np.where(p[:, -1] >= p[:, 0], 1.0, -1.0)[:, None]
        back = (np.maximum.accumulate(p, axis=1) - p).max(1)
        n_tot.append(p.shape[0]); ext.append(p.max(1) - p.min(1))
        for t in tols:
            n_ov[t].append(int((back > t).sum()))
    return np.array(n_tot), {t: np.array(v) for t, v in n_ov.items()}, ext


def base_events(y, x, theta0_deg, tol):
    """Back steps > tol of the projection on the base direction; returns (events, projected length)."""
    t0 = np.radians(theta0_deg)
    s = y * np.cos(t0) + x * np.sin(t0)
    if s[-1] < s[0]:
        s = -s
    deficit = np.maximum.accumulate(s) - s
    inside = deficit > tol
    return int(np.count_nonzero(np.diff(inside.astype(np.int8)) == 1) + inside[0]), float(s.max() - s.min())


def half_point(xv, pv, level=0.5):
    """First x where pv crosses level upwards (log interpolation); nan if never."""
    above = np.nonzero(pv >= level)[0]
    if above.size == 0 or not np.isfinite(pv).any():
        return np.nan
    j = above[0]
    if j == 0:
        return xv[0]
    f = (level - pv[j - 1]) / (pv[j] - pv[j - 1])
    return float(np.exp(np.log(xv[j - 1]) + f * (np.log(xv[j]) - np.log(xv[j - 1]))))


def analyse(conts, theta0, step, lams):
    """conts: list of (y, x) per frame. Returns per-frame counts and extents."""
    NT, NO, EX, EV, PL = [], {t: [] for t in TOLS}, [[] for _ in lams], [], []
    for y, x in conts:
        nt, no, ex = segments(y, x, lams, step, TOLS)
        NT.append(nt)
        for t in TOLS:
            NO[t].append(no[t])
        for j, e in enumerate(ex):
            EX[j].append(e)
        ev, pl = base_events(y, x, theta0, args.tol)
        EV.append(ev); PL.append(pl)
    NT = np.array(NT); NO = {t: np.array(v) for t, v in NO.items()}
    E = np.array([np.median(np.concatenate(e)) if sum(a.size for a in e) else np.nan for e in EX])
    return NT, NO, E, np.array(EV), np.array(PL)


def l0_geom(NT, NOt, E, pick=None):
    s = slice(None) if pick is None else pick
    P = NOt[s].sum(0) / np.maximum(NT[s].sum(0), 1)
    P[NT[s].sum(0) < 20] = np.nan
    ok = np.isfinite(P) & np.isfinite(E)
    return half_point(E[ok], P[ok]), P


# ------------------------------------------------------------------ self test
if args.selftest:
    from scipy.ndimage import gaussian_filter
    rng = np.random.default_rng(1)
    H, W, theta = 1400, 900, -12.0
    Y, X = np.mgrid[0:H, 0:W]
    lams = np.unique(np.round(np.logspace(np.log10(8), np.log10(H / 2), 24)))
    k = np.fft.rfftfreq(H)

    def wall(D):
        amp = np.zeros_like(k); amp[1:] = k[1:] ** (-(1 + 2 * 1.0) / 2)
        u = np.fft.irfft(amp * np.exp(2j * np.pi * rng.random(k.size)), n=H); u *= 8 / u.std()
        f = 400 + (Y[:, 0] - H / 2) * np.tan(np.radians(theta)) + u
        dom = X < f[:, None]
        if D:
            for yc in np.arange(rng.uniform(0, D), H, D).astype(int):       # mushrooms: stem + cap
                x0 = int(f[yc]); depth = int(rng.uniform(25, 45)); capw = int(rng.uniform(12, 25))
                dom[max(yc - 3, 0):yc + 4, x0 - 2:x0 + depth] = True
                dom[max(yc - capw, 0):yc + capw + 1, x0 + depth - 8:x0 + depth] = True
        return contour_points(gaussian_filter(dom.astype(float), 1.0))

    print(f"selftest: tilted wall ({theta} deg), zeta = 1, rms 8 px, tol = {args.tol:g} px; 8 walls each")
    for D in (0, 100, 200, 400):
        conts = [wall(D) for _ in range(8)]
        NT, NO, E, EV, PL = analyse(conts, theta, 0.5, lams)
        l0, P = l0_geom(NT, NO[args.tol], E)
        dens = EV.sum() / PL.sum()
        lab = "no overhangs" if D == 0 else f"mushrooms every {D} px"
        print(f"  {lab:24s}: P_ov max {np.nanmax(P):.3f} (at E = {E[np.nanargmax(P)]:.0f} px), l0_geom = {l0:6.1f} px, "
              f"l0_geom / D = {l0 / D if D else np.nan:.2f}; base-plane spacing = {1 / dens if dens else np.inf:.0f} px")
    raise SystemExit


# ------------------------------------------------------------------ studies
if not args.dirs:
    ap.error("give study folders or --selftest")
rows = []
for d in args.dirs:
    c = np.load(os.path.join(d, "contours.npz"))
    off, step, theta0 = c["offsets"], float(c["step"]), float(c["theta0"])
    conts = [(c["y"][a:b], c["x"][a:b]) for a, b in zip(off[:-1], off[1:])]
    import json
    with open(os.path.join(d, "meta.json")) as f:
        k_um = json.load(f).get("um_per_px") or 0.0
    kk, unit = (k_um, "µm") if k_um else (1.0, "px")
    arc_min = min(y.size for y, _ in conts) * step
    lams = np.unique(np.round(np.logspace(np.log10(8), np.log10(0.95 * arc_min), 32)))
    NT, NO, E, EV, PL = analyse(conts, theta0, step, lams)
    blocks = np.array_split(np.arange(NT.shape[0]), 12)
    rng = np.random.default_rng(5)
    res = {}
    for t in TOLS:
        l0, P = l0_geom(NT, NO[t], E)
        bs = [l0_geom(NT, NO[t], E, np.concatenate([blocks[j] for j in rng.integers(0, 12, 12)]))[0]
              for _ in range(args.n_boot)]
        bs = np.array(bs)
        frac_nan = np.isnan(bs).mean()
        ci = np.nanpercentile(bs, [16, 84]) if np.isfinite(bs).any() else np.array([np.nan, np.nan])
        res[t] = dict(l0=l0, P=P, ci=ci, frac_nan=frac_nan)
    dens = EV.sum() / PL.sum()
    cx = {}
    cpath = os.path.join(d, "crossover.npz")
    if os.path.exists(cpath):
        z = np.load(cpath)
        for e in ("a", "b"):
            lx = float(z[f"l0_{e}"]) / kk                         # back to px
            P = res[args.tol]["P"]; ok = np.isfinite(P)
            cx[e] = (lx * kk, float(np.interp(np.log(lx), np.log(E[ok]), P[ok], left=np.nan, right=np.nan)))
    r = res[args.tol]
    print(f"\n=== {d}: {len(conts)} frames, theta0 = {theta0:.1f} deg, segments up to {E[np.isfinite(E)].max() * kk:.0f} {unit}")
    for t in TOLS:
        rt = res[t]
        bnd = "" if np.isfinite(rt["l0"]) else f"  (P_ov < 1/2 up to the largest scale: l0_geom > {np.nanmax(E) * kk:.0f})"
        nr = "" if rt["frac_nan"] == 0 else f", {rt['frac_nan']:.0%} of resamples never reach 1/2"
        print(f"   tol {t:g} px: l0_geom = {rt['l0'] * kk:6.1f} {unit} (bootstrap {rt['ci'][0] * kk:.1f}-{rt['ci'][1] * kk:.1f}"
              f"{nr}){bnd}")
    for j in np.unique(np.searchsorted(E * kk, [2, 5, 10, 20, 40])).clip(0, len(E) - 1):
        print(f"   E = {E[j] * kk:6.1f} {unit} (arc {lams[j] * kk:6.1f}): P_ov = {r['P'][j]:.2f}  ({int(NT[:, j].sum())} segments)")
    print(f"   base plane: {EV.sum()} overhangs > {args.tol:g} px in {PL.sum() * kk:.0f} {unit} of wall, "
          f"mean spacing {1 / dens * kk if dens else np.inf:.1f} {unit}")
    for e, (lx, px) in cx.items():
        print(f"   step8 l0_{e} = {lx:.1f} {unit}: P_ov there = {px:.2f}")
    np.savez(os.path.join(d, "overhangs.npz"), lams=lams * kk, E=E * kk, n_seg=NT.sum(0), tol=args.tol, tols=TOLS,
             **{f"P_tol{t:g}": res[t]["P"] for t in TOLS}, **{f"l0_geom_tol{t:g}": res[t]["l0"] * kk for t in TOLS},
             **{f"l0_geom_ci_tol{t:g}": res[t]["ci"] * kk for t in TOLS},
             base_events=EV.sum(), base_length=PL.sum() * kk, base_spacing=(1 / dens * kk if dens else np.inf),
             um_per_px=k_um, **{f"P_at_l0_{e}": v[1] for e, v in cx.items()})
    rows.append((d, r["l0"] * kk, r["ci"] * kk, res[2.0]["l0"] * kk, res[5.0]["l0"] * kk, 1 / dens * kk if dens else np.inf,
                 cx.get("a", (np.nan, np.nan)), cx.get("b", (np.nan, np.nan)), unit))

    fig, ax = plt.subplots(1, 2, figsize=(13, 4.6))
    for t, st in zip(TOLS, (":", "-", "--")):
        ax[0].semilogx(E * kk, res[t]["P"], st, marker="o", ms=3, label=f"tol {t:g} px: l0_geom = {res[t]['l0'] * kk:.1f} {unit}")
    ax[0].axhline(0.5, color="0.5", lw=0.8)
    for e, (lx, px), col in zip(cx, cx.values(), ("C3", "C4")):
        ax[0].axvline(lx, color=col, ls="-.", lw=1.2, label=f"step8 l0_{e} = {lx:.1f} {unit} (P_ov = {px:.2f})")
    ax[0].set_xlabel(f"segment extent along its own axis E [{unit}]"); ax[0].set_ylabel("P_ov: fraction with an overhang")
    ax[0].set_ylim(-0.02, 1.02); ax[0].legend(fontsize=8); ax[0].set_title(os.path.basename(os.path.normpath(d)), fontsize=10)
    ax[1].loglog(lams * kk, E * kk, "o-", ms=3)
    ax[1].loglog(lams * kk, lams * kk, "k:", lw=0.8)
    ax[1].set_xlabel(f"arc length λ [{unit}]"); ax[1].set_ylabel(f"median extent E [{unit}]")
    ax[1].set_title("Extent along the own axis vs arc length (dotted: E = λ)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(d, "fig_overhangs.png"), dpi=110); plt.close(fig)

print("\nsummary (tol = %g px; l0_geom: P_ov = 1/2)" % args.tol)
print(f"{'study':34s} {'l0_geom':>8s} {'16-84 %':>13s} {'tol 2':>7s} {'tol 5':>7s} {'spacing':>8s} {'l0_a':>6s} {'P_ov(l0_a)':>10s} {'l0_b':>6s} {'P_ov(l0_b)':>10s}")
for d, l0, ci, l2, l5, sp, (la, pa), (lb, pb), unit in rows:
    print(f"{os.path.basename(os.path.normpath(d)):34s} {l0:8.1f} {ci[0]:6.1f}-{ci[1]:<6.1f} {l2:7.1f} {l5:7.1f} {sp:8.1f} "
          f"{la:6.1f} {pa:10.2f} {lb:6.1f} {pb:10.2f}   [{unit}]")
