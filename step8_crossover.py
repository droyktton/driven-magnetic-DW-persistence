"""Step 8: crossover length l0 between the small-scale (qEW-like) and large-scale roughness regimes.

Usage: python step8_crossover.py DIR [DIR2 ...] [--m 2] [--n-boot 200]
       python step8_crossover.py --selftest            (synthetic walls with a known l0)
Input: DIR/local_width.npz and DIR/roughness.npz from step7 (with the per-time-block arrays).
Outputs: DIR/crossover.npz, DIR/fig_crossover.png and a table on stdout; with several DIRs also
fig_crossover_compare.png in the first one.

Picture (Kolton, Ferrero & Rosso, arXiv:2306.13415): below l0 the wall is qEW-like (zeta1 ~ 1.25),
above it zeta_eff ~ 0.5. Three estimates of l0:
  a. local width with local rotation (step7, PCA per segment), smooth broken power law
       w2(l) = A l^(2 z1) [1 + (l/l0)^m]^((2 z2 - 2 z1)/m)
     with z1 = 1.25, z2 = 0.5 fixed ("fixed") or free ("free"); on the contour, the contour without
     folded segments and u_area (or h); lower cutoff 2, 3 or 5 um (16, 25, 40 px without a scale);
  b. global S(q) of the wall on its base plane, the same crossover in q:
       S(q) = A q^-(1 + 2 z2) [1 + (q l0 / 2 pi)^m]^(-(2 z1 - 2 z2)/m)
     fitted between 2 pi / L and 2 pi / lmin;
  c. model free: the local slope zeta_eff(l) of w2 (over a factor 3 in l) crosses the midpoint
     between its values at small and large l.
Statistical errors: block bootstrap over time (12 blocks of frames, refitting each resample).
Systematic error: spread of the central values over the variants. When l0 is not bracketed by the
data (bootstrap fits beyond the largest scale), a lower bound is reported.
Disorder proxies from step7: l_fold (folded fraction down to half its small-scale value) and
l_tilt (sd of the local tilt down to half its plateau).
"""
import argparse
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import least_squares

ap = argparse.ArgumentParser()
ap.add_argument("dirs", nargs="*", help="study folders (step7 output)")
ap.add_argument("--m", type=float, default=2, help="sharpness of the crossover in the main variant")
ap.add_argument("--n-boot", type=int, default=200, help="bootstrap resamples")
ap.add_argument("--selftest", action="store_true", help="synthetic walls with a known l0")
args = ap.parse_args()
Z1, Z2 = 1.25, 0.5


# ------------------------------------------------------------------ models and fits
def model_w2(p, l, m):
    lnA, z1, z2, ll0 = p
    return lnA + 2 * z1 * np.log(l) + (2 * z2 - 2 * z1) / m * np.log1p((l / np.exp(ll0)) ** m)


def model_S(p, q, m):
    lnA, z1, z2, ll0 = p
    return lnA - (1 + 2 * z2) * np.log(q) - (2 * z1 - 2 * z2) / m * np.log1p((q * np.exp(ll0) / (2 * np.pi)) ** m)


def fit_cross(x, y, sig, kind, m, free, lo_l0, hi_l0):
    """Fit ln y with the broken power law; returns (l0, z1, z2) or nans. kind: 'w2' or 'S'."""
    ok = np.isfinite(y) & (y > 0) & np.isfinite(sig)
    x, ly, sig = x[ok], np.log(y[ok]), np.maximum(sig[ok], 0.02)
    if x.size < 5:
        return np.nan, np.nan, np.nan
    f = model_w2 if kind == "w2" else model_S
    best = None
    for g in np.log(np.geomspace(lo_l0 * 1.2, hi_l0 / 1.2, 5)):          # several starts for l0
        if free:
            p0, lb, ub = [0.0, Z1, Z2, g], [-np.inf, 0.0, 0.0, np.log(lo_l0)], [np.inf, 3.0, 2.0, np.log(hi_l0)]
            fun = lambda p: (f(p, x, m) - ly) / sig
        else:
            p0, lb, ub = [0.0, g], [-np.inf, np.log(lo_l0)], [np.inf, np.log(hi_l0)]
            fun = lambda p: (f([p[0], Z1, Z2, p[1]], x, m) - ly) / sig
        p0[0] = np.median(ly - f([0.0, Z1, Z2, g], x, m))
        r = least_squares(fun, p0, bounds=(lb, ub))
        if best is None or r.cost < best.cost:
            best = r
    p = best.x
    return (np.exp(p[3]), p[1], p[2]) if free else (np.exp(p[1]), Z1, Z2)


def local_slope(l, w2, fac=3.0):
    out = np.full(l.size, np.nan)
    for i, li in enumerate(l):
        s = (l >= li / np.sqrt(fac)) & (l <= li * np.sqrt(fac)) & np.isfinite(w2) & (w2 > 0)
        if s.sum() >= 3:
            out[i] = np.polyfit(np.log(l[s]), np.log(w2[s]), 1)[0] / 2
    return out


def l0_model_free(l, w2, lmin):
    """l where zeta_eff(l) first falls below the midpoint of its small- and large-scale values."""
    ze = local_slope(l, w2)
    sel = (l >= lmin) & np.isfinite(ze)
    if sel.sum() < 6:
        return np.nan, ze
    li, zi = l[sel], ze[sel]
    mid = (zi[:3].mean() + zi[-3:].mean()) / 2
    below = np.nonzero(zi < mid)[0]
    below = below[below > 0]
    if not below.size:
        return np.nan, ze
    k = below[0]
    t = (mid - zi[k - 1]) / (zi[k] - zi[k - 1])
    return float(np.exp(np.log(li[k - 1]) + t * np.log(li[k] / li[k - 1]))), ze


def logbin_edges(q):
    return np.logspace(np.log10(q.min()), np.log10(q.max() * 1.0001), int(np.log10(q.max() / q.min()) * 8) + 2)


def binned(q, S, e):
    i = np.digitize(q, e)
    u = np.unique(i)
    return np.array([np.exp(np.log(q[i == j]).mean()) for j in u]), np.array([S[..., i == j].mean(-1) for j in u]).T


# ------------------------------------------------------------------ analysis of one data set
def analyse(D, n_boot, m_main, rng):
    """D: dict with ells, blk_sum/blk_cnt per curve (12 x nl), q, Sblk (12 x nq), L, kk (um/px), unit."""
    kk, L = D["kk"], D["L"]
    ells = D["ells"]
    cuts = [2, 3, 5] if D["unit"] == "µm" else [16, 25, 40]
    cuts_px = [c / kk for c in cuts] if D["unit"] == "µm" else cuts
    lo_l0, hi_l0 = 2.0, 4 * L
    nb = D["Sblk"].shape[0]
    picks = [np.arange(nb)] + [rng.integers(0, nb, nb) for _ in range(n_boot)]

    def w2_of(curve, pk):
        s_, c_ = D["blk_sum"][curve][pk].sum(0), D["blk_cnt"][curve][pk].sum(0)
        return np.where(c_ > 0, s_ / np.maximum(c_, 1), np.nan)

    eS = logbin_edges(D["q"])
    qb, Sb_blk = binned(D["q"], D["Sblk"], eS)                       # nb x nbins
    rows = []
    for curve in D["blk_sum"]:
        allw = np.array([w2_of(curve, pk) for pk in picks])         # (1+n_boot) x nl
        sig = np.nanstd(np.log(allw[1:]), 0)
        for free in (False, True):
            for m in (1, 2, 4):
                for c_px, c_lab in zip(cuts_px, cuts):
                    sel = (ells >= c_px) & (ells <= L)
                    if m != m_main and c_lab != cuts[1]:
                        continue                                    # m varied only at the middle cutoff
                    res = [fit_cross(ells[sel], w[sel], sig[sel], "w2", m, free, lo_l0, hi_l0) for w in allw]
                    rows.append(dict(est="a", curve=curve, free=free, m=m, cut=c_lab, fits=np.array(res)))
    allS = np.array([Sb_blk[pk].mean(0) for pk in picks])
    sigS = np.nanstd(np.log(allS[1:]), 0)
    for free in (False, True):
        for m in (1, 2, 4):
            for c_px, c_lab in zip(cuts_px, cuts):
                if m != m_main and c_lab != cuts[1]:
                    continue
                sel = (qb <= 2 * np.pi / c_px) & (qb >= 2 * np.pi / L)
                res = [fit_cross(qb[sel], s[sel], sigS[sel], "S", m, free, lo_l0, hi_l0) for s in allS]
                rows.append(dict(est="b", curve="global S", free=free, m=m, cut=c_lab, fits=np.array(res)))
    mf_curve = "unfolded" if "unfolded" in D["blk_sum"] else list(D["blk_sum"])[0]
    mf = []
    for pk in picks:
        l0c, _ = l0_model_free(ells, w2_of(mf_curve, pk), cuts_px[1])
        mf.append(l0c)
    rows.append(dict(est="c", curve=mf_curve, free=None, m=None, cut=cuts[1],
                     fits=np.column_stack([mf, np.full(len(mf), np.nan), np.full(len(mf), np.nan)])))
    lmax = ells[ells <= L].max()
    for r in rows:
        f0, bt = r["fits"][0], r["fits"][1:, 0]
        bt = bt[np.isfinite(bt)]
        r["l0"], r["z1"], r["z2"] = f0
        r["lo"], r["hi"] = (np.percentile(bt, [16, 84]) if bt.size else (np.nan, np.nan))
        r["frac_beyond"] = np.mean(bt > lmax) if bt.size else np.nan
        r["bound"] = bool((np.isfinite(f0[0]) and f0[0] > lmax) or r["frac_beyond"] > 0.16)
    return rows, dict(qb=qb, S=allS[0], ells=ells, w2={c: w2_of(c, picks[0]) for c in D["blk_sum"]}, lmax=lmax,
                      cuts=cuts, cuts_px=cuts_px)


def summarize(rows):
    """Main value per estimator: fixed exponents, m = --m, middle cutoff, preferred curve; plus the
    spread of the central values over all variants (systematic)."""
    pref = ("unfolded", "contour", "pca", "global S")
    out = {}
    for est in ("a", "b", "c"):
        rr = [r for r in rows if r["est"] == est and np.isfinite(r["l0"])]
        if not rr:
            continue
        cuts = sorted({r["cut"] for r in rr})
        cand = [r for r in rr if r["free"] in (False, None) and r["m"] in (args.m, None) and r["cut"] == cuts[len(cuts) // 2]]
        r0 = next((r for c in pref for r in cand if r["curve"] == c), cand[0] if cand else rr[0])
        cen = np.array([r["l0"] for r in rr])
        out[est] = dict(l0=r0["l0"], lo=r0["lo"], hi=r0["hi"], bound=r0["bound"], curve=r0["curve"],
                        sys_lo=cen.min(), sys_hi=cen.max(), n_var=len(rr), n_bound=sum(r["bound"] for r in rr))
    return out


def print_rows(rows, kk, unit):
    print(f"   {'est':3s} {'curve':9s} {'exps':5s} {'m':>3s} {'cut':>4s}   l0 [{unit}] (16-84 % bootstrap)      z1    z2")
    for r in rows:
        ex = "-" if r["free"] is None else ("free" if r["free"] else "fixed")
        flag = "  lower bound" if r["bound"] else ""
        print(f"   {r['est']:3s} {r['curve']:9s} {ex:5s} {str(r['m'] or '-'):>3s} {r['cut']:>4g}   "
              f"{r['l0'] * kk:7.1f} ({r['lo'] * kk:6.1f}-{r['hi'] * kk:6.1f})"
              f"{'':6s}{r['z1']:5.2f} {r['z2']:5.2f}{flag}")


def half_point(l, y, plateau):
    """First l (log-interpolated) beyond the maximum where y drops below plateau / 2."""
    if not np.isfinite(plateau) or plateau <= 0:
        return np.nan
    k0 = int(np.nanargmax(y))
    for k in range(k0 + 1, len(y)):
        if y[k] < plateau / 2:
            t = (plateau / 2 - y[k - 1]) / (y[k] - y[k - 1])
            return float(np.exp(np.log(l[k - 1]) + t * np.log(l[k] / l[k - 1])))
    return np.nan


# ------------------------------------------------------------------ self test
if args.selftest:
    from scipy.ndimage import gaussian_filter
    from wall_tools import contour_points, pca_width, structure_factor

    rng = np.random.default_rng(7)
    H, W, nwall = 920, 700, 48
    ells = np.unique(np.round(np.logspace(np.log10(4), np.log10(H), 28)).astype(int)).astype(float)
    print("selftest: walls with S ~ q^-(1+2*0.5) [1 + (q l0/2pi)^2]^-(1.25-0.5), rms width 45 px, "
          f"{nwall} walls in 12 blocks; local width = PCA per segment on the contour")
    for l0_true in (100, 300):
        for blur in (0, 3):
            bs, bc, Sblk = np.zeros((12, ells.size)), np.zeros((12, ells.size)), []
            us = []
            for i in range(nwall):
                k = np.fft.rfftfreq(H)
                amp = np.zeros_like(k)
                amp[1:] = np.sqrt(k[1:] ** (-(1 + 2 * Z2)) * (1 + (k[1:] * l0_true) ** 2) ** (-(Z1 - Z2)))
                u = np.fft.irfft(amp * np.exp(2j * np.pi * rng.random(k.size)), n=H)
                u *= 45 / u.std()
                Y, X = np.mgrid[0:H, 0:W]
                dom = (X < W / 2 + u[:, None]).astype(float)
                if blur:
                    dom = gaussian_filter(dom, blur)
                domb = dom > 0.5
                yy, xx = contour_points(domb)
                w2l, _, _ = pca_width(yy, xx, 0.0, ells)
                b = i * 12 // nwall
                bs[b] += [v.sum() for v in w2l]; bc[b] += [v.size for v in w2l]
                ua = domb.sum(1).astype(float)
                us.append(ua - np.polyval(np.polyfit(np.arange(H), ua, 1), np.arange(H)))
            us = np.array(us)
            for b, idx in enumerate(np.array_split(np.arange(nwall), 12)):
                q, S = structure_factor(us[idx])
                Sblk.append(S)
            D = dict(ells=ells, blk_sum={"pca": bs}, blk_cnt={"pca": bc}, q=q, Sblk=np.array(Sblk), L=H, kk=1.0,
                     unit="px")
            rows, _ = analyse(D, 100, args.m, rng)
            s = summarize(rows)
            txt = ", ".join(f"{e}: {v['l0']:.0f} ({v['lo']:.0f}-{v['hi']:.0f}) [variants {v['sys_lo']:.0f}-{v['sys_hi']:.0f}]"
                            for e, v in s.items())
            free = [r for r in rows if r["est"] == "a" and r["free"] and r["m"] == args.m and r["cut"] == 25][0]
            print(f"   l0 = {l0_true} px, blur sigma {blur} px -> {txt}; free fit z1 = {free['z1']:.2f}, z2 = {free['z2']:.2f}")
    raise SystemExit


# ------------------------------------------------------------------ data
rng = np.random.default_rng(11)
results = {}
for d in args.dirs:
    lw = np.load(os.path.join(d, "local_width.npz"))
    rg = np.load(os.path.join(d, "roughness.npz"))
    kk = float(lw["um_per_px"]) or 1.0
    unit = "µm" if float(lw["um_per_px"]) else "px"
    pre = "area" if "area_q" in rg.files else "h"
    curves = [c for c in ("pca", "unfolded", "ctrl") if f"blk_sum_{c}" in lw.files]
    if not curves or f"{pre}_Sblk" not in rg.files:
        raise SystemExit(f"{d}: rerun step7 (per-block arrays missing)")
    ells = lw["ells"]
    # in --from-h studies the contour IS h and there are no folds: keep only one curve
    if np.allclose(lw["blk_sum_pca"], lw["blk_sum_unfolded"]) and np.allclose(lw["blk_sum_pca"], lw["blk_sum_ctrl"]):
        curves = ["pca"]
    names = {"pca": "contour", "unfolded": "unfolded", "ctrl": f"u_{pre}"}
    D = dict(ells=ells, blk_sum={names[c]: lw[f"blk_sum_{c}"] for c in curves},
             blk_cnt={names[c]: lw[f"blk_cnt_{c}"] for c in curves},
             q=rg[f"{pre}_q"], Sblk=rg[f"{pre}_Sblk"], L=float(ells.max()), kk=kk, unit=unit)
    rows, extra = analyse(D, args.n_boot, args.m, rng)
    summ = summarize(rows)
    ff, ts = lw["frac_folded"], lw["tilt_sd"]
    small = ells * kk <= (5 if unit == "µm" else 40)
    l_fold = half_point(ells, ff, ff[small].mean()) if ff[small].mean() > 0.01 else np.nan
    l_tilt = half_point(ells, ts, np.nanmax(ts))
    print(f"\n{d}: crossover length l0 (fits up to l = {extra['lmax'] * kk:.0f} {unit})")
    print_rows(rows, kk, unit)
    for e, lab in (("a", "local width, broken power law"), ("b", "global S(q), broken power law"),
                   ("c", "local slope of w2, midpoint")):
        if e in summ:
            v = summ[e]
            b = " (lower bound: not bracketed)" if v["bound"] else ""
            print(f"   [{e}] {lab:32s}: l0 = {v['l0'] * kk:6.1f} {unit} (bootstrap {v['lo'] * kk:.1f}-{v['hi'] * kk:.1f}; "
                  f"variants {v['sys_lo'] * kk:.1f}-{v['sys_hi'] * kk:.1f}, {v['n_bound']}/{v['n_var']} not bracketed){b}")
    print(f"   proxies: l_fold = {l_fold * kk:.1f} {unit}, l_tilt = {l_tilt * kk:.1f} {unit}")
    np.savez(os.path.join(d, "crossover.npz"), um_per_px=float(lw["um_per_px"]),
             table=np.array([(r["est"], r["curve"], str(r["free"]), str(r["m"]), r["cut"], r["l0"] * kk, r["lo"] * kk,
                              r["hi"] * kk, r["z1"], r["z2"], r["bound"]) for r in rows], dtype=object),
             **{f"l0_{e}": v["l0"] * kk for e, v in summ.items()}, **{f"l0_{e}_ci": (v["lo"] * kk, v["hi"] * kk) for e, v in summ.items()},
             **{f"l0_{e}_sys": (v["sys_lo"] * kk, v["sys_hi"] * kk) for e, v in summ.items()},
             **{f"l0_{e}_bound": v["bound"] for e, v in summ.items()}, l_fold=l_fold * kk, l_tilt=l_tilt * kk)
    results[d] = dict(rows=rows, summ=summ, extra=extra, kk=kk, unit=unit, l_fold=l_fold, l_tilt=l_tilt, lw=lw)

    # figure
    fig, ax = plt.subplots(1, 4, figsize=(22, 5))
    ex = extra
    lsel = ex["ells"] <= ex["lmax"]
    a_main = [r for r in rows if r["est"] == "a" and not r["free"] and r["m"] == args.m and r["cut"] == ex["cuts"][1]]
    a_free = [r for r in rows if r["est"] == "a" and r["free"] and r["m"] == args.m and r["cut"] == ex["cuts"][1]]
    for j, (c, w2) in enumerate(ex["w2"].items()):
        col = f"C{j}"
        x = ex["ells"][lsel]
        ax[0].loglog(x * kk, w2[lsel] / x ** (2 * Z1), "o", ms=3, color=col, label=f"{c}")
        for rr, ls in ((a_main, "-"), (a_free, "--")):
            r = [r_ for r_ in rr if r_["curve"] == c]
            if not r:
                continue
            r = r[0]
            sel = (x >= ex["cuts_px"][1])
            lnA = np.nanmedian(np.log(w2[lsel][sel]) - model_w2([0, r["z1"], r["z2"], np.log(r["l0"])], x[sel], args.m))
            ax[0].loglog(x * kk, np.exp(model_w2([lnA, r["z1"], r["z2"], np.log(r["l0"])], x, args.m)) / x ** (2 * Z1),
                         ls, color=col, lw=1)
    ax[0].axvline(ex["cuts_px"][1] * kk, color="0.6", ls=":", lw=0.8)
    ax[0].set_xlabel(f"ℓ [{unit}]"); ax[0].set_ylabel(r"$\langle w^2\rangle / \ell^{2.5}$")
    ax[0].set_title("Local width compensated by the qEW exponent\n(solid: z1 = 1.25, z2 = 0.5 fixed; dashed: free)", fontsize=10)
    ax[0].legend(fontsize=8)
    qb, S = ex["qb"], ex["S"]
    ax[1].loglog(qb / kk, S * qb ** (1 + 2 * Z1), "ko", ms=3, label="global S(q)")
    for free, ls in ((False, "-"), (True, "--")):
        r = [r_ for r_ in rows if r_["est"] == "b" and r_["free"] == free and r_["m"] == args.m and r_["cut"] == ex["cuts"][1]][0]
        sel = (qb <= 2 * np.pi / ex["cuts_px"][1])
        lnA = np.nanmedian(np.log(S[sel]) - model_S([0, r["z1"], r["z2"], np.log(r["l0"])], qb[sel], args.m))
        ax[1].loglog(qb / kk, np.exp(model_S([lnA, r["z1"], r["z2"], np.log(r["l0"])], qb, args.m)) * qb ** (1 + 2 * Z1),
                     ls, color="C3", lw=1, label=f"{'free' if free else 'fixed'}: l0 = {r['l0'] * kk:.1f} {unit}")
        ax[1].axvline(2 * np.pi / r["l0"] / kk, color="C3", ls=":", lw=0.8)
    ax[1].axvline(2 * np.pi / ex["cuts_px"][1] / kk, color="0.6", ls=":", lw=0.8)
    ax[1].set_xlabel(f"q [{unit}$^{{-1}}$]"); ax[1].set_ylabel(r"$S(q)\, q^{3.5}$")
    ax[1].set_title("Global S(q) compensated by the qEW exponent", fontsize=10); ax[1].legend(fontsize=8)
    for j, (c, w2) in enumerate(ex["w2"].items()):
        ax[2].semilogx(ex["ells"][lsel] * kk, local_slope(ex["ells"], w2)[lsel], "o-", ms=3, color=f"C{j}", label=c)
    for e, col in (("a", "C0"), ("b", "C3"), ("c", "C2")):
        if e in summ and np.isfinite(summ[e]["l0"]):
            v = summ[e]
            ax[2].axvspan(v["lo"] * kk, v["hi"] * kk, color=col, alpha=0.12)
            ax[2].axvline(v["l0"] * kk, color=col, lw=1, label=f"l0 ({e})")
    for zz in (Z1, Z2):
        ax[2].axhline(zz, color="0.8", lw=0.8, zorder=0)
    ax[2].set_xlabel(f"ℓ [{unit}]"); ax[2].set_ylabel("ζ_eff(ℓ) = ½ d ln w²/d ln ℓ")
    ax[2].set_title("Local slope of the local width, l0 estimates (band: bootstrap 16-84 %)", fontsize=10)
    ax[2].legend(fontsize=8)
    labs, vals = [], []
    for e, lab in (("a", "a: local width"), ("b", "b: global S(q)"), ("c", "c: local slope")):
        if e in summ:
            v = summ[e]
            labs.append(lab); vals.append(v)
    for i, v in enumerate(vals):
        ax[3].errorbar(v["l0"] * kk, i, xerr=[[max(v["l0"] - v["lo"], 0) * kk], [max(v["hi"] - v["l0"], 0) * kk]],
                       fmt="o" if not v["bound"] else ">", color="k", capsize=3)
        ax[3].plot([v["sys_lo"] * kk, v["sys_hi"] * kk], [i + 0.15] * 2, "-", color="0.6", lw=3, alpha=0.6)
    yy = len(vals)
    for val, lab in ((l_fold, "l_fold"), (l_tilt, "l_tilt")):
        if np.isfinite(val):
            ax[3].plot(val * kk, yy, "s", color="C1"); labs.append(lab); yy += 1
    ax[3].set_yticks(range(len(labs))); ax[3].set_yticklabels(labs); ax[3].set_xscale("log")
    ax[3].set_xlabel(f"length [{unit}]"); ax[3].axvline(ex["lmax"] * kk, color="0.6", ls=":", lw=0.8)
    ax[3].set_title("l0 estimates (black: main variant ± bootstrap; grey: spread over variants)\n"
                    "and disorder proxies (orange); dotted: largest ℓ", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(d, "fig_crossover.png"), dpi=110); plt.close(fig)
    print(f"   outputs: {d}/crossover.npz, fig_crossover.png")

if len(results) > 1:
    fig, ax = plt.subplots(figsize=(7, 1.2 + 0.9 * len(results) * 5 / 4))
    y = 0
    ticks, labs = [], []
    for d, R in results.items():
        kk = R["kk"]
        for e, lab, col in (("a", "local width", "C0"), ("b", "global S(q)", "C3"), ("c", "local slope", "C2")):
            if e in R["summ"]:
                v = R["summ"][e]
                ax.errorbar(v["l0"] * kk, y, xerr=[[max(v["l0"] - v["lo"], 0) * kk], [max(v["hi"] - v["l0"], 0) * kk]],
                            fmt="o" if not v["bound"] else ">", color=col, capsize=3)
                ax.plot([v["sys_lo"] * kk, v["sys_hi"] * kk], [y + 0.2] * 2, "-", color=col, lw=3, alpha=0.3)
                ticks.append(y); labs.append(f"{os.path.basename(d.rstrip('/'))}: {lab}"); y += 1
        for val, lab in ((R["l_fold"], "l_fold"), (R["l_tilt"], "l_tilt")):
            if np.isfinite(val):
                ax.plot(val * kk, y, "s", color="C1"); ticks.append(y); labs.append(f"{os.path.basename(d.rstrip('/'))}: {lab}"); y += 1
        y += 0.5
    ax.set_yticks(ticks); ax.set_yticklabels(labs, fontsize=8); ax.set_xscale("log"); ax.set_xlabel("length [µm]")
    ax.set_title("Crossover length l0 in the studies (> : lower bound)", fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(args.dirs[0], "fig_crossover_compare.png"), dpi=110); plt.close(fig)
    print(f"comparison: {args.dirs[0]}/fig_crossover_compare.png")
