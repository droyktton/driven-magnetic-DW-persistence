"""Step 10: qEW or qKPZ below the crossover length l0? Tests on the measured walls.

Usage: python step10_universality.py DIR [DIR2 ...] [--n-boot 300] [--lmin 4]
Input (per DIR): local_width.npz, crossover.npz (step7, step8); for TIFF studies also contours.npz,
tiff_extra.npz (step7, step1) and avalanches_tm*.npz (step4); h_xt_sub.npy and meta.json.
Outputs: DIR/universality.npz, DIR/fig_universality.png, a table on stdout; with several DIRs also
fig_universality_compare.png in the first one.

Between the resolution and l0 the creep avalanches are depinning-like. The two candidate depinning
classes of a line in d = 2 differ in:
  A. roughness:   qEW zeta ~ 1.25,  qKPZ zeta ~ 0.63 (and the creep equilibrium regime 2/3);
  B. KPZ term:    the projected velocity of a segment with slope m relative to the mean motion,
                  v(m) = v0 [1 + (Lambda/2) m^2 + ...].  Lambda = lambda / v0.
                    Lambda = 0   elastic line driven along h (qEW, statistical tilt symmetry);
                    Lambda = 1   isotropic medium, the wall moves along its normal with a speed that
                                 does not depend on its orientation (geometric term, lambda = v -> 0);
                    Lambda >> 1  qKPZ (lambda stays finite when v -> 0, so lambda/v diverges).
                  Measured on the normal velocity v_n of contour segments of arc length L:
                  v_n(dphi) = a [1 + b dphi + (kappa/2) dphi^2], dphi = local normal - mean normal,
                  and since v = v_n / cos(dphi), Lambda = 1 + kappa.
                  R(L) = |Lambda| <m_L^2> / 2: relative velocity change produced by the tilts of
                  scale L; the KPZ term matters at scale L only if R(L) is not small.
  C. avalanche geometry: mean advance S/l_y ~ l_y^zeta for clusters with 3 um <= l_y < l0.
  D. up-down asymmetry: skewness of the local heights (from step7, sign: + = bumps ahead).
A. Free power law of the local width with local rotation (contour, folds excluded, u_area) over
   [--lmin, l0] (l0 from step8, local-width estimate), with lmin = 3, 4, 5 um as a check; and the
   broken power law of step8 (zeta2 = 0.5, m = 2, l0 free) over [lmin, L] with zeta1 fixed to 1.25,
   0.63 or 2/3, or free: chi2 per model. Errors: block bootstrap over time (step7 blocks).
B. TIFF studies: contour of every 5th frame (contours.npz, image coordinates), cut into segments of
   arc length L; orientation from PCA; normal towards the unswitched side. The pixels that switch in
   (t, t + dt] (arrival map) are assigned to the nearest contour point at t; v_n = their number /
   (L dt). Segments closer than --d-def px to a static defect, or near the image border, are
   dropped (also run without the defect mask). dt is chosen so that the mean advance is 2, 5 and
   15 px. Errors: block bootstrap over 12 time blocks. All studies: the same from h(x, t) (segments
   of L columns, v_n = v_x cos(phi)), which also covers video studies.
"""
import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view as swv
from scipy.ndimage import distance_transform_edt, binary_dilation
from scipy.optimize import least_squares
from scipy.spatial import cKDTree

ap = argparse.ArgumentParser()
ap.add_argument("dirs", nargs="+", help="study folders")
ap.add_argument("--n-boot", type=int, default=300, help="bootstrap resamples")
ap.add_argument("--lmin", type=float, default=4.0, help="lower end of the fit range in um (PSF bias, see step9)")
ap.add_argument("--d-def", type=float, default=8.0, help="drop segments closer than this to a defect, px")
args = ap.parse_args()
MODELS = {"qEW": 1.25, "qKPZ": 0.63, "equilibrium": 2 / 3}
rng = np.random.default_rng(10)


# ------------------------------------------------------------------ A. roughness exponent below l0
def model_w2(p, l, z1, z2=0.5, m=2.0):
    lnA, ll0 = p
    return lnA + 2 * z1 * np.log(l) + (2 * z2 - 2 * z1) / m * np.log1p((l / np.exp(ll0)) ** m)


def fit_broken(l, w2, sig, z1, lo, hi):
    """Broken power law with zeta2 = 0.5, m = 2; zeta1 fixed (number) or free (None). Returns chi2, l0, z1."""
    ok = np.isfinite(w2) & (w2 > 0)
    l, y, sig = l[ok], np.log(w2[ok]), np.maximum(sig[ok], 0.02)
    best = None
    for g in np.log(np.geomspace(lo * 1.2, hi / 1.2, 5)):
        if z1 is None:
            fun = lambda p: (model_w2(p[:2], l, p[2]) - y) / sig
            p0, lb, ub = [np.median(y - model_w2([0, g], l, 1.0)), g, 1.0], [-np.inf, np.log(lo), 0.0], [np.inf, np.log(hi), 3.0]
        else:
            fun = lambda p: (model_w2(p, l, z1) - y) / sig
            p0, lb, ub = [np.median(y - model_w2([0, g], l, z1)), g], [-np.inf, np.log(lo)], [np.inf, np.log(hi)]
        r = least_squares(fun, p0, bounds=(lb, ub))
        if best is None or r.cost < best.cost:
            best = r
    return 2 * best.cost, float(np.exp(best.x[1])), float(best.x[2]) if z1 is None else z1


def part_a(D):
    lw = np.load(os.path.join(D, "local_width.npz"))
    cr = np.load(os.path.join(D, "crossover.npz"))
    kk = float(lw["um_per_px"]) or 1.0
    unit = "µm" if float(lw["um_per_px"]) else "px"
    ells = lw["ells"].astype(float) * kk
    l0 = float(cr["l0_a"])                          # step8 stores l0 in the study unit
    l0_b = float(cr["l0_b"])
    curves = [c for c in ("pca", "unfolded", "ctrl") if f"blk_sum_{c}" in lw.files]
    nb = lw[f"blk_sum_{curves[0]}"].shape[0]
    picks = [np.arange(nb)] + [rng.integers(0, nb, nb) for _ in range(args.n_boot)]
    out = {"ells": ells, "l0": l0, "l0_b": l0_b, "unit": unit, "curves": curves}
    lmins = (3.0, 4.0, 5.0) if unit == "µm" else (24.0, 32.0, 40.0)
    lmin_main = args.lmin if unit == "µm" else 8 * args.lmin
    for c in curves:
        S, N = lw[f"blk_sum_{c}"], lw[f"blk_cnt_{c}"]
        W = np.array([np.where(N[p].sum(0) > 0, S[p].sum(0) / np.maximum(N[p].sum(0), 1), np.nan) for p in picks])
        sig = np.nanstd(np.log(W[1:]), 0)
        out[f"w2_{c}"], out[f"sig_{c}"] = W[0], sig
        for lm in lmins:
            for lmax_name, lmax in (("l0", l0), ("l0min", min(l0, l0_b))):
                sel = (ells >= lm) & (ells <= lmax)
                if sel.sum() < 3:
                    out[f"zfree_{c}_{lm:g}_{lmax_name}"] = (np.nan, np.nan, np.nan, int(sel.sum()))
                    continue
                z = np.array([np.polyfit(np.log(ells[sel]), np.log(w[sel]), 1)[0] / 2 for w in W])
                out[f"zfree_{c}_{lm:g}_{lmax_name}"] = (z[0], *np.nanpercentile(z[1:], [16, 84]), int(sel.sum()))
        sel = ells >= lmin_main
        Lmax = ells[np.isfinite(W[0])].max()
        chis = {name: fit_broken(ells[sel], W[0][sel], sig[sel], z1, lmin_main, 4 * Lmax) for name, z1 in MODELS.items()}
        chis["free"] = fit_broken(ells[sel], W[0][sel], sig[sel], None, lmin_main, 4 * Lmax)
        zb = [fit_broken(ells[sel], w[sel], sig[sel], None, lmin_main, 4 * Lmax)[2] for w in W[1:101]]
        out[f"broken_{c}"] = chis
        out[f"broken_zfree_ci_{c}"] = np.nanpercentile(zb, [16, 84])
        out[f"broken_npts_{c}"] = int((sel & np.isfinite(W[0])).sum())
    return out


# ------------------------------------------------------------------ B. KPZ term: v_n vs orientation
def fit_quadratic(dphi, v, blocks):
    """v/<v> = a + a b dphi + a kappa/2 dphi^2 by least squares; Lambda = 1 + kappa. Block bootstrap."""
    def one(idx):
        d, y = dphi[idx], v[idx] / v[idx].mean()
        A = np.vstack([np.ones_like(d), d, d * d]).T
        c = np.linalg.lstsq(A, y, rcond=None)[0]
        return c[1] / c[0], 2 * c[2] / c[0]
    b0, k0 = one(np.concatenate(blocks))
    bs = []
    for _ in range(args.n_boot):
        idx = np.concatenate([blocks[j] for j in rng.integers(0, len(blocks), len(blocks))])
        bs.append(one(idx))
    bs = np.array(bs)
    return dict(b=b0, kappa=k0, Lambda=1 + k0, kappa_ci=np.percentile(bs[:, 1], [16, 84]),
                b_ci=np.percentile(bs[:, 0], [16, 84]))


def binned_curve(dphi, v, edges):
    i = np.digitize(dphi, edges) - 1
    m = np.array([v[i == k].mean() if (i == k).sum() > 20 else np.nan for k in range(len(edges) - 1)])
    e = np.array([v[i == k].std() / np.sqrt((i == k).sum()) if (i == k).sum() > 20 else np.nan for k in range(len(edges) - 1)])
    return m / np.nanmean(v), e / np.nanmean(v)


def segments_contour(D, Ls_px, dts, mask_defects):
    c = np.load(os.path.join(D, "contours.npz"))
    e = np.load(os.path.join(D, "tiff_extra.npz"))
    ta, defects = e["arrival"].astype(float), e["defects"]
    H, W = ta.shape
    ddef = distance_transform_edt(~binary_dilation(defects, iterations=1)) if defects.any() else np.full(ta.shape, np.inf)
    ys, xs, off, frames = c["y"], c["x"], c["offsets"], c["frames"]
    step = float(c["step"])
    res = {(L, dt): ([], [], []) for L in Ls_px for dt in dts}       # dphi, v_n, frame
    for i, t in enumerate(frames):
        y, x = ys[off[i]:off[i + 1]], xs[off[i]:off[i + 1]]
        if y.size < 50:
            continue
        tree = cKDTree(np.c_[y, x])
        for dt in dts:
            if t + dt > frames[-1]:
                continue
            sw = (ta > t) & (ta <= t + dt)
            py, px = np.nonzero(sw)
            if py.size == 0:
                near = np.empty(0, int)
            else:
                dist, near = tree.query(np.c_[py, px], distance_upper_bound=60)
                near = near[np.isfinite(dist)]
            cnt = np.bincount(near, minlength=y.size)
            for L in Ls_px:
                n = int(round(L / step))
                ns = y.size // n
                if ns < 2:
                    continue
                Y, X, C = y[:ns * n].reshape(ns, n), x[:ns * n].reshape(ns, n), cnt[:ns * n].reshape(ns, n)
                my, mx = Y.mean(1), X.mean(1)
                syy, sxx = ((Y - my[:, None]) ** 2).mean(1), ((X - mx[:, None]) ** 2).mean(1)
                sxy = ((Y - my[:, None]) * (X - mx[:, None])).mean(1)
                ang = 0.5 * np.arctan2(2 * sxy, sxx - syy)            # tangent direction (x towards y)
                nx, ny = -np.sin(ang), np.cos(ang)                      # a normal
                # normal towards the unswitched side: probe 4 px along it at the segment centre
                cy = np.clip(np.round(my + 4 * ny).astype(int), 0, H - 1)
                cx = np.clip(np.round(mx + 4 * nx).astype(int), 0, W - 1)
                flip = ta[cy, cx] <= t
                nx, ny = np.where(flip, -nx, nx), np.where(flip, -ny, ny)
                phi = np.arctan2(ny, nx)
                ok = (Y.min(1) > 6) & (Y.max(1) < H - 7) & (X.min(1) > 6) & (X.max(1) < W - 7)
                if mask_defects:
                    yi = np.clip(np.round(Y).astype(int), 0, H - 1)
                    xi = np.clip(np.round(X).astype(int), 0, W - 1)
                    ok &= ddef[yi, xi].min(1) > args.d_def
                v = C.sum(1) / (n * step) / dt
                r = res[(L, dt)]
                r[0].append(phi[ok]); r[1].append(v[ok]); r[2].append(np.full(ok.sum(), t))
    return {k: tuple(np.concatenate(a) for a in v) for k, v in res.items() if v[0]}


def segments_h(D, Ls_px, dts):
    """From h(x, t): advance = -dh (both conventions of step1). Segments of L columns."""
    h = np.load(os.path.join(D, "h_xt_sub.npy"))
    T = h.shape[0]
    res = {}
    for dt in dts:
        for L in Ls_px:
            L = int(L)
            y = np.arange(L) - (L - 1) / 2
            P, V, F = [], [], []
            for t in range(0, T - dt, max(dt // 2, 1)):
                w0, w1 = swv(-h[t], L)[::max(L // 2, 1)], swv(-h[t + dt], L)[::max(L // 2, 1)]
                ok = np.isfinite(w0).all(1) & np.isfinite(w1).all(1)
                m = ((w0 - w0.mean(1, keepdims=True)) * y).sum(1) / (y * y).sum()
                vx = (w1 - w0).mean(1) / dt
                phi = -np.arctan(m)                    # normal angle (sign irrelevant for kappa)
                P.append(phi[ok]); V.append((vx * np.cos(phi))[ok]); F.append(np.full(ok.sum(), t))
            res[(L, dt)] = tuple(np.concatenate(a) for a in (P, V, F))
    return res


def analyse_b(res, kk):
    out = {}
    for (L, dt), (phi, v, fr) in res.items():
        if phi.size < 200:
            continue
        dphi = np.angle(np.exp(1j * (phi - np.angle(np.mean(np.exp(1j * phi))))))
        keep = np.abs(dphi) < np.radians(60)
        dphi, v, fr = dphi[keep], v[keep], fr[keep]
        order = np.argsort(fr)
        blocks = np.array_split(order, 12)
        f = fit_quadratic(dphi, v, blocks)
        m2 = np.mean(np.tan(dphi) ** 2)
        f.update(R=abs(f["Lambda"]) * m2 / 2, m2=m2, n=int(v.size), vmean=float(v.mean()), L_um=L * kk, dt=dt)
        f["curve"] = binned_curve(dphi, v, np.radians(np.arange(-60, 61, 10)))
        out[(L, dt)] = f
    return out


# ------------------------------------------------------------------ C. avalanche geometry
def part_c(D, kk, l0):
    out = {}
    for tm in (3, 10, 30):
        f = os.path.join(D, f"avalanches_tm{tm}.npz")
        if not os.path.exists(f):
            continue
        a = np.load(f)
        ly, S = a["l_y"] * kk, a["S"] * kk * kk
        sel = (ly >= 3.0) & (ly < l0) & ~a["edge"].astype(bool)
        if sel.sum() < 30:
            continue
        x, y = np.log(ly[sel]), np.log(S[sel] / ly[sel])
        z0 = np.polyfit(x, y, 1)[0]
        zb = [np.polyfit(x[i], y[i], 1)[0] for i in (rng.integers(0, x.size, x.size) for _ in range(args.n_boot))]
        out[tm] = (z0, *np.percentile(zb, [16, 84]), int(sel.sum()))
    return out


# ------------------------------------------------------------------ main
summary = []
for D in args.dirs:
    meta = json.load(open(os.path.join(D, "meta.json")))
    kk = meta.get("um_per_px") or 1.0
    tiff = os.path.exists(os.path.join(D, "contours.npz"))
    A = part_a(D)
    l0 = A["l0"]
    print(f"\n=== {D}: l0 (step8, local width) = {l0:.1f} {A['unit']}, (global S) = {A['l0_b']:.1f}")
    print("A. zeta of the local width (rotation per segment), free power law; 16-84 % block bootstrap")
    for c in A["curves"]:
        name = {"pca": "contour", "unfolded": "folds excluded", "ctrl": "u_area / h"}[c]
        line = []
        for k in sorted(k for k in A if k.startswith(f"zfree_{c}_")):
            z, lo, hi, n = A[k]
            line.append(f"[{k.split('_')[-2]}–{k.split('_')[-1]}] {z:.2f} ({lo:.2f}–{hi:.2f}, {n} pts)")
        print(f"   {name:15s} " + "  ".join(line))
    print(f"A. broken power law over [{args.lmin:g} {A['unit']}, L] (zeta2 = 0.5, m = 2, l0 free): chi2 per model")
    for c in A["curves"]:
        ch = A[f"broken_{c}"]
        lo, hi = A[f"broken_zfree_ci_{c}"]
        print(f"   {c:9s} n={A[f'broken_npts_{c}']:2d}  " + "  ".join(f"{k}: chi2={v[0]:7.1f} (l0={v[1]:5.1f})" for k, v in ch.items() if k != "free")
              + f"  free: chi2={ch['free'][0]:6.1f} zeta1={ch['free'][2]:.2f} ({lo:.2f}–{hi:.2f}) l0={ch['free'][1]:.1f}")
    # B
    Ls_um = [l for l in (2.0, 4.0, 8.0, 16.0) if l < l0] if kk != 1.0 else [16, 32, 64, 128]
    Ls_px = [l / kk for l in Ls_um] if kk != 1.0 else Ls_um
    h = np.load(os.path.join(D, "h_xt_sub.npy"))
    adv = np.nanmean(h[0]) - np.nanmean(h[-1])
    v_mean = abs(adv) / (h.shape[0] - 1)
    B = {}
    if tiff:
        dts = sorted({max(5, int(round(a / v_mean / 5)) * 5) for a in (2, 5, 15)})
        B["contour"] = analyse_b(segments_contour(D, Ls_px, dts, True), kk)
        B["contour, no defect mask"] = analyse_b(segments_contour(D, Ls_px, dts, False), kk)
    dts_h = sorted({max(1, int(round(a / v_mean))) for a in (2, 5, 15)})
    B["h(x)"] = analyse_b(segments_h(D, [max(4, int(round(l))) for l in Ls_px], dts_h), kk)
    print("B. normal velocity vs local orientation: Lambda = lambda/v (0 elastic qEW, 1 isotropic, >>1 qKPZ),"
          " R = |Lambda|<m_L^2>/2")
    for meth, res in B.items():
        for (L, dt), f in sorted(res.items()):
            print(f"   {meth:24s} L={f['L_um']:5.1f} dt={dt:4d}  n={f['n']:6d}  Lambda={f['Lambda']:+.2f} "
                  f"({1 + f['kappa_ci'][0]:+.2f}..{1 + f['kappa_ci'][1]:+.2f})  b={f['b']:+.2f}  "
                  f"<m^2>={f['m2']:.2f}  R={f['R']:.2f}")
    C = part_c(D, kk, l0) if tiff and kk != 1.0 else {}
    for tm, (z, lo, hi, n) in C.items():
        print(f"C. avalanches tau_m={tm:3d}: S/l_y ~ l_y^zeta, 3 µm <= l_y < l0: zeta = {z:.2f} ({lo:.2f}–{hi:.2f}), n = {n}")
    hd = os.path.join(D, "height_distribution.npz")
    Dsk = {}
    if os.path.exists(hd):
        z = np.load(hd)
        sgn = 1.0 if tiff else -1.0              # u_area grows with the advance; h of video studies decreases
        for k in z.files:
            if k.endswith("|skew") and "local" in k:
                Dsk[k[:-5]] = (sgn * float(z[k]), float(z[k[:-5] + "|dskew"]))
        print("D. skewness of local heights (+ = bumps ahead): " + "  ".join(f"{k}: {v:+.2f}±{e:.2f}" for k, (v, e) in Dsk.items()))

    # figure
    fig, ax = plt.subplots(1, 4, figsize=(22, 5))
    e = A["ells"]
    for c, col in zip(A["curves"], ("C0", "C3", "C1")):
        ax[0].errorbar(e, A[f"w2_{c}"] * (kk if A["unit"] == "µm" else 1) ** 2, A[f"w2_{c}"] * A[f"sig_{c}"] * (kk if A["unit"] == "µm" else 1) ** 2,
                       fmt="o", ms=3, color=col, label={"pca": "contour", "unfolded": "folds excluded", "ctrl": "u_area / h"}[c])
    ref = A["w2_unfolded"] if "unfolded" in A["curves"] else A[f"w2_{A['curves'][0]}"]
    j = np.nanargmin(np.abs(e - (args.lmin if A["unit"] == "µm" else 32)))
    for (name, z1), ls in zip(MODELS.items(), ("-", "--", ":")):
        ax[0].plot(e, ref[j] * (kk if A["unit"] == "µm" else 1) ** 2 * (e / e[j]) ** (2 * z1), "k" + ls, lw=1, label=f"{name}: ζ = {z1:.2f}")
    ax[0].axvline(l0, color="g", lw=0.8); ax[0].axvline(args.lmin if A["unit"] == "µm" else 32, color="gray", lw=0.8, ls="--")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel(f"ℓ [{A['unit']}]"); ax[0].set_ylabel("w²")
    ax[0].set_title("A. local width (rotation per segment); green: l0"); ax[0].legend(fontsize=7)
    for c, col in zip(A["curves"], ("C0", "C3", "C1")):
        ks = sorted(k for k in A if k.startswith(f"zfree_{c}_") and k.endswith("_l0"))
        lm = [float(k.split("_")[-2]) for k in ks]
        zz = np.array([A[k][:3] for k in ks])
        ax[1].errorbar(lm, zz[:, 0], [zz[:, 0] - zz[:, 1], zz[:, 2] - zz[:, 0]], fmt="o-", color=col, capsize=3)
    for (name, z1), ls in zip(MODELS.items(), ("-", "--", ":")):
        ax[1].axhline(z1, color="k", ls=ls, lw=0.8)
    ax[1].set_xlabel(f"lower end of the fit [{A['unit']}] (upper end l0)"); ax[1].set_ylabel("ζ"); ax[1].set_title("A. ζ below l0")
    ph = np.arange(-55, 56, 10)
    meth = "contour" if "contour" in B else "h(x)"
    for (L, dt), f in sorted(B[meth].items()):
        if dt == sorted({k[1] for k in B[meth]})[len({k[1] for k in B[meth]}) // 2]:
            ax[2].errorbar(ph, f["curve"][0], f["curve"][1], fmt="o-", ms=3, label=f"L = {f['L_um']:.0f}")
    g = np.radians(ph)
    ax[2].plot(ph, np.ones_like(ph), "k-", lw=0.8, label="isotropic (Λ = 1)")
    ax[2].plot(ph, np.cos(g) / np.mean(np.cos(g)), "k--", lw=0.8, label="elastic along h (Λ = 0)")
    ax[2].set_xlabel("local normal − mean normal [deg]"); ax[2].set_ylabel("v_n / ⟨v_n⟩")
    ax[2].set_title(f"B. normal velocity vs orientation ({meth}, middle dt)"); ax[2].legend(fontsize=7)
    for (mname, res), mk in zip(B.items(), ("o", "s", "^")):
        for dt in sorted({k[1] for k in res}):
            fs = [res[k] for k in sorted(res) if k[1] == dt]
            ax[3].errorbar([f["L_um"] for f in fs], [f["Lambda"] for f in fs],
                           [[f["Lambda"] - 1 - f["kappa_ci"][0] for f in fs], [1 + f["kappa_ci"][1] - f["Lambda"] for f in fs]],
                           fmt=mk + "-", ms=4, capsize=2, label=f"{mname}, dt={dt}")
    ax[3].axhline(1, color="k", lw=0.8); ax[3].axhline(0, color="k", lw=0.8, ls="--")
    ax[3].set_xscale("log"); ax[3].set_xlabel(f"segment length L [{A['unit']}]"); ax[3].set_ylabel("Λ = λ / v")
    ax[3].set_title("B. KPZ coefficient (0 elastic, 1 isotropic, ≫1 qKPZ)"); ax[3].legend(fontsize=6)
    fig.suptitle(D); fig.tight_layout(); fig.savefig(os.path.join(D, "fig_universality.png"), dpi=100); plt.close(fig)
    np.savez(os.path.join(D, "universality.npz"), A=np.array(A, dtype=object), B=np.array(B, dtype=object),
             C=np.array(C, dtype=object), D=np.array(Dsk, dtype=object), lmin=args.lmin, d_def=args.d_def)
    summary.append((D, A, B, C))
    print(f"outputs: {D}/universality.npz, fig_universality.png")

if len(summary) > 1:
    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    for i, (D, A, B, C) in enumerate(summary):
        c = "unfolded" if "unfolded" in A["curves"] else A["curves"][0]
        z, lo, hi, _ = A[f"zfree_{c}_{args.lmin:g}_l0"] if A["unit"] == "µm" else A[f"zfree_{c}_32_l0"]
        ax[0].errorbar(i, z, [[z - lo], [hi - z]], fmt="o", capsize=3)
        meth = "contour" if "contour" in B else "h(x)"
        Ls = [f["Lambda"] for f in B[meth].values()]
        ax[1].plot(np.full(len(Ls), i), Ls, "o", alpha=0.6)
    for (name, z1), ls in zip(MODELS.items(), ("-", "--", ":")):
        ax[0].axhline(z1, color="k", ls=ls, lw=0.8, label=name)
    ax[1].axhline(1, color="k", lw=0.8); ax[1].axhline(0, color="k", lw=0.8, ls="--")
    for a in ax:
        a.set_xticks(range(len(summary))); a.set_xticklabels([s[0] for s in summary], rotation=20, fontsize=7)
    ax[0].set_ylabel(f"ζ below l0 (folds excluded, from {args.lmin:g} µm)"); ax[0].legend(fontsize=7)
    ax[1].set_ylabel("Λ = λ / v (all L, dt)")
    fig.tight_layout(); fig.savefig(os.path.join(args.dirs[0], "fig_universality_compare.png"), dpi=100); plt.close(fig)
    print(f"comparison: {args.dirs[0]}/fig_universality_compare.png")
