"""Step 4: avalanche statistics from the arrival-time map (TIFF studies).

Usage: python step4_avalanches.py DIR [--smin 20] [--tau-m 1 2 5 10] [--roi X0 X1 Y0 Y1]
                                  [--compare FILE] [--um-per-px X] [--tag SUFFIX]
Input: DIR/tiff_extra.npz (arrival map, defects, rows) and DIR/meta.json, written by step1.
Outputs in DIR: avalanches<tag>_tm<m>.npz, fig_avalanches<tag>.png, fig_avalanche_map<tag>.png,
  summary on stdout.

An avalanche is a connected patch (8-neighbours) of pixels that switch within the same
measurement window: arrival frame in (j*m, (j+1)*m] for a window of m frames (--tau-m;
m = 1 is one frame). For every avalanche:
  S      area (px, and µm^2 with the scale)
  l_y    extent along y (approximately along the wall)
  l_x    extent along x (direction of motion)
  frame  last frame of its window; x, y: centroid (registered coordinates, rows as in the image)
Avalanches smaller than --smin px are discarded (default 20 px, the resolution limit used in
the thesis of this measurement, ~0.3 µm^2); the swept area they leave uncovered are the
"lagunas" of the thesis. Static defects (step1 inpainted their arrival time) are excluded.
Patches touching the image border, the right crop, the never-reached side or a defect are
kept but flagged (they may be cut).

Statistics:
  P(S), P(l_y): log-binned densities;
  tau, S_cut: maximum-likelihood fit of P(S) ~ S^-tau exp(-S/S_cut) above smin (the form used
         in the thesis), errors from the inverse Hessian. When S_cut is not much larger than
         smin there is no power-law range and tau is poorly determined;
  tau_S: pure power law (no cutoff) by maximum likelihood (Clauset et al. 2009), printed for
         reference only: with a cutoff in the data it overestimates the exponent;
  S ~ l_y^(1+zeta): power-law fit of the median S in log bins of l_y.
--compare FILE: a text file whose second column is the avalanche area in px (thesis format,
  Data/Histogramas_lineal/Avalanchas_15s.txt), overlaid on P(S); --roi restricts our
  avalanches to the same region (camera px, centroid inside X0<=x<X1, Y0<=y<Y1).
"""
import argparse
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import binary_dilation, center_of_mass, find_objects, label
from scipy.optimize import minimize

ap = argparse.ArgumentParser()
ap.add_argument("dir", help="step1 output folder of a TIFF study")
ap.add_argument("--smin", type=float, default=20, help="minimum avalanche area in px")
ap.add_argument("--tau-m", type=int, nargs="+", default=[1], help="measurement windows in frames")
ap.add_argument("--roi", type=float, nargs=4, metavar=("X0", "X1", "Y0", "Y1"), default=None,
                help="keep avalanches whose centroid is in this region (camera px)")
ap.add_argument("--compare", default=None, help="thesis-format avalanche list (col 2 = area in px)")
ap.add_argument("--um-per-px", type=float, default=None, help="default: the value in meta.json")
ap.add_argument("--tag", default="", help="suffix for the output files (e.g. _roi)")
args = ap.parse_args()

with open(os.path.join(args.dir, "meta.json")) as f:
    meta = json.load(f)
fps = meta["fps"]
k_um = meta.get("um_per_px") if args.um_per_px is None else args.um_per_px
e = np.load(os.path.join(args.dir, "tiff_extra.npz"))
ta, y_off = e["arrival"], int(e["rows"][0])
H, W = ta.shape
# static defects inside the swept area carry no magnetic contrast: their arrival time was
# inpainted in step1, so they are excluded here (not part of any avalanche nor of the lagunas)
defects = e["defects"] if "defects" in e.files else np.zeros(ta.shape, bool)
valid = np.isfinite(ta) & (ta >= 0) & ~defects
K = np.where(valid, np.round(ta), -1).astype(int)
never = np.isinf(ta)
T = int(K.max()) + 1
EIGHT = np.ones((3, 3), bool)
near_edge = np.zeros((H, W), bool)
near_edge[[0, -1], :] = near_edge[:, [0, -1]] = True
near_edge = binary_dilation(near_edge | never | defects, structure=EIGHT)


def avalanches(m):
    win = np.where(K >= 1, (K - 1) // m, -1)        # window index of each pixel
    rows = []
    for j in range(win.max() + 1):
        mask = win == j
        if not mask.any():
            continue
        lab, n = label(mask, structure=EIGHT)
        sizes = np.bincount(lab.ravel())[1:]
        keep = np.nonzero(sizes >= args.smin)[0] + 1
        if keep.size == 0:
            continue
        cms = center_of_mass(mask, lab, keep)
        slices = find_objects(lab)
        for i, (cy, cx) in zip(keep, cms):
            sl = slices[i - 1]
            ys, xs = sl
            # touches the image edges or the never-reached region -> may be cut
            border = bool((lab[sl] == i)[near_edge[sl]].any())
            rows.append((min((j + 1) * m, T - 1), sizes[i - 1], ys.stop - ys.start, xs.stop - xs.start,
                         cx, cy + y_off, border))
    a = np.array(rows, float).reshape(-1, 7)
    out = dict(frame=a[:, 0].astype(int), S=a[:, 1], l_y=a[:, 2], l_x=a[:, 3], x=a[:, 4], y=a[:, 5],
               edge=a[:, 6].astype(bool))
    if args.roi:
        x0, x1, y0, y1 = args.roi
        sel = (out["x"] >= x0) & (out["x"] < x1) & (out["y"] >= y0) & (out["y"] < y1)
        out = {k: v[sel] for k, v in out.items()}
    return out


def mle_tau(S, smin):
    """Continuous power-law MLE above smin (discrete areas: smin - 0.5 as lower bound)."""
    s = S[S >= smin]
    tau = 1 + s.size / np.sum(np.log(s / (smin - 0.5)))
    return tau, (tau - 1) / np.sqrt(s.size), s.size


def fit_cutoff(S, smin):
    """MLE of P(S) ~ S^-tau exp(-S/S_cut) above smin; errors from the inverse Hessian."""
    s = S[S >= smin]
    a = smin - 0.5

    def nll(p):
        tau, lsc = p
        sc = np.exp(lsc)
        x = np.logspace(np.log10(a), np.log10(a + 60 * sc), 4000)
        Z = np.trapz(x ** -tau * np.exp(-x / sc), x)
        return -(np.sum(-tau * np.log(s) - s / sc) - s.size * np.log(Z))

    best = min((minimize(nll, [t0, np.log(c0)], method="Nelder-Mead",
                         options=dict(xatol=1e-5, fatol=1e-6, maxiter=4000))
                for t0 in (0.5, 1.0, 1.5) for c0 in (2 * smin, 10 * smin, 50 * smin)), key=lambda r: r.fun)
    p0, h = best.x, np.array([1e-3, 1e-3])
    Hs = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            ei, ej = np.eye(2)[i] * h[i], np.eye(2)[j] * h[j]
            Hs[i, j] = (nll(p0 + ei + ej) - nll(p0 + ei - ej) - nll(p0 - ei + ej) + nll(p0 - ei - ej)) / (4 * h[i] * h[j])
    try:
        cov = np.linalg.inv(Hs)
        dtau, dlsc = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        dtau = dlsc = np.nan
    tau, sc = p0[0], np.exp(p0[1])
    return tau, dtau, sc, sc * dlsc


def logbin(v, vmin, nb=15):
    edges = np.logspace(np.log10(vmin), np.log10(v.max() * 1.0001), nb + 1)
    n, _ = np.histogram(v, edges)
    dens = n / np.diff(edges) / v.size
    c = np.sqrt(edges[1:] * edges[:-1])
    ok = n > 0
    return c[ok], dens[ok], n[ok]


area_unit = "µm²" if k_um else "px"
a2 = k_um ** 2 if k_um else 1.0
L_unit = "µm" if k_um else "px"
a1 = k_um if k_um else 1.0
results = {}
for m in args.tau_m:
    av = avalanches(m)
    np.savez(os.path.join(args.dir, f"avalanches{args.tag}_tm{m}.npz"), tau_m=m, smin=args.smin,
             um_per_px=k_um or 0.0, fps=fps, roi=args.roi or [], **av)
    tau, dtau, n = mle_tau(av["S"], args.smin)
    tau_in, dtau_in, n_in = mle_tau(av["S"][~av["edge"]], args.smin)
    tau_c, dtau_c, scut, dscut = fit_cutoff(av["S"], args.smin)
    # S ~ l_y^(1+zeta): median S in log bins of l_y
    ly, S = av["l_y"], av["S"]
    eb = np.unique(np.round(np.logspace(np.log10(max(ly.min(), 2)), np.log10(ly.max() + 1), 12)))
    med = [(np.sqrt(eb[i] * eb[i + 1]), np.median(S[(ly >= eb[i]) & (ly < eb[i + 1])]))
           for i in range(len(eb) - 1) if ((ly >= eb[i]) & (ly < eb[i + 1])).sum() >= 10]
    med = np.array(med)
    p = np.polyfit(np.log(med[:, 0]), np.log(med[:, 1]), 1) if len(med) >= 3 else [np.nan, np.nan]
    results[m] = dict(av=av, tau=tau, dtau=dtau, n=n, tau_in=tau_in, dtau_in=dtau_in, n_in=n_in,
                      tau_c=tau_c, dtau_c=dtau_c, scut=scut, dscut=dscut,
                      zeta=p[0] - 1, med=med, p=p)
    swept_mask = K >= 1
    if args.roi:
        x0, x1, ya, yb = args.roi
        yy, xx = np.mgrid[y_off:y_off + H, 0:W]
        swept_mask &= (xx >= x0) & (xx < x1) & (yy >= ya) & (yy < yb)
    swept = swept_mask.sum()
    print(f"tau_m = {m} frame(s) = {m / fps:g} s: {n} avalanches >= {args.smin:g} px "
          f"({av['edge'].sum()} touching an edge or a defect); they cover {av['S'].sum() / swept:.0%} "
          f"of the swept area (defects excluded); the rest, {1 - av['S'].sum() / swept:.0%}, are lagunas"
          + (" (approx., centroid in ROI)" if args.roi else ""))
    print(f"   S: median {np.median(S):.0f} px = {np.median(S) * a2:.2f} {area_unit}, "
          f"max {S.max():.0f} px = {S.max() * a2:.2f} {area_unit}; "
          f"l_y: median {np.median(ly):.0f} px, max {ly.max():.0f} px = {ly.max() * a1:.1f} {L_unit}")
    print(f"   tau_S (MLE, S >= {args.smin:g} px) = {tau:.3f} ± {dtau:.3f}; "
          f"without edge avalanches {tau_in:.3f} ± {dtau_in:.3f}")
    print(f"   with cutoff, P ~ S^-tau exp(-S/S_cut): tau = {tau_c:.2f} ± {dtau_c:.2f}, "
          f"S_cut = {scut:.0f} ± {dscut:.0f} px = {scut * a2:.2f} {area_unit}")
    print(f"   S ~ l_y^(1+zeta): 1+zeta = {p[0]:.2f} -> zeta = {p[0] - 1:.2f}")

# --- figures
ms = list(results)
cmap = plt.cm.viridis(np.linspace(0, 0.85, len(ms)))
fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
for c, m in zip(cmap, ms):
    r = results[m]
    x, d, _ = logbin(r["av"]["S"], args.smin)
    ax[0].loglog(x * a2, d / a2, "o", ms=4, color=c,
                 label=rf"$\tau_m$={m / fps:g} s: $\tau$={r['tau_c']:.2f}±{r['dtau_c']:.2f}, "
                       rf"$S_{{cut}}$={r['scut'] * a2:.2g} {area_unit}")
    xs = np.logspace(np.log10(args.smin), np.log10(r["av"]["S"].max()), 200)
    f = xs ** -r["tau_c"] * np.exp(-xs / r["scut"])
    xn = np.logspace(np.log10(args.smin - 0.5), np.log10(args.smin + 60 * r["scut"]), 4000)
    f /= np.trapz(xn ** -r["tau_c"] * np.exp(-xn / r["scut"]), xn)
    ax[0].loglog(xs * a2, f / a2, "-", lw=1, color=c)
    x, d, _ = logbin(r["av"]["l_y"], 1, 12)
    ax[1].loglog(x * a1, d / a1, "o-", ms=4, color=c, label=rf"$\tau_m$={m / fps:g} s")
    ax[2].loglog(r["av"]["l_y"] * a1, r["av"]["S"] * a2, ".", ms=1.5, color=c, alpha=0.25)
    if len(r["med"]):
        ax[2].loglog(r["med"][:, 0] * a1, r["med"][:, 1] * a2, "o", color=c, mec="k",
                     label=rf"$\tau_m$={m / fps:g} s: $1+\zeta$={r['p'][0]:.2f}")
        xx = np.array([r["med"][0, 0], r["med"][-1, 0]])
        ax[2].loglog(xx * a1, np.exp(np.polyval(r["p"], np.log(xx))) * a2, "-", color=c)
if args.compare:
    ref = np.loadtxt(args.compare)[:, 1]
    x, d, _ = logbin(ref, args.smin)
    tau_r, dtau_r, sc_r, dsc_r = fit_cutoff(ref, args.smin)
    ax[0].loglog(x * a2, d / a2, "s--", mfc="none", color="C3",
                 label=rf"comparison list: $\tau$={tau_r:.2f}±{dtau_r:.2f}, $S_{{cut}}$={sc_r * a2:.2g} {area_unit}")
    print(f"comparison list: {ref.size} avalanches, median S {np.median(ref):.0f} px; with cutoff "
          f"tau = {tau_r:.2f} ± {dtau_r:.2f}, S_cut = {sc_r:.0f} px")
ax[0].set_xlabel(f"S [{area_unit}]"); ax[0].set_ylabel("P(S)"); ax[0].legend(fontsize=7)
ax[0].set_title(rf"Avalanche areas (S ≥ {args.smin:g} px); lines: MLE $S^{{-\tau}}e^{{-S/S_{{cut}}}}$", fontsize=10)
ax[1].set_xlabel(rf"$\ell_y$ [{L_unit}]"); ax[1].set_ylabel(r"P($\ell_y$)"); ax[1].legend(fontsize=7)
ax[1].set_title("Lateral extent along the wall", fontsize=10)
ax[2].set_xlabel(rf"$\ell_y$ [{L_unit}]"); ax[2].set_ylabel(f"S [{area_unit}]"); ax[2].legend(fontsize=7)
ax[2].set_title(r"$S \sim \ell_y^{1+\zeta}$ (dots: avalanches; circles: median)", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(args.dir, f"fig_avalanches{args.tag}.png"), dpi=130); plt.close(fig)

# map of the avalanches of the first window, coloured by frame
m = ms[0]
win = np.where(K >= 1, (K - 1) // m, -1)
lab = np.zeros_like(K)
img = np.full(K.shape, np.nan)
for j in range(win.max() + 1):
    mask = win == j
    if not mask.any():
        continue
    l, n = label(mask, structure=EIGHT)
    sizes = np.bincount(l.ravel())
    big = sizes[l] >= args.smin
    img[mask & big] = j * m
fig, ax = plt.subplots(figsize=(9, 7.5))
ax.imshow(np.where(valid, 1.0, np.nan), cmap="gray", vmin=0, vmax=1.6, extent=(0, W, H + y_off, y_off))
ax.imshow(np.where(defects, 0.0, np.nan), cmap="gray", vmin=0, vmax=1, extent=(0, W, H + y_off, y_off),
          interpolation="nearest")
im = ax.imshow(img, cmap="jet", extent=(0, W, H + y_off, y_off), interpolation="nearest")
plt.colorbar(im, ax=ax, label="frame")
if args.roi:
    x0, x1, y0, y1 = args.roi
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], "k--", lw=1)
ax.set_title(f"Avalanches ≥ {args.smin:g} px (τ_m = {m / fps:g} s); grey: lagunas (smaller steps); "
             f"black: defects", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(args.dir, f"fig_avalanche_map{args.tag}.png"), dpi=110); plt.close(fig)
print(f"outputs in {args.dir}/")
