"""Walls before/after each field pulse in the TbCo pulsed-field images (V. Jeudy and J. Rodriguez, LPS, Nov 2023, X20).

Usage: python tbco_pulses.py SRC_DIR PREFIX --outdir DIR --origin-csv CSV [--row-offset 2] [--max-row 51]
                             [--um-per-px 0.136] [--edge 30] [--sigma 2]

Raw 16-bit images (not difference images), one per pulse, files PREFIX_<n>.tif. The dodecagon is the microscope field
diaphragm (not a sample edge): only the diaphragm minus --edge px is analysed, and contour pieces on its border are
dropped; domains entering there come from outside the field of view.
Segmentation: frames are grouped by diaphragm setting (size and position of the illuminated field) and registered on the
sample features (defects, shifts <= 60 px); per group, each image is divided by an illumination factor chained between
consecutive images; bright reference B = 90th percentile per pixel; dark-state contrast map c(x) from the pixels that
switch in a pulse (1 - min/max of the before/after pair, fixed pattern cancels), smoothed (sigma 25 px), D = B (1 - c);
normalised signal n = (G_sigma(I) - D)/(B - D), wall = n = 0.5 (subpixel contour of the continuous field). Static
defects (dark spots in B) and areas with contrast < max(0.03, 0.4 median) are masked.
Pulses: image k-1 -> k for the rows of the Origin table (CSV from the worksheet PDFs) with a displacement, file = row +
--row-offset (rows up to --max-row). Per pulse:
 - swept area (|mask_k xor mask_k-1|) and the same from the consecutive difference (n_k - n_k-1 beyond +-0.5);
 - mean displacement = area / mean wall length; local displacement of each point of wall k-1 = distance to wall k
   (median, 90th percentile, max), to compare with the Origin hand-measured value (Criterio 3, a local forward advance);
 - local width w2(E) of arc windows rotated to their own axis (wall_tools.overhang_segments) on wall k-1;
 - KPZ term as step10 B on the contour: swept pixels assigned to the nearest point of wall k-1, segments of arc length L,
   v_n against the tilt dphi relative to the 8L window, v/<v> = a[1 + b dphi + kappa/2 dphi^2], Lambda = 1 + kappa.
Outputs in DIR: pulses.json (per pulse), walls.npz (contour pieces per file, px), roughness.npz, fig_overlay_*.png,
fig_pulses.png.
"""
import argparse
import csv
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.spatial import cKDTree
from skimage.filters import threshold_otsu
from skimage.measure import find_contours
from skimage.morphology import convex_hull_image, remove_small_holes, remove_small_objects

from wall_tools import fit_power, overhang_segments

ap = argparse.ArgumentParser()
ap.add_argument("src"); ap.add_argument("prefix")
ap.add_argument("--outdir", required=True); ap.add_argument("--origin-csv", required=True)
ap.add_argument("--row-offset", type=int, default=2); ap.add_argument("--max-row", type=int, default=10 ** 6)
ap.add_argument("--page-offsets", default="", help="e.g. '0:2,1:65': file = row + offset for rows on worksheet page block 0, 1")
ap.add_argument("--file-col", action="store_true", help="use the File column of the Origin table")
ap.add_argument("--skip-files", default="", help="comma list of files not used as pulses")
ap.add_argument("--um-per-px", type=float, default=0.136); ap.add_argument("--edge", type=int, default=30)
ap.add_argument("--sigma", type=float, default=2.0); ap.add_argument("--n-boot", type=int, default=300)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)
out = lambda f: os.path.join(a.outdir, f)
rng = np.random.default_rng(5)
STEP = 0.5

files = {int(f.rsplit("_", 1)[1][:-4]): f for f in glob.glob(os.path.join(a.src, a.prefix + "_*.tif"))
         if f.rsplit("_", 1)[1][:-4].isdigit()}
ks = sorted(files)
ld = lambda k: np.asarray(Image.open(files[k])).astype(np.float32)
G = lambda x: ndi.gaussian_filter(x, a.sigma)

# ---------------------------------------------------------------- diaphragm states, sample registration, references
raw = {k: ld(k) for k in ks}


def diaphragm(im, nray=360):
    """field diaphragm from its sharp edge. 4x-downsampled log image; along rays from the bright-field centre, the first
    local maximum of the radial drop d(log I)/dr above 30 % of the typical edge drop (a dim halo or vignetted rim beyond
    the edge cannot win over it); outlying radii replaced by the angular running median; convex hull (the diaphragm is
    a convex dodecagon). Candidates must be steps, not dips (dark scratches)."""
    g = ndi.gaussian_filter(im[::4, ::4].astype(float), 1.5)
    lg = np.log(np.maximum(g, np.percentile(g, 1)))
    gy, gx = ndi.sobel(lg, 0), ndi.sobel(lg, 1)
    lo, hi = np.percentile(g, 5), np.percentile(g, 99)
    cy, cx = ndi.center_of_mass(g > lo + 0.5 * (hi - lo))
    th = np.linspace(0, 2 * np.pi, nray, endpoint=False); r = np.arange(5, 400, 0.5)
    Y = cy + np.outer(np.sin(th), r); X = cx + np.outer(np.cos(th), r)
    ok = (Y >= 0) & (Y < g.shape[0] - 1) & (X >= 0) & (X < g.shape[1] - 1)
    dr = -(ndi.map_coordinates(gy, [Y, X], order=1) * np.sin(th)[:, None] + ndi.map_coordinates(gx, [Y, X], order=1) * np.cos(th)[:, None])
    dr[~ok] = -np.inf
    S = dr.max(1)
    lm = (dr >= np.roll(dr, 1, 1)) & (dr >= np.roll(dr, -1, 1)) & (dr > 0.3 * np.median(S))
    # a step, not a dip: the mean log level 3-20 px (downsampled) beyond the candidate stays 0.2 below the level just
    # inside it (a dark scratch or a domain edge recovers or drops less)
    P = ndi.map_coordinates(lg, [Y, X], order=1, mode="nearest")
    cs = np.c_[np.zeros(len(th)), np.cumsum(P, 1)]
    i_ = np.arange(r.size)
    win = lambda a_, b_: (cs[:, np.clip(i_ + b_, 0, r.size)] - cs[:, np.clip(i_ + a_, 0, r.size)]) / np.maximum(
        np.clip(i_ + b_, 0, r.size) - np.clip(i_ + a_, 0, r.size), 1)
    lm &= win(6, 40) < win(-8, -2) - 0.2
    R = np.where(lm.any(1), r[np.argmax(lm, 1)], r[np.argmax(dr, 1)])
    Rm = ndi.median_filter(R, 15, mode="wrap"); R = np.where(np.abs(R - Rm) > 6, Rm, R)
    pts = np.zeros(g.shape, bool)
    pts[np.clip(np.round(cy + R * np.sin(th)).astype(int), 0, g.shape[0] - 1),
        np.clip(np.round(cx + R * np.cos(th)).astype(int), 0, g.shape[1] - 1)] = True
    m = convex_hull_image(pts)
    lab, n = ndi.label(ndi.binary_fill_holes(m)); m = lab == 1 + np.argmax(ndi.sum(m, lab, range(1, n + 1)))
    return np.kron(m, np.ones((4, 4), bool))


DIA = {k: diaphragm(raw[k]) for k in ks}
state, reps = {}, []
for k in ks:                                     # same diaphragm state: overlap (IoU) > 0.9 with the state's first frame
    for i, r in enumerate(reps):
        if (DIA[k] & DIA[r]).sum() / (DIA[k] | DIA[r]).sum() > 0.9:
            state[k] = i; break
    else:
        reps.append(k); state[k] = len(reps) - 1
cnt_ = np.bincount([state[k] for k in ks], minlength=len(reps))
for k in ks:                                     # frames of a small state (detection jitter) join the best-overlapping
    if cnt_[state[k]] < 8:                       # large state if IoU > 0.75
        iou = [(DIA[k] & DIA[r]).sum() / (DIA[k] | DIA[r]).sum() if cnt_[i] >= 8 else 0 for i, r in enumerate(reps)]
        if max(iou) > 0.75:
            state[k] = int(np.argmax(iou))
bp = lambda x: ndi.gaussian_filter(x, 1) - ndi.gaussian_filter(x, 6)


def shift_to(ref, im, mask, smax=60):
    """sample shift of im relative to ref from band-pass features inside mask; argmax restricted to |s| <= smax.
    The camera fixed pattern (static in camera coordinates) always gives a peak at s = 0, so the sample peak is the
    maximum outside |s| <= 2, accepted if it reaches 0.25 of the s = 0 peak; shifts of 1-2 px are not resolved (otherwise s = 0)."""
    A = np.fft.fft2(bp(ref) * mask); Bf = np.fft.fft2(bp(im) * mask)
    c = np.fft.ifft2(A * np.conj(Bf) / np.maximum(np.abs(A * np.conj(Bf)), 1e-9)).real
    c = np.fft.fftshift(c); cy, cx = np.array(c.shape) // 2
    w = c[cy - smax:cy + smax + 1, cx - smax:cx + smax + 1].copy()
    c0 = w[smax - 2:smax + 3, smax - 2:smax + 3].max()
    w[smax - 2:smax + 3, smax - 2:smax + 3] = -np.inf
    iy, ix = np.unravel_index(np.argmax(w), w.shape)
    if w[iy, ix] < 0.25 * c0:                     # no sample peak away from the fixed-pattern peak: no shift
        return 0.0, 0.0
    return float(iy - smax), float(ix - smax)


SH = {}
for st_ in sorted(set(state.values())):
    mem = [k for k in ks if state[k] == st_]
    r0 = raw[mem[len(mem) // 2]]; mk = ndi.gaussian_filter(ndi.binary_erosion(DIA[mem[len(mem) // 2]], iterations=40).astype(np.float32), 10)
    for k in mem:
        SH[k] = shift_to(r0, raw[k], mk)
        if np.hypot(*SH[k]) > 0.5:
            raw[k] = ndi.shift(raw[k], SH[k], order=1, mode="nearest")
            DIA[k] = ndi.shift(DIA[k].astype(float), SH[k], order=0) > 0.5
STATES = {}
for st_ in sorted(set(state.values())):
    mem = [k for k in ks if state[k] == st_ and k != ks[0]]     # first file of each sample: as-grown state
    if len(mem) < 8:
        print(f"diaphragm state {st_}: files {mem} -- fewer than 8 frames, skipped", flush=True)
        continue
    dia_s = diaphragm(np.max([raw[k] for k in mem[:: max(1, len(mem) // 12)]], axis=0))   # brightest state: domains cannot cut it
    core0 = ndi.binary_erosion(dia_s, iterations=a.edge)
    Gs = {}
    scale = {mem[0]: 1.0}; prev = G(raw[mem[0]]); Gs[mem[0]] = prev
    for i, k in enumerate(mem[1:], 1):
        # illumination step = ratio of the UNCHANGED pixels: the histogram peak of log(g/prev) closest to 0 (a pulse
        # can switch most of the field, so the median may sit on the switched peak)
        g = G(raw[k]); lr = np.log((g / prev)[core0])
        hst, e = np.histogram(lr, bins=300, range=(-0.3, 0.3)); hst = ndi.gaussian_filter1d(hst.astype(float), 2)
        pk = [j for j in range(1, 299) if hst[j] >= hst[j - 1] and hst[j] >= hst[j + 1] and hst[j] > 0.05 * hst.max()]
        c0 = 0.5 * (e[1:] + e[:-1])[min(pk, key=lambda j: abs(e[j] + 0.001))] if pk else float(np.median(lr))
        sel_ = np.abs(lr - c0) < 0.015
        scale[k] = scale[mem[i - 1]] * float(np.exp(np.median(lr[sel_]) if sel_.sum() > 100 else c0)); prev = g; Gs[k] = g
    for k in mem:
        Gs[k] = Gs[k] / scale[k]
    # dark-state contrast from the pixels that switch in a pulse: 1 - min/max of the before/after pair (fixed pattern
    # cancels), pooled over all pulses of the state and smoothed (normalised convolution, sigma 25 px)
    csum = np.zeros(core0.shape, np.float32); cnt = np.zeros(core0.shape, np.float32)
    for i in range(1, len(mem)):
        k0_, k1_ = mem[i - 1], mem[i]
        d = Gs[k1_] - Gs[k0_]
        sd = 1.4826 * np.median(np.abs(d[core0] - np.median(d[core0])))
        sw = remove_small_objects((np.abs(d) > 6 * sd) & core0, 300)
        if sw.sum() < 500:
            continue
        lo_, hi_ = np.minimum(Gs[k0_], Gs[k1_]), np.maximum(Gs[k0_], Gs[k1_])
        csum[sw] += np.log(hi_[sw] / lo_[sw]); cnt[sw] += 1
    cmed = float(np.median(csum[cnt > 0] / cnt[cnt > 0])) if cnt.sum() > 0 else 0.0
    # per-pixel two-level reconstruction from the time series (no reference image): log signal L_k(x); a pixel whose
    # range over the state exceeds 0.6 of the local contrast has switched, and its own bright and dark levels are the
    # means of its values above / below the mid level (fixed pattern and vignetting included exactly). Smooth bright and
    # dark maps from the switched pixels (normalised convolution, sigma 25 px, 100 px fallback); a pixel that never
    # switched takes the level (bright or dark) closest to its median and the other one shifted by the local contrast;
    # unswitched pixels darker than dark (or brighter than bright) by more than the contrast are static defects.
    Lst = np.stack([np.log(np.maximum(Gs[k], 1.0)) for k in mem]).astype(np.float32); Gs = None
    num, den = ndi.gaussian_filter(csum, 25), ndi.gaussian_filter(cnt, 25)
    numL, denL = ndi.gaussian_filter(csum, 100), ndi.gaussian_filter(cnt, 100)
    ctil = np.where(den > 0.05, num / np.maximum(den, 1e-6), np.where(denL > 1e-3, numL / np.maximum(denL, 1e-6), cmed))
    ctil = np.clip(ctil, 0.0, 0.8).astype(np.float32)
    # levels and per-frame illumination are fitted together: the LED intensity jumps between frames by up to ~0.5 in log
    # and a pulse can reverse most of the field, so the chained ratio can be off; each frame is shifted by the median
    # residual of its pixels to their nearest level (bright or dark) and the levels re-estimated, 4 times
    sub = (slice(None), slice(0, None, 4), slice(0, None, 4)); c0s = core0[::4, ::4]
    for it_ in range(5):
        nf = Lst.shape[0]
        q = max(1, int(round(0.03 * nf)))
        srt = np.sort(Lst, axis=0)
        p_lo, p_hi, M = srt[min(q, nf - 1) if nf > 12 else 0], srt[max(nf - 1 - q, 0) if nf > 12 else nf - 1], srt[nf // 2]
        del srt
        swp = ((p_hi - p_lo) > 0.6 * ctil) & core0
        mid = 0.5 * (p_hi + p_lo)
        up = Lst > mid[None]
        nu = up.sum(0); Bl = np.where(nu > 0, (Lst * up).sum(0) / np.maximum(nu, 1), p_hi)
        Dl = np.where(nf - nu > 0, (Lst * ~up).sum(0) / np.maximum(nf - nu, 1), p_lo)
        del up
        def ncv(v, w, s_):
            return ndi.gaussian_filter(v * w, s_), ndi.gaussian_filter(w, s_)
        w_ = swp.astype(np.float32)
        nb_, db_ = ncv(np.where(swp, Bl, 0), w_, 25); nd_, dd_ = ncv(np.where(swp, Dl, 0), w_, 25)
        nbL, dbL = ncv(np.where(swp, Bl, 0), w_, 100); ndL, ddL = ncv(np.where(swp, Dl, 0), w_, 100)
        Mloc = ndi.gaussian_filter(M, 25)
        Bt = np.where(db_ > 0.05, nb_ / np.maximum(db_, 1e-6), np.where(dbL > 1e-3, nbL / np.maximum(dbL, 1e-6), Mloc))
        Dt = np.where(dd_ > 0.05, nd_ / np.maximum(dd_, 1e-6), np.where(ddL > 1e-3, ndL / np.maximum(ddL, 1e-6), Mloc - ctil))
        # the smooth maps carry the illumination at 25-100 px; the step between them is the local contrast
        cl = np.clip(Bt - Dt, 0.0, 0.8)
        bright = np.abs(M - Bt) < np.abs(M - Dt)
        B = np.where(swp, Bl, np.where(bright, M, M + cl)).astype(np.float32)
        D = np.where(swp, Dl, np.where(bright, M - cl, M)).astype(np.float32)
        if it_ == 4:
            break
        Bs, Ds = B[::4, ::4][c0s], D[::4, ::4][c0s]
        corr = np.array([np.median(np.where(np.abs(Li[c0s] - Bs) < np.abs(Li[c0s] - Ds), Li[c0s] - Bs, Li[c0s] - Ds)) for Li in Lst[sub]])
        Lst -= corr[:, None, None].astype(np.float32)
        for i, k in enumerate(mem):
            scale[k] *= float(np.exp(corr[i]))
        print(f"   illumination refit {it_}: |correction| max {np.abs(corr).max():.3f}, median {np.median(np.abs(corr)):.3f}; "
              f"switched {swp[core0].mean():.1%}", flush=True)
    defect = ~swp & core0 & ((M < Dt - cl) | (M > Bt + cl))
    defect = ndi.binary_dilation(remove_small_objects(defect, 20), iterations=6)
    lowc = cl < max(0.03, 0.4 * cmed)                                    # too little domain contrast to segment
    core = core0 & ~defect & ~lowc
    print(f"diaphragm state {st_}: {len(mem)} files ({mem[0]}-{mem[-1]}), field {core0.mean():.1%} of the image; switched "
          f"pixels {swp[core0].mean():.1%}; log contrast median {np.median(cl[core0]):.3f} (p5 {np.percentile(cl[core0], 5):.3f}); "
          f"defects {defect[core0].mean():.1%}, low contrast {(lowc & core0).sum() / max(core0.sum(), 1):.1%}; sample shifts up to "
          f"{max(np.hypot(*SH[k]) for k in mem):.1f} px", flush=True)
    Nst = {}
    den_ = np.maximum(B - D, 1e-3)
    for i, k in enumerate(mem):
        Nst[k] = ((Lst[i] - D) / den_).astype(np.float32)
    del Lst
    STATES[st_] = dict(mem=mem, core0=core0, core=core, dist_in=ndi.distance_transform_edt(core), N=Nst, scale=scale)
Nf, mask, CORE, CORE0, DIST = {}, {}, {}, {}, {}
for st_, S_ in STATES.items():
    for k in S_["mem"]:
        Nf[k] = S_["N"][k]
        mask[k] = remove_small_holes(remove_small_objects((Nf[k] < 0.5) & S_["core"], 300), 300)
        CORE[k], CORE0[k], DIST[k] = S_["core"], S_["core0"], S_["dist_in"]
    S_["N"] = None
scale = {k: S_["scale"][k] for S_ in STATES.values() for k in S_["mem"]}
del raw


def pieces(k, minlen=60):
    res = []
    for c in find_contours(np.where(CORE0[k], Nf[k], np.nan), 0.5):
        seg = np.hypot(*np.diff(c, axis=0).T); arc = np.r_[0, np.cumsum(seg)]
        if arc[-1] < minlen:
            continue
        s = np.arange(0, arc[-1], STEP); y, x = np.interp(s, arc, c[:, 0]), np.interp(s, arc, c[:, 1])
        yi, xi = np.clip(np.round(y).astype(int), 0, 2047), np.clip(np.round(x).astype(int), 0, 2047)
        ok = DIST[k][yi, xi] > 3
        d = np.diff(np.r_[0, ok.astype(int), 0]); st, en = np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]
        res += [(y[i:j], x[i:j]) for i, j in zip(st, en) if (j - i) * STEP >= minlen]
    return res


W = {k: pieces(k) for k in Nf}
Lw = {k: sum(p[0].size for p in W[k]) * STEP for k in W}

# ---------------------------------------------------------------- Origin rows -> pulses
org = {}
for r in csv.DictReader(open(a.origin_csv)):
    off, row = map(int, r["page_row"].split(":"))
    po = dict(tuple(map(int, t.split(":"))) for t in a.page_offsets.split(",")) if a.page_offsets else {0: a.row_offset}
    if a.file_col:
        if not r["file"].strip().isdigit():
            continue
        fk = int(r["file"])
    else:
        if off not in po or (off == 0 and row > a.max_row):
            continue
        fk = row + po[off]
    try:
        org[fk] = dict(H=float(r["H"]), dur_us=float(r["dur_us"]), disp_px=float(r["disp_px"]),
                                       err_px=float(r["err_px"]), V=float(r["V_HP"]))
    except ValueError:
        pass
skip = {int(t) for t in a.skip_files.split(",") if t}
pulses = [k for k in ks if k in org and k - 1 in Nf and k in Nf and state[k] == state[k - 1] and k not in skip
          and np.isfinite(org[k]["disp_px"])]
LAMS = np.unique(np.round(np.logspace(np.log10(4), np.log10(400), 20)))
LK = [16, 32, 64, 128]


def kpz_rows(k, L):
    """(dphi, v_n / <v_n>) per segment of arc length L of wall k-1 for pulse k-1 -> k."""
    sw = (mask[k] ^ mask[k - 1]) & CORE[k]
    py, px = np.nonzero(sw)
    rows = []
    n = int(round(L / STEP)); nb = 8 * n
    allp = [(y, x) for y, x in W[k - 1] if y.size >= nb + n]
    if not allp or py.size == 0:
        return np.empty((0, 2))
    Y = np.concatenate([p[0] for p in allp]); X = np.concatenate([p[1] for p in allp])
    pid = np.concatenate([np.full(p[0].size, i) for i, p in enumerate(allp)])
    dist, near = cKDTree(np.c_[Y, X]).query(np.c_[py, px], distance_upper_bound=4 * max(1.0, sw.sum() / max(Lw[k - 1], 1)) + 20)
    cnt = np.bincount(near[np.isfinite(dist)], minlength=Y.size)
    o = 0
    for i, (y, x) in enumerate(allp):
        c_ = cnt[o:o + y.size]; o += y.size
        for c0 in range(nb // 2, y.size - nb // 2 - n, n):
            sl = slice(c0, c0 + n); bl = slice(c0 + n // 2 - nb // 2, c0 + n // 2 + nb // 2)
            def ang(Yv, Xv):
                my, mx = Yv.mean(), Xv.mean()
                return 0.5 * np.arctan2(2 * ((Yv - my) * (Xv - mx)).mean(), ((Xv - mx) ** 2).mean() - ((Yv - my) ** 2).mean())
            d = np.angle(np.exp(2j * (ang(y[sl], x[sl]) - ang(y[bl], x[bl])))) / 2
            rows.append((d, c_[sl].sum() / (n * STEP)))
    r = np.array(rows)
    if r.size:
        r[:, 1] /= max(r[:, 1].mean(), 1e-12)
    return r


res = []
for k in pulses:
    o = org[k]
    core = CORE[k]
    sw = (mask[k] ^ mask[k - 1]) & core
    dn = ndi.gaussian_filter(Nf[k] - Nf[k - 1], 1)
    swd = remove_small_objects((dn > 0.5) & core, 300) | remove_small_objects((dn < -0.5) & core, 300)
    Lm = 0.5 * (Lw[k] + Lw[k - 1])
    loc = np.empty(0)
    if W[k] and W[k - 1]:
        tree = cKDTree(np.c_[np.concatenate([p[0] for p in W[k]]), np.concatenate([p[1] for p in W[k]])])
        loc = tree.query(np.c_[np.concatenate([p[0] for p in W[k - 1]]), np.concatenate([p[1] for p in W[k - 1]])])[0]
    w2 = [[] for _ in LAMS]
    for y, x in W[k - 1]:
        _, _, _, ww = overhang_segments(y, x, LAMS, STEP, [3.0])
        for j in range(LAMS.size):
            w2[j].append(np.asarray(ww[j]))
    rec = dict(file=k, H=o["H"], dur_us=o["dur_us"], V=o["V"], disp_origin=o["disp_px"], err_origin=o["err_px"],
               L_before=Lw[k - 1], L_after=Lw[k], n_pieces=len(W[k - 1]),
               A_mask=int(sw.sum()), A_diff=int(swd.sum()), overlap=float((sw & swd).sum() / max(1, (sw | swd).sum())),
               disp_mean=float(sw.sum() / Lm) if Lm > 50 else np.nan,
               loc_median=float(np.median(loc)) if loc.size else np.nan, loc_p90=float(np.percentile(loc, 90)) if loc.size else np.nan,
               loc_max=float(loc.max()) if loc.size else np.nan,
               w2_sum=[float(np.concatenate(v).sum()) if v else 0.0 for v in w2], w2_n=[int(sum(len(u) for u in v)) for v in w2],
               kpz={L: kpz_rows(k, L).tolist() for L in LK})
    res.append(rec)
    print(f"file {k:3d}  H {o['H']:6.1f} mT  {o['dur_us']:4.0f} us | wall {Lw[k-1]:6.0f} px in {len(W[k-1])} pieces | "
          f"mean disp {rec['disp_mean']:6.1f} px (area overlap mask/diff {rec['overlap']:.2f}) | local median {rec['loc_median']:6.1f} "
          f"p90 {rec['loc_p90']:6.1f} max {rec['loc_max']:6.1f} | Origin {o['disp_px']:5.0f} +- {o['err_px']:3.0f}", flush=True)

json.dump(dict(um_per_px=a.um_per_px, lams_px=LAMS.tolist(), L_kpz_px=LK, scale={str(k): v for k, v in scale.items()},
               state={str(k): int(v) for k, v in state.items()}, shift={str(k): v for k, v in SH.items()},
               pulses=res), open(out("pulses.json"), "w"))
np.savez_compressed(out("walls.npz"), **{f"f{k}_{i}": np.c_[y, x] for k in W for i, (y, x) in enumerate(W[k])})

# ---------------------------------------------------------------- figures
for blk in range(0, len(W), 14):
    sel = [k for k in ks if k in W][blk:blk + 14]
    fig, ax = plt.subplots(2, 7, figsize=(21, 6.6))
    for axx, k in zip(ax.flat, sel):
        im = ndi.shift(ld(k), SH[k], order=1)[::2, ::2] if np.hypot(*SH[k]) > 0.5 else ld(k)[::2, ::2]
        axx.imshow(im, cmap="gray", vmin=np.percentile(im, 2), vmax=np.percentile(im, 99.8))
        axx.contour(CORE[k][::2, ::2], [0.5], colors="y", linewidths=0.4)
        if k - 1 in W and state[k - 1] == state[k]:
            for y, x in W[k - 1]:
                axx.plot(x / 2, y / 2, "c-", lw=0.6)
        for y, x in W[k]:
            axx.plot(x / 2, y / 2, "r-", lw=0.7)
        t = f"file {k} s{state[k]}" + (f"  {org[k]['H']:.0f} mT {org[k]['dur_us']:.0f} us" if k in org else "")
        axx.set_title(t, fontsize=9); axx.axis("off")
    for axx in ax.flat[len(sel):]:
        axx.axis("off")
    plt.tight_layout(); plt.savefig(out(f"fig_overlay_{blk // 14}.png"), dpi=55); plt.close()
print("done", flush=True)
