"""Shared helpers for the roughness steps (step7_roughness.py, step8_crossover.py): detrending,
structure factor, height differences, local width, log binning and power-law fits, rotation onto the
base plane, subpixel contour, local width with local rotation (PCA per segment) and spectra of
rotated segments. See step7_roughness.py for the definitions of the observables."""
import numpy as np
from scipy.ndimage import rotate


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
    best line, i.e. after rotating the segment by its local tilt). Returns per-l lists of w2,
    local tilt (deg, relative to theta0) and contour length inside the window divided by l."""
    t0 = np.radians(theta0_deg)
    s = y * np.cos(t0) + x * np.sin(t0)
    o = np.argsort(s)
    s, y, x = s[o], y[o], x[o]
    C = [np.concatenate([[0], np.cumsum(v)]) for v in (np.ones_like(y), y, x, y * y, x * x, x * y)]
    out_w2, out_tilt, out_arc = [], [], []
    for l in ells:
        s0 = np.arange(s[0], s[-1] - l, l / 2)
        a = np.searchsorted(s, s0)
        b = np.searchsorted(s, s0 + l)
        n = b - a
        ok = n >= max(4, 0.5 * l / step)
        a, b, n = a[ok], b[ok], n[ok].astype(float)
        if n.size == 0:
            out_w2.append(np.empty(0)); out_tilt.append(np.empty(0)); out_arc.append(np.empty(0)); continue
        my, mx = (C[1][b] - C[1][a]) / n, (C[2][b] - C[2][a]) / n
        syy = (C[3][b] - C[3][a]) / n - my ** 2
        sxx = (C[4][b] - C[4][a]) / n - mx ** 2
        sxy = (C[5][b] - C[5][a]) / n - mx * my
        lam = (syy + sxx) / 2 - np.sqrt(((syy - sxx) / 2) ** 2 + sxy ** 2)
        alpha = 0.5 * np.arctan2(2 * sxy, syy - sxx)          # major axis, from y towards x
        out_w2.append(np.clip(lam, 0, None)); out_tilt.append(np.degrees(alpha - t0)); out_arc.append(n * step / l)
    return out_w2, out_tilt, out_arc


def seg_spectra(y, x, theta0_deg, l, step, edges, dx=1.0):
    """Spectra of the segments of projected length l (as in pca_width), each on its own axis.
    Returns per-segment arrays: binned S sums (nseg x nbins), counts, contour length / l."""
    t0 = np.radians(theta0_deg)
    s = y * np.cos(t0) + x * np.sin(t0)
    o = np.argsort(s)
    s, y, x = s[o], y[o], x[o]
    nb_q = len(edges) - 1
    out_S, out_C, out_arc = [], [], []
    for s0 in np.arange(s[0], s[-1] - l, l / 2):
        a0, b0 = np.searchsorted(s, s0), np.searchsorted(s, s0 + l)
        if b0 - a0 < max(4, 0.5 * l / step):
            continue
        py, px = y[a0:b0], x[a0:b0]
        my, mx = py.mean(), px.mean()
        cyy, cxx, cxy = ((py - my) ** 2).mean(), ((px - mx) ** 2).mean(), ((py - my) * (px - mx)).mean()
        al = 0.5 * np.arctan2(2 * cxy, cyy - cxx)                       # own axis, from y towards x
        ca, sa = np.cos(al), np.sin(al)
        ax_ = (py - my) * ca + (px - mx) * sa
        nn = -(py - my) * sa + (px - mx) * ca
        ib = np.floor((ax_ - ax_.min()) / dx).astype(int)
        nb = ib.max() + 1
        if nb < 16:
            continue
        cnt = np.bincount(ib, minlength=nb)
        prof = np.bincount(ib, nn, minlength=nb) / np.maximum(cnt, 1)
        good = cnt > 0
        prof = np.interp(np.arange(nb), np.nonzero(good)[0], prof[good])
        # line through the end points (made periodic): unbiased from the first mode for short segments,
        # whereas a least-squares line + Hann window overestimates zeta by ~0.05-0.2 when l < 400 px
        e0, e1 = prof[:3].mean(), prof[-3:].mean()
        prof = prof - (e0 + (e1 - e0) * (np.arange(nb) - 1) / max(nb - 3, 1))
        S = np.abs(np.fft.rfft(prof)) ** 2 / nb
        q = 2 * np.pi * np.arange(S.size) / (nb * dx)
        k = np.digitize(q[1:], edges) - 1
        ok = (k >= 0) & (k < nb_q)
        out_S.append(np.bincount(k[ok], S[1:][ok], minlength=nb_q))
        out_C.append(np.bincount(k[ok], minlength=nb_q))
        out_arc.append((b0 - a0) * step / l)
    if not out_S:
        return np.empty((0, nb_q)), np.empty((0, nb_q)), np.empty(0)
    return np.array(out_S), np.array(out_C), np.array(out_arc)


def q_edges(L):
    return np.logspace(np.log10(2 * np.pi / L / 1.5), np.log10(np.pi * 1.01), int(np.log10(1.5 * L / 2) * 8) + 2)


def zeta_band(qc, S, lo, hi):
    return -(fit_power(qc, S, lo, hi)[0] + 1) / 2


