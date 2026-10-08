#!/usr/bin/env python3
"""
wall_exponents.py -- test harness for local roughness-exponent estimators on
synthetic tortuous (possibly overhanging) walls with known zeta.

Wall model:  r(s) = r0(s) + u(s) * n(s)
  r0 : smooth base curve, sine-generated: tangent angle theta(s) = th0*sin(2 pi s/Lam + phase)
         flat      th0 = 0
         curved    th0 = 0.9 rad  (max tilt 52 deg, y(x) still single valued)
         overhang  th0 = 2.0 rad  (max tilt 115 deg, y(x) multivalued)
  u  : self-affine profile, Fourier amplitudes ~ k^-(zeta+1/2) for 1/xi <= k <= 1/2,
       i.e. S(k) ~ k^-(2 zeta+1); normalised so that the linearly detrended local
       width at l=64 equals EPS*64 (rms slope 2% there).

Estimators (they see only the ordered point set r_j, never u or h(x)):
  pca    rms perpendicular distance to the local total-least-squares line
  dfa1-3 rms residual of a local polynomial fit of order 1..3 in the local frame
  dog2-4 wavelet coefficient at the window centre, psi = derivative of Gaussian
         of order m (m vanishing moments), L1-normalised so that W ~ a^zeta
All of them work in a local frame (PCA axis of the window), after resampling the
contour at uniform arclength, and discard windows where t is not monotonic.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse
import csv
import pickle
import time
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view
from numpy.polynomial import hermite_e as He

# ------------------------------------------------------------------ settings
N_POINTS = 65536
GEOMS = {"flat": 0.0, "curved": 0.9, "overhang": 2.0}
LAMBDA = 8192.0                  # period of the sine-generated base curve
EPS, L0 = 0.02, 64.0             # amplitude normalisation (slope 2% at l=64)
SCALES_N = [12, 17, 24, 34, 48, 68, 96, 136, 192]   # window sizes (points) for pca / dfa
SCALES_A = [3, 4, 6, 8, 12, 17, 24, 34]             # wavelet scales (l ~ 4a)
FIT = slice(1, 7)                # fit n = 17..96 and a = 4..24 (l ~ 16..96); larger wavelet
                                 # windows (~15a points) stop being graph-like near the sharpest bends
METHODS = ["pca", "dfa1", "dfa2", "dfa3", "dog2", "dog3", "dog4"]
ZETAS = [0.3, 0.6, 0.9, 1.3, 1.7, 2.3]


# ---------------------------------------------------------------- generation
def spectral_profile(n, zeta, xi, rng):
    k = np.fft.rfftfreq(n)
    amp = np.zeros_like(k)
    sel = k >= 1.0 / xi
    amp[sel] = k[sel] ** (-(zeta + 0.5))
    c = (rng.standard_normal(k.size) + 1j * rng.standard_normal(k.size)) * amp
    c[0] = 0.0
    return np.fft.irfft(c, n=n)


def lin_width(u, n):
    """rms residual of an OLS line in windows of n+1 samples (used for calibration)."""
    w = sliding_window_view(u, n + 1)[:: max(n // 2, 1)]
    t = np.arange(n + 1) - n / 2.0
    wc = w - w.mean(axis=1, keepdims=True)
    b = (wc * t).sum(axis=1, keepdims=True) / (t * t).sum()
    return np.sqrt(np.mean((wc - b * t) ** 2))


def resample_arclength(X, Y, h=1.0):
    s = np.concatenate([[0.0], np.cumsum(np.hypot(np.diff(X), np.diff(Y)))])
    sn = np.arange(0.0, s[-1], h)
    return np.interp(sn, s, X), np.interp(sn, s, Y)


def make_wall(geom, zeta, seed, n=N_POINTS):
    rng = np.random.default_rng(seed)
    xi = 512.0 if zeta > 2 else 1024.0
    u = spectral_profile(n, zeta, xi, rng)
    u *= EPS * L0 / lin_width(u, int(L0))
    th0 = GEOMS[geom]
    phase = rng.uniform(0, 2 * np.pi)
    s = np.arange(n)
    th = th0 * np.sin(2 * np.pi * s / LAMBDA + phase)
    x, y = np.cumsum(np.cos(th)), np.cumsum(np.sin(th))
    X = x - u * np.sin(th)
    Y = y + u * np.cos(th)
    return resample_arclength(X, Y)


def diagnostics(X, Y):
    from scipy.spatial import cKDTree
    d = {"frac_back": float(np.mean(np.diff(X) < 0)), "npts": len(X)}
    pairs = cKDTree(np.c_[X, Y]).query_pairs(r=60.0, output_type="ndarray")
    if len(pairs):
        far = np.abs(pairs[:, 0] - pairs[:, 1]) > 1500
        pairs = pairs[far]
    if len(pairs):
        dist = np.hypot(X[pairs[:, 0]] - X[pairs[:, 1]], Y[pairs[:, 0]] - Y[pairs[:, 1]])
        d["min_far_dist"] = float(dist.min())
    else:
        d["min_far_dist"] = np.inf
    return d


# ------------------------------------------------------------ local frames
def local_frame_windows(X, Y, n, stride):
    """Windows of n+1 consecutive points -> (T, U) in each window's PCA frame.
    Keeps only windows where t is strictly monotonic. Returns T, U, kept fraction."""
    Xw = sliding_window_view(X, n + 1)[::stride]
    Yw = sliding_window_view(Y, n + 1)[::stride]
    xc = Xw - Xw.mean(axis=1, keepdims=True)
    yc = Yw - Yw.mean(axis=1, keepdims=True)
    cxx, cyy, cxy = (xc * xc).mean(1), (yc * yc).mean(1), (xc * yc).mean(1)
    phi = 0.5 * np.arctan2(2 * cxy, cxx - cyy)
    c, s = np.cos(phi)[:, None], np.sin(phi)[:, None]
    T = xc * c + yc * s
    U = -xc * s + yc * c
    flip = np.where(T[:, -1] - T[:, 0] < 0, -1.0, 1.0)[:, None]
    T, U = T * flip, U * flip
    valid = np.all(np.diff(T, axis=1) > 0, axis=1)
    return T[valid], U[valid], float(valid.mean())


def dfa_rms(T, U, order):
    t0, t1 = T[:, :1], T[:, -1:]
    tt = 2 * (T - t0) / (t1 - t0) - 1
    V = tt[:, :, None] ** np.arange(order + 1)
    Q, _ = np.linalg.qr(V)
    coef = np.einsum("mnp,mn->mp", Q, U)
    R = U - np.einsum("mnp,mp->mn", Q, coef)
    return float(np.sqrt(np.mean(R ** 2)))


def dog_kernel(a, m, H):
    x = np.arange(-H, H + 1) / a
    c = np.zeros(m + 1)
    c[m] = 1.0
    return He.hermeval(x, c) * np.exp(-x ** 2 / 2) / a      # L1 normalisation => W ~ a^zeta


def wavelet_rms(X, Y, a, ms):
    H = int(np.ceil(6 * a))
    # arclength of the window must exceed the needed t-extent 2H by enough margin that
    # coverage never depends on how wiggly the window is (otherwise windows are selected)
    n = int(np.ceil(2 * H * 2.2))
    T, U, vf = local_frame_windows(X, Y, n, max(1, int(a)))
    c = n // 2
    tc = T[:, c]
    ok = (T[:, 0] <= tc - H) & (T[:, -1] >= tc + H)
    kept = float(ok.mean()) * vf
    T, U, tc = T[ok], U[ok], tc[ok]
    M = T.shape[0]
    off = (np.arange(M) * (4.0 * n + 10))[:, None]
    Tf, Uf = (T + off).ravel(), U.ravel()
    Gf = (tc[:, None] + np.arange(-H, H + 1)[None, :] + off).ravel()
    idx = np.clip(np.searchsorted(Tf, Gf), 1, Tf.size - 1)
    x0, x1 = Tf[idx - 1], Tf[idx]
    w = (Gf - x0) / (x1 - x0)
    val = (Uf[idx - 1] * (1 - w) + Uf[idx] * w).reshape(M, -1)
    return {m: float(np.sqrt(np.mean((val @ dog_kernel(a, m, H)) ** 2))) for m in ms}, kept


def analyse(X, Y):
    out = {k: [] for k in ["ell", "pca", "dfa1", "dfa2", "dfa3", "dog2", "dog3", "dog4", "valid_n", "valid_a"]}
    for n in SCALES_N:
        T, U, vf = local_frame_windows(X, Y, n, max(1, n // 2))
        out["ell"].append(float(np.mean(T[:, -1] - T[:, 0])))
        out["pca"].append(float(np.sqrt(np.mean(U ** 2))))
        for p in (1, 2, 3):
            out[f"dfa{p}"].append(dfa_rms(T, U, p))
        out["valid_n"].append(vf)
    for a in SCALES_A:
        r, kept = wavelet_rms(X, Y, a, (2, 3, 4))
        for m in (2, 3, 4):
            out[f"dog{m}"].append(r[m])
        out["valid_a"].append(kept)
    return {k: np.array(v) for k, v in out.items()}


def fit_slope(x, F):
    x, F = np.asarray(x)[FIT], np.asarray(F)[FIT]
    return float(np.polyfit(np.log(x), np.log(F), 1)[0])


def zeta_hats(res):
    z = {m: fit_slope(res["ell"], res[m]) for m in ["pca", "dfa1", "dfa2", "dfa3"]}
    a = np.array(SCALES_A, float)
    z.update({m: fit_slope(a, res[m]) for m in ["dog2", "dog3", "dog4"]})
    return z


# ------------------------------------------------------------------- runner
def one_task(args):
    geom, zeta, seed = args
    X, Y = make_wall(geom, zeta, seed)
    res = analyse(X, Y)
    return {"geom": geom, "zeta": zeta, "seed": seed, "zhat": zeta_hats(res),
            "curves": res, "diag": diagnostics(X, Y)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nreal", type=int, default=16)
    ap.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    ap.add_argument("--out", default="results")
    ap.add_argument("--zetas", type=float, nargs="*", default=ZETAS)
    ap.add_argument("--geoms", nargs="*", default=list(GEOMS))
    a = ap.parse_args()

    tasks = [(g, z, 1000 * i + int(100 * z) + 7 * len(g))
             for g in a.geoms for z in a.zetas for i in range(a.nreal)]
    t0 = time.time()
    with ProcessPoolExecutor(a.workers) as ex:
        results = list(ex.map(one_task, tasks, chunksize=1))
    print(f"{len(tasks)} walls analysed in {time.time() - t0:.0f} s with {a.workers} workers")

    with open(a.out + ".pkl", "wb") as f:
        pickle.dump(results, f)

    rows = []
    for g in a.geoms:
        for z in a.zetas:
            sel = [r for r in results if r["geom"] == g and r["zeta"] == z]
            for m in METHODS:
                v = np.array([r["zhat"][m] for r in sel])
                rows.append([g, z, m, v.mean(), v.std(ddof=1), len(v)])
    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["geometry", "zeta_true", "method", "zeta_hat_mean", "zeta_hat_std", "n_walls"])
        w.writerows(rows)

    print(f"{'geometry':9s} {'zeta':>5s} " + " ".join(f"{m:>11s}" for m in METHODS))
    for g in a.geoms:
        for z in a.zetas:
            line = []
            for m in METHODS:
                r = [x for x in rows if x[0] == g and x[1] == z and x[2] == m][0]
                line.append(f"{r[3]:6.2f}±{r[4]:4.2f}")
            print(f"{g:9s} {z:5.2f} " + " ".join(line))


if __name__ == "__main__":
    main()
