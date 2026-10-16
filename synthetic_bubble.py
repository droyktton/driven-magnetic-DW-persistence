"""Synthetic arrival map of a growing rough bubble with a known roughness exponent (calibration of the estimators).

Usage: python synthetic_bubble.py ZETA SEED --out FILE[.csv|.npz] [--w50 1.0] [--size 1040] [--r0 106] [--r1 400]
ZETA: "Z" (pure power law) or "Z1:Z2:LC" (zeta Z1 below LC px of arc, Z2 above; continuous broken spectrum).
Wall r(theta) = R + (R/Rcal)^zeta u(theta): u has modes n = 2..nmax with S(q) ~ q^-(1+2 zeta) (random amplitudes and
phases) and is fixed in theta, rescaled during the growth so that its amplitude is stationary in arc length and the
domains are nested (each seed is one realisation of the wall shape). Amplitude set so that the detrended local width
over 50 px of arc is --w50 px at Rcal = 0.7 r1. A pixel switches when R reaches r - (R/Rcal)^zeta u(theta).
Output: .csv = event map as bubble_event_map.py reads it (index = R level x 64 angular sectors, -1 = initial bubble
and outside); .npz = arrival map as contour_analysis.py reads it (ta = R level, 1 px per frame; 0 = initial bubble; -1 never).
"""
import argparse

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("zeta")
ap.add_argument("seed", type=int)
ap.add_argument("--out", required=True)
ap.add_argument("--w50", type=float, default=1.0, help="local width over 50 px of arc, px")
ap.add_argument("--size", type=int, default=1040)
ap.add_argument("--r0", type=float, default=106, help="initial radius, px")
ap.add_argument("--r1", type=float, default=400, help="final radius, px")
args = ap.parse_args()

parts = [float(v) for v in args.zeta.split(":")]
zeta = parts[1] if len(parts) == 3 else parts[0]          # large-scale exponent, used for the growth rescaling
N, c = args.size, (args.size - 1) / 2
R0, R1 = args.r0, args.r1
Rcal = 0.7 * R1
nmax = int(np.pi * R1)
rng = np.random.default_rng(args.seed)
n = np.arange(2, nmax + 1)
if len(parts) == 3:
    z1, z2, lc = parts
    nc = 2 * np.pi * Rcal / lc
    amp = np.where(n < nc, (n / nc) ** (-(1 + 2 * z2) / 2), (n / nc) ** (-(1 + 2 * z1) / 2))
else:
    amp = n ** (-(1 + 2 * zeta) / 2)
ca, sa = amp * rng.standard_normal(n.size), amp * rng.standard_normal(n.size)

M = 2 ** 15                                               # dense profile in theta
th = 2 * np.pi * np.arange(M) / M
uc = np.zeros(M)
for k0 in range(0, M, 4096):
    t = th[k0:k0 + 4096]
    uc[k0:k0 + 4096] = np.cos(np.outer(t, n)) @ ca + np.sin(np.outer(t, n)) @ sa
ds = 2 * np.pi * Rcal / M
m = int(round(50 / ds))
x = np.arange(m) * ds
v = [np.var(uc[i:i + m] - np.polyval(np.polyfit(x, uc[i:i + m], 1), x)) for i in range(0, M - m, m // 2)]
A = args.w50 / np.sqrt(np.mean(v))

yy, xx = np.mgrid[0:N, 0:N]
r = np.hypot(yy - c, xx - c)
t = np.arctan2(yy - c, xx - c)
up = A * np.interp(t.ravel() % (2 * np.pi), th, uc, period=2 * np.pi).reshape(N, N)
Rneed = r - up
for _ in range(30):
    Rneed = r - (np.clip(Rneed, 1, None) / Rcal) ** zeta * up
inside = (Rneed >= R0) & (Rneed < R1)
lev = np.floor(Rneed - R0).astype(int)
if args.out.endswith(".csv"):
    sec = ((t + np.pi) / (2 * np.pi) * 64).astype(int).clip(0, 63)
    np.savetxt(args.out, np.where(inside, lev * 64 + sec, -1), fmt="%d", delimiter=",")
else:
    ta = np.where(inside, lev + 1, -1)
    ta[Rneed < R0] = 0                                     # initial bubble: switched before the first frame
    np.savez_compressed(args.out, ta=ta.astype(np.int16))
print(f"zeta {args.zeta} seed {args.seed}: amplitude {A:.3g}, rms u = {A * uc.std():.2f} px, max |u| = {A * np.abs(uc).max():.1f} px")
