"""Study 18, full frames (after lps_frames.py): KPZ term, overhangs and shape memory of the moving wall.

Usage: python lps_frames_post.py DIR [--n-boot 300]
DIR = <series>/frames (frames.npz, contours.npz of lps_frames.py); crossover.npz of the parent folder (step8 on Lucas's walls).
The images are differences to a reference frame: one edge of the bright stripe is the wall of the reference (it never moves),
the other is the wall. The moving edge is the one with the larger frame-to-frame motion; the protocol sweeps it back and forth
(4-11 um per pulse, reversing every 7-17 pulses).
  B. KPZ term as step10 (h(x, t) version): pairs of consecutive frames (one pulse), segments of L um, slope m by least squares,
     v_x = mean advance (sign of the pulse: median advance > 0), v_n = v_x cos(phi), v_n/<v_n> = a[1 + b dphi + kappa/2 dphi^2],
     Lambda = 1 + kappa; 16-84 % bootstrap over pulses. Also forward and backward pulses separately.
  Overhangs: writes DIR/contours.npz in the step11 format (moving wall only: stripe contour minus the points within 1.5 um of the
     reference edge and the columns where the stripe leaves the field of view; longest run; 0.5 px steps, px units, theta0 = 90)
     and DIR/meta.json, and copies crossover.npz; then run step11_overhangs.py DIR.  The original contours go to contours_stripe.npz.
  Memory: for every pair of frames (i, j), r(l) = <w2_l(h_i - h_j)> / (<w2_l(h_i)> + <w2_l(h_j)>), w2_l = variance about the
     linear fit in windows of l columns (common valid columns): 0 = same shape, 1 = independent shapes. Binned against the distance
     of the mean positions |<h_i> - <h_j>| and against l; pairs from different sweeps at the same position test the memory of
     the pinning landscape, pairs from the same sweep the decorrelation by motion.
Outputs: DIR/post.npz, DIR/fig_post.png, summary on stdout.
"""
import argparse
import json
import os
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from numpy.lib.stride_tricks import sliding_window_view as swv

ap = argparse.ArgumentParser()
ap.add_argument("dir")
ap.add_argument("--n-boot", type=int, default=300)
a = ap.parse_args()
D = a.dir
rng = np.random.default_rng(18)
F = np.load(os.path.join(D, "frames.npz"))
k = float(F["um_per_px"]); nums = F["frame"]
ht, hb = F["h_top"], F["h_bot"]
mov_t = np.nanmedian(np.abs(np.diff(ht, axis=0))); mov_b = np.nanmedian(np.abs(np.diff(hb, axis=0)))
h, href = (hb, ht) if mov_b > mov_t else (ht, hb)
nf, nx = h.shape
# drop the 5 um next to the ends of the valid range of every frame (end caps of the stripe)
trim = int(round(5 / k))
for i in range(nf):
    ok = np.nonzero(np.isfinite(h[i]))[0]
    if ok.size:
        h[i, :ok[0] + trim] = np.nan; h[i, ok[-1] - trim + 1:] = np.nan
pos = np.nanmedian(h, 1)

# ------------------------------------------------------------------ B. KPZ term from consecutive pulses
pairs = [i for i in range(nf - 1) if nums[i + 1] - nums[i] == 1]
sgn = {i: np.sign(np.nanmedian(h[i + 1] - h[i])) for i in pairs}
Ls = [2.0, 4.0, 8.0, 16.0]


def seg_b(L_um, sel):
    L = int(round(L_um / k)); y = (np.arange(L) - (L - 1) / 2) * k
    P, V, B = [], [], []
    for i in sel:
        w0 = swv(h[i] * sgn[i], L)[::max(L // 2, 1)]; w1 = swv(h[i + 1] * sgn[i], L)[::max(L // 2, 1)]
        ok = np.isfinite(w0).all(1) & np.isfinite(w1).all(1)
        m = ((w0 - w0.mean(1, keepdims=True)) * y).sum(1) / (y * y).sum()
        phi = -np.arctan(m); vx = (w1 - w0).mean(1)
        P.append(phi[ok]); V.append((vx * np.cos(phi))[ok]); B.append(np.full(ok.sum(), i))
    return np.concatenate(P), np.concatenate(V), np.concatenate(B)


def fit(dphi, v):
    y = v / v.mean(); A = np.vstack([np.ones_like(dphi), dphi, dphi * dphi]).T
    c = np.linalg.lstsq(A, y, rcond=None)[0]
    return 2 * c[2] / c[0]


def lam(L_um, sel):
    phi, v, blk = seg_b(L_um, sel)
    if phi.size < 200:
        return None
    dphi = phi - np.angle(np.mean(np.exp(1j * phi))); keep = np.abs(dphi) < np.radians(60)
    dphi, v, blk = dphi[keep], v[keep], blk[keep]
    k0 = fit(dphi, v)
    ub = np.unique(blk); idx = {b: np.nonzero(blk == b)[0] for b in ub}
    bs = [fit(dphi[j], v[j]) for j in (np.concatenate([idx[b] for b in rng.choice(ub, ub.size)]) for _ in range(a.n_boot))]
    m2 = np.mean(np.tan(dphi) ** 2)
    edges = np.radians(np.arange(-60, 61, 10)); ib = np.digitize(dphi, edges) - 1
    curve = np.array([v[ib == q].mean() / v.mean() if (ib == q).sum() > 20 else np.nan for q in range(len(edges) - 1)])
    return dict(Lambda=1 + k0, ci=1 + np.percentile(bs, [16, 84]), n=int(v.size), vmean=float(v.mean()),
                R=abs(1 + k0) * m2 / 2, curve=curve)


B = {}
for name, sel in (("all", pairs), ("forward", [i for i in pairs if sgn[i] > 0]), ("backward", [i for i in pairs if sgn[i] < 0])):
    for L in Ls:
        r = lam(L, sel)
        if r:
            B[(name, L)] = r
adv = np.array([np.nanmedian((h[i + 1] - h[i]) * sgn[i]) for i in pairs])

# ------------------------------------------------------------------ overhangs: step11 input
C = np.load(os.path.join(D, "contours.npz"))
if "x" not in C.files:   # first run: keep the stripe contours
    shutil.copy(os.path.join(D, "contours.npz"), os.path.join(D, "contours_stripe.npz"))
C = np.load(os.path.join(D, "contours_stripe.npz"))
X0 = np.arange(nx) * k
ys, xs, off, fr = [], [], [0], []
for i, n in enumerate(nums):
    c = C[f"c{n:05d}"]
    if c.shape[0] < 50:
        continue
    cx, cy = c[:, 0], c[:, 1]
    ok = np.isfinite(h[i]); xr = X0[ok]
    if xr.size < 50:
        continue
    keep = (cx > xr.min()) & (cx < xr.max())
    yref = np.interp(cx, X0[np.isfinite(href[i])], href[i][np.isfinite(href[i])])
    keep &= np.abs(cy - yref) > 1.5
    # longest run of kept points (the contour is closed: roll so that it starts at a dropped point)
    if keep.all() or not keep.any():
        continue
    s = np.argmin(keep); keep, cx, cy = np.roll(keep, -s), np.roll(cx, -s), np.roll(cy, -s)
    d = np.diff(np.r_[0, keep.astype(int), 0]); st, en = np.nonzero(d == 1)[0], np.nonzero(d == -1)[0]
    j = np.argmax(en - st); px, py = cx[st[j]:en[j]] / k, cy[st[j]:en[j]] / k
    sarc = np.r_[0, np.cumsum(np.hypot(np.diff(px), np.diff(py)))]
    if sarc[-1] < 100:
        continue
    g = np.arange(0, sarc[-1], 0.5)
    ys.append(np.interp(g, sarc, py)); xs.append(np.interp(g, sarc, px)); off.append(off[-1] + g.size); fr.append(n)
np.savez_compressed(os.path.join(D, "contours.npz"), y=np.concatenate(ys), x=np.concatenate(xs), offsets=np.array(off),
                    frames=np.array(fr), theta0=90.0, step=0.5, smoothing_px=float(F["sigma"]))
json.dump(dict(um_per_px=k, source="lps_frames.py", moving_wall_only=True), open(os.path.join(D, "meta.json"), "w"), indent=1)
par = os.path.dirname(os.path.normpath(D))
if os.path.exists(os.path.join(par, "crossover.npz")):
    shutil.copy(os.path.join(par, "crossover.npz"), os.path.join(D, "crossover.npz"))

# ------------------------------------------------------------------ shape memory
ells = np.array([1, 2, 4, 8, 16, 32]) / k


def w2l(p, L):
    L = int(L); w = swv(p, L)[::max(L // 2, 1)]; ok = np.isfinite(w).all(1)
    if ok.sum() < 3:
        return np.nan
    w = w[ok]; y = np.arange(L) - (L - 1) / 2
    m = ((w - w.mean(1, keepdims=True)) * y).sum(1) / (y * y).sum()
    r = w - w.mean(1, keepdims=True) - m[:, None] * y
    return float((r ** 2).mean())


# sweep index: a new sweep starts where the pulse sign changes
sw = np.zeros(nf, int); last = None
for i in range(1, nf):
    s_ = sgn.get(i - 1, None)
    sw[i] = sw[i - 1] + (1 if (s_ is not None and last is not None and s_ != last) or nums[i] - nums[i - 1] > 1 else 0)
    last = s_ if s_ is not None else last
rows = []
for i in range(nf):
    for j in range(i + 1, nf):
        com = np.isfinite(h[i]) & np.isfinite(h[j])
        if com.sum() < 2 * ells.max():
            continue
        hi, hj = np.where(com, h[i], np.nan), np.where(com, h[j], np.nan)
        r = [w2l(hi - hj, L) / (w2l(hi, L) + w2l(hj, L)) for L in ells]
        rows.append([abs(pos[i] - pos[j]), sw[i] != sw[j], *r])
rows = np.array(rows, float)
dbins = np.array([0, 1, 2, 4, 8, 16, 32, 64])
mem = np.full((2, len(dbins) - 1, len(ells)), np.nan); cnt = np.zeros((2, len(dbins) - 1), int)
for q, other in enumerate((0, 1)):
    for b in range(len(dbins) - 1):
        s_ = (rows[:, 1] == other) & (rows[:, 0] >= dbins[b]) & (rows[:, 0] < dbins[b + 1])
        cnt[q, b] = s_.sum()
        if s_.sum() >= 3:
            mem[q, b] = np.nanmedian(rows[s_, 2:], 0)

# ------------------------------------------------------------------ output
name = os.path.basename(os.path.dirname(os.path.normpath(D)))
print(f"=== {name}: {nf} frames, {len(pairs)} single pulses, median advance {np.median(np.abs(adv)):.1f} um/pulse, "
      f"{sw.max() + 1} sweeps")
for (nm, L), r in sorted(B.items()):
    print(f"   B {nm:8s} L = {L:4.0f} um: Lambda = {r['Lambda']:+.2f} ({r['ci'][0]:+.2f} to {r['ci'][1]:+.2f}), R = {r['R']:.3f}, n = {r['n']}")
print("   memory r(l) (0 same shape, 1 independent); rows |d<h>| [um], columns l = " + " ".join(f"{e * k:.0f}" for e in ells) + " um")
for q, lab in enumerate(("same sweep", "other sweep")):
    for b in range(len(dbins) - 1):
        if np.isfinite(mem[q, b]).any():
            print(f"     {lab:11s} {dbins[b]:3d}-{dbins[b + 1]:3d} um (n={cnt[q, b]:4d}): " + " ".join(f"{v:5.2f}" for v in mem[q, b]))
np.savez(os.path.join(D, "post.npz"), Ls=Ls, **{f"B_{nm}_L{L:g}": np.r_[r["Lambda"], r["ci"], r["R"], r["n"]] for (nm, L), r in B.items()},
         **{f"curve_{nm}_L{L:g}": r["curve"] for (nm, L), r in B.items()},
         adv=adv, ells=ells * k, dbins=dbins, mem=mem, mem_count=cnt, pairs_rows=rows, sweep=sw, pos=pos)

fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
ax[0].plot(nums, pos, ".-"); ax[0].set_xlabel("frame (pulse)"); ax[0].set_ylabel("mean wall position [um]")
ax[0].set_title(f"{name}: wall sweeps (reference edge at {np.nanmedian(href):.0f} um)")
cen = np.arange(-55, 56, 10)
for L, col in zip(Ls, ("C0", "C1", "C2", "C3")):
    if ("all", L) in B:
        r = B[("all", L)]; ax[1].plot(cen, r["curve"], "o-", color=col, label=f"L = {L:g} um: Λ = {r['Lambda']:+.2f}")
ax[1].set_xlabel("dphi [deg]"); ax[1].set_ylabel("v_n / <v_n>"); ax[1].legend(fontsize=8)
ax[1].set_title("B. normal velocity vs orientation (one pulse)")
dc = np.sqrt(dbins[:-1] * dbins[1:]); dc[0] = 0.5
for q, (lab, mk) in enumerate((("same sweep", "o-"), ("other sweep", "s--"))):
    for e, col in zip((1, 3, 5), ("C0", "C2", "C3")):
        ax[2].semilogx(dc, mem[q, :, e], mk, color=col, label=f"{lab}, l = {ells[e] * k:.0f} um")
ax[2].axhline(1, color="0.5", lw=0.8); ax[2].set_xlabel("|<h_i> - <h_j>| [um]"); ax[2].set_ylabel("r(l): 0 same shape, 1 independent")
ax[2].set_title("shape memory between frames"); ax[2].legend(fontsize=7)
fig.tight_layout(); fig.savefig(os.path.join(D, "fig_post.png"), dpi=100)
