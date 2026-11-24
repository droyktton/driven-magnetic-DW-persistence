"""Forward model of the TbCo X20 wall measurement + S(q) estimator on arc windows (data and synthetic, same code).

Synthetic wall h(x) (NX px, Gaussian, S(q) ~ q^-(1+2 z1) above q_x = 2 pi / Lx, q^-(1+2 z2) below) -> domain image
(pixel area fraction) -> Gaussian PSF sigma_psf -> white noise sd NOISE (contrast units) -> Gaussian sigma 2 -> contour
at 0.5 (subpixel), resampled at 0.5 px -> same local width (wall_tools.overhang_segments) and same S(q) as the data.
Amplitude: true wall scaled so that the MEASURED w2 at E ~ 10 um matches the data (creep group).
"""
import json, sys
import numpy as np
from scipy import ndimage as ndi
from skimage.measure import find_contours
sys.path.insert(0, '/home/koltona/Codigos/driven-magnetic-DW-persistence')
from wall_tools import overhang_segments, fit_power

R = '/home/koltona/Codigos/driven-magnetic-DW-persistence/TbCo_Paris_Nov2023/AV6517_1_0_#2_0.83nm/'
S_UM = 0.136
STEP = 0.5
NOISE = 0.18
J = json.load(open(R + 'pulses.json'))
LAMS = np.array(J['lams_px'])
rng = np.random.default_rng(11)
SPEC_L = 256                      # arc window for S(q) [px] (35 um)
QE = np.logspace(np.log10(2 * np.pi / SPEC_L), np.log10(np.pi), 22)


def spec_arc(y, x, l=SPEC_L):
    """S(q) of windows of arc length l along a contour piece, each rotated to its own principal axis."""
    n = int(round(l / STEP)); out = []
    for a0 in range(0, y.size - n, n // 2):
        py, px = y[a0:a0 + n], x[a0:a0 + n]
        my, mx = py.mean(), px.mean()
        cyy, cxx, cxy = ((py - my) ** 2).mean(), ((px - mx) ** 2).mean(), ((py - my) * (px - mx)).mean()
        al = 0.5 * np.arctan2(2 * cxy, cyy - cxx)
        ax_ = (py - my) * np.cos(al) + (px - mx) * np.sin(al)
        nn = -(py - my) * np.sin(al) + (px - mx) * np.cos(al)
        ib = np.floor(ax_ - ax_.min()).astype(int); nb = ib.max() + 1
        if nb < 64:
            continue
        cnt = np.bincount(ib, minlength=nb); prof = np.bincount(ib, nn, minlength=nb) / np.maximum(cnt, 1)
        g = cnt > 0
        prof = np.interp(np.arange(nb), np.nonzero(g)[0], prof[g])
        e0, e1 = prof[:3].mean(), prof[-3:].mean()
        prof = prof - (e0 + (e1 - e0) * (np.arange(nb) - 1) / max(nb - 3, 1))
        S = np.abs(np.fft.rfft(prof)) ** 2 / nb
        q = 2 * np.pi * np.arange(S.size) / nb
        k = np.digitize(q[1:], QE) - 1; ok = (k >= 0) & (k < QE.size - 1)
        out.append((np.bincount(k[ok], S[1:][ok], minlength=QE.size - 1), np.bincount(k[ok], minlength=QE.size - 1)))
    return out


def estimators(pieces):
    w2s = np.zeros(LAMS.size); w2n = np.zeros(LAMS.size); Ss = np.zeros(QE.size - 1); Sn = np.zeros(QE.size - 1)
    for y, x in pieces:
        _, _, _, ww = overhang_segments(y, x, LAMS, STEP, [3.0])
        for j in range(LAMS.size):
            w2s[j] += np.sum(ww[j]); w2n[j] += len(ww[j])
        for s_, c_ in spec_arc(y, x):
            Ss += s_; Sn += c_
    return w2s, w2n, Ss, Sn


def wall(nx, z1, z2, lx):
    m = 4 * nx; q = 2 * np.pi * np.fft.rfftfreq(m); S = np.zeros_like(q)
    qx = 2 * np.pi / lx if lx else 0.0
    hi = q >= qx; S[1:][hi[1:]] = q[1:][hi[1:]] ** (-(1 + 2 * z1))
    lo = (~hi) & (q > 0)
    if lo.any():
        S[lo] = qx ** (-(1 + 2 * z1)) * (q[lo] / qx) ** (-(1 + 2 * z2))
    a = np.sqrt(S) * (rng.normal(size=q.size) + 1j * rng.normal(size=q.size))
    h = np.fft.irfft(a, m)[:nx]
    return h - h.mean()


def measure(h, psf):
    nx = h.size; ny = int(np.ptp(h) + 120)
    xw = h - h.min() + 60
    yy = np.arange(ny)[None, :]
    img = np.clip(xw[:, None] - yy + 0.5, 0, 1)                    # rows = along the wall, columns = across
    img = ndi.gaussian_filter(img, psf) + rng.normal(0, NOISE, img.shape)
    s = ndi.gaussian_filter(img, 2.0)
    pieces = []
    for c in find_contours(s, 0.5):
        keep = (c[:, 0] > 20) & (c[:, 0] < nx - 20)
        c = c[keep]
        if c.shape[0] < 200:
            continue
        seg = np.hypot(*np.diff(c, axis=0).T); arc = np.r_[0, np.cumsum(seg)]
        if arc[-1] < 300:
            continue
        t = np.arange(0, arc[-1], STEP)
        pieces.append((np.interp(t, arc, c[:, 0]), np.interp(t, arc, c[:, 1])))
    return pieces


def true_pieces(h):
    t = np.arange(20, h.size - 20, STEP)
    return [(t, np.interp(t, np.arange(h.size), h))]


def synth(z1, z2, lx, psf, amp, nrep=12, nx=1200):
    acc_m = [np.zeros(LAMS.size), np.zeros(LAMS.size), np.zeros(QE.size - 1), np.zeros(QE.size - 1)]
    acc_t = [np.zeros(LAMS.size), np.zeros(LAMS.size), np.zeros(QE.size - 1), np.zeros(QE.size - 1)]
    for _ in range(nrep):
        h = amp * wall(nx, z1, z2, lx)
        for acc, p in ((acc_m, measure(h, psf)), (acc_t, true_pieces(h))):
            e = estimators(p)
            for i in range(4):
                acc[i] += e[i]
    f = lambda A: (A[0] / np.maximum(A[1], 1), A[2] / np.maximum(A[3], 1))
    return f(acc_m), f(acc_t)


def data_curves(sel):
    W = np.load(R + 'walls.npz')
    pcs = []
    for k in sel:
        for f in W.files:
            if f.startswith(f'f{k - 1}_'):
                c = W[f]; pcs.append((c[:, 0], c[:, 1]))
    e = estimators(pcs)
    return e[0] / np.maximum(e[1], 1), e[2] / np.maximum(e[3], 1), e[1]


def slopes(lams, w2, lo, hi):
    ok = (lams >= lo) & (lams <= hi) & (w2 > 0)
    return fit_power(lams[ok], w2[ok], lo, hi)[0] / 2


def zS(S, lo_um, hi_um):
    qc = np.sqrt(QE[1:] * QE[:-1]); lam = 2 * np.pi / qc * S_UM
    ok = (lam >= lo_um) & (lam <= hi_um) & (S > 0)
    return -(np.polyfit(np.log(qc[ok]), np.log(S[ok]), 1)[0] + 1) / 2


if __name__ == '__main__':
    groups = {'creep': [p['file'] for p in J['pulses'] if p['H'] < 82],
              'nearHd': [p['file'] for p in J['pulses'] if 82 <= p['H'] < 104]}
    res = dict(lams=LAMS.tolist(), qe=QE.tolist(), data={}, models=[])
    for g, sel in groups.items():
        w2, S, n = data_curves(sel)
        res['data'][g] = dict(w2=w2.tolist(), S=S.tolist())
        print(f"DATA {g}: w2 zeta 2-5/5-10/2-20 um = {slopes(LAMS*S_UM, w2, 2, 5):.2f}/{slopes(LAMS*S_UM, w2, 5, 10):.2f}/"
              f"{slopes(LAMS*S_UM, w2, 2, 20):.2f};  S(q) zeta lambda 2-6/6-30 um = {zS(S, 2, 6):.2f}/{zS(S, 6, 30):.2f}", flush=True)
    wref = np.array(res['data']['creep']['w2']); jref = int(np.argmin(abs(LAMS * S_UM - 10)))
    models = [(z, z, 0) for z in (0.6, 0.8, 1.0, 1.25, 1.5)] + [(1.25, z2, L) for z2 in (0.6, 0.8) for L in (3, 6, 12)] + [(1.0, 0.6, 6)]
    for (z1, z2, Lum) in models:
        for psf in (2.5, 3.3):
            amp = 1.0
            for it in range(3):
                (wm, Sm), (wt, St) = synth(z1, z2, Lum / S_UM if Lum else 0, psf, amp, nrep=4 if it < 2 else 16)
                if it < 2:
                    amp *= np.sqrt(wref[jref] / wm[jref])
            rec = dict(z1=z1, z2=z2, Lx_um=Lum, psf=psf, amp=amp, w2=wm.tolist(), w2_true=wt.tolist(), S=Sm.tolist(), S_true=St.tolist())
            res['models'].append(rec)
            print(f"z1 {z1:.2f} z2 {z2:.2f} Lx {Lum:>2} um psf {psf}: measured w2 zeta 2-5/5-10/2-20 = "
                  f"{slopes(LAMS*S_UM, wm, 2, 5):.2f}/{slopes(LAMS*S_UM, wm, 5, 10):.2f}/{slopes(LAMS*S_UM, wm, 2, 20):.2f} "
                  f"(true {slopes(LAMS*S_UM, wt, 2, 5):.2f}/{slopes(LAMS*S_UM, wt, 5, 10):.2f}/{slopes(LAMS*S_UM, wt, 2, 20):.2f}); "
                  f"S(q) 2-6/6-30 = {zS(Sm, 2, 6):.2f}/{zS(Sm, 6, 30):.2f} (true {zS(St, 2, 6):.2f}/{zS(St, 6, 30):.2f})", flush=True)
            json.dump(res, open('fwd_results.json', 'w'))
