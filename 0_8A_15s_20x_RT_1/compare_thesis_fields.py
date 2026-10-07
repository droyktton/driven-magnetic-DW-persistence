"""Compare Study 2 avalanches with the thesis data at 46.1 Oe (this measurement) and 24.2 Oe.

The thesis (M. Grassi) analysed 0_8A_15s_20x_RT_1 (46.1 Oe, room temperature) and
0_42A_20s_20x_52C_1 (24.2 Oe, ~52 °C) with its own method; only the processed results of the
second measurement are available. This script puts them next to ours:
  (a) mean avalanche area vs measurement window, (b) number of avalanches vs window,
  (c) area distributions at one frame rescaled by the mean area, as the fraction of
      avalanches per logarithmic bin (the thesis files are normalised densities p(S), so the
      fraction per log bin is proportional to S p(S)).
Thesis windows for 0_8A are labelled 15, 30, ... s, i.e. 15 s per frame; the real interval
was 20 s, so they are converted to frames with 15 s/frame. For 0_42A the labels are
20, 40, ... s (20 s per frame).

Usage: python 0_8A_15s_20x_RT_1/compare_thesis_fields.py THESIS_DATA_DIR
       (THESIS_DATA_DIR = …/ResumenMaestría_Mati/Data)
Output: 0_8A_15s_20x_RT_1/fig_thesis_fields.png
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import center_of_mass, label

D = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
e = np.load(os.path.join(here, "tiff_extra.npz"))
k_um = 0.1175
X0, X1, Y0, Y1 = 430, 790, 528, 790      # thesis ROI (camera px)
WINDOWS = (1, 2, 3, 4, 5, 6, 8, 10, 12, 20, 30, 40, 60)


def ours():
    ta, y0, defects = e["arrival"], int(e["rows"][0]), e["defects"]
    M = 60
    sub = np.s_[Y0 - y0 - M:Y1 - y0 + M, X0 - M:X1 + M]
    t = ta[sub]
    valid = np.isfinite(t) & (t >= 1) & ~defects[sub]
    K = np.where(valid, np.round(t), -1).astype(int)
    rows, S1 = [], None
    for m in WINDOWS:
        win = np.where(K >= 1, (K - 1) // m, -1)
        S = []
        for j in range(win.max() + 1):
            mk = win == j
            if not mk.any():
                continue
            lab, n = label(mk, structure=np.ones((3, 3), bool))
            sz = np.bincount(lab.ravel())[1:]
            keep = np.nonzero(sz >= 20)[0] + 1
            if keep.size == 0:
                continue
            for i, (cy, cx) in zip(keep, center_of_mass(mk, lab, keep)):
                if M <= cx < M + (X1 - X0) and M <= cy < M + (Y1 - Y0):
                    S.append(sz[i - 1])
        S = np.array(S) * k_um ** 2
        rows.append((m, S.size, S.mean()))
        if m == 1:
            S1 = S
    return np.array(rows), S1


def logbin(v, vmin, nb=9):
    edges = np.logspace(np.log10(vmin), np.log10(v.max() * 1.0001), nb + 1)
    n, _ = np.histogram(v, edges)
    c = np.sqrt(edges[1:] * edges[:-1])
    ok = n > 0
    return c[ok], (n / np.diff(edges) / v.size)[ok]


o, S1 = ours()
m08 = np.loadtxt(os.path.join(D, "S_medio", "Media_0_8A.txt"))
n08 = np.loadtxt(os.path.join(D, "S_medio", "N_0_8A.txt"))
m42 = np.loadtxt(os.path.join(D, "S_medio", "Media_0_42A.txt"))
h08 = np.loadtxt(os.path.join(D, "Histogramas_log", "histograma_0_8.txt"))
h42 = np.loadtxt(os.path.join(D, "Histogramas_log", "histograma_0_42.txt"))

fig, ax = plt.subplots(1, 3, figsize=(17, 4.8))
ax[0].loglog(m08[:, 0] / 15, m08[:, 1], "s--", mfc="none", color="C3", label="thesis, 46.1 Oe RT")
ax[0].loglog(o[:, 0], o[:, 2], "o-", color="C0", label="ours, 46.1 Oe RT (thesis ROI)")
ax[0].loglog(m42[:, 0] / 20, m42[:, 1], "^--", mfc="none", color="C2", label="thesis, 24.2 Oe ~52 °C")
ax[0].set_xlabel(r"window $\tau_m$ [frames]"); ax[0].set_ylabel(r"mean area $\bar S$ [µm²]")
ax[0].legend(fontsize=8); ax[0].set_title("Mean avalanche area", fontsize=10)
ax[1].loglog(n08[:, 0] / 15, n08[:, 1], "s--", mfc="none", color="C3", label="thesis, 46.1 Oe RT")
ax[1].loglog(o[:, 0], o[:, 1], "o-", color="C0", label="ours, 46.1 Oe RT (thesis ROI)")
ax[1].set_xlabel(r"window $\tau_m$ [frames]"); ax[1].set_ylabel("number of avalanches N")
ax[1].legend(fontsize=8); ax[1].set_title("Number of avalanches (same ROI)", fontsize=10)
x, d = logbin(S1, 20 * k_um ** 2)
s08, s42, sO = m08[0, 1], m42[0, 1], S1.mean()
ax[2].loglog(h08[:, 0] / s08, h08[:, 0] * h08[:, 1] / np.sum(h08[:, 0] * h08[:, 1]), "s--", mfc="none", color="C3",
             label=rf"thesis 46.1 Oe ($\bar S$={s08:.2f} µm²)")
ax[2].loglog(h42[:, 0] / s42, h42[:, 0] * h42[:, 1] / np.sum(h42[:, 0] * h42[:, 1]), "^--", mfc="none", color="C2",
             label=rf"thesis 24.2 Oe ($\bar S$={s42:.2f} µm²)")
ax[2].loglog(x / sO, d * x / np.sum(d * x), "o-", color="C0", label=rf"ours 46.1 Oe ($\bar S$={sO:.2f} µm²)")
ax[2].set_xlabel(r"$S/\bar S$"); ax[2].set_ylabel("fraction per log bin")
ax[2].legend(fontsize=8); ax[2].set_title("Area distributions at 1 frame, rescaled", fontsize=10)
fig.tight_layout(); fig.savefig(os.path.join(here, "fig_thesis_fields.png"), dpi=130)

print("frames  N_thesis  N_ours   S_thesis  S_ours  difference")
for m, n, s in o:
    t = np.interp(m, m08[:, 0] / 15, m08[:, 1])
    nt = np.interp(m, n08[:, 0] / 15, n08[:, 1])
    print(f"{int(m):5d}  {nt:7.0f}  {int(n):6d}   {t:7.2f}  {s:6.2f}  {t - s:+6.2f}")
print("thesis S_mean at 24.2 Oe vs 46.1 Oe (same method), by frames:")
for fr in (1, 2, 4, 8, 16):
    a, b = np.interp(fr, m42[:, 0] / 20, m42[:, 1]), np.interp(fr, m08[:, 0] / 15, m08[:, 1])
    print(f"  {fr:3d} frames: {a:.2f} vs {b:.2f} µm²  (ratio {a / b:.2f})")
