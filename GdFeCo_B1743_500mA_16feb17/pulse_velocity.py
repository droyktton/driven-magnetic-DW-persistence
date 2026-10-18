"""Velocity vs coil current from the short pulse movies of 16 Feb 2017 (B1743, same field of view as Study 17).

Usage: python GdFeCo_B1743_500mA_16feb17/pulse_velocity.py [DATA_DIR] [--long I_mA V_nm_s LABEL ...]
DATA_DIR (default ~/MOKEdata_local/GdFeCo_B1743_2017_02_16) holds the folders B1743_Bg1iz<I>mA<tp>s_<n>/Pos0,
one image after each field pulse of duration tp (taken from the folder name; Interval_ms = tp + 8-10 s).
For every movie: log(frame k / frame 0), smoothed 2 px and thresholded at --thr (half the domain contrast,
0.15 in log units from the 900-1000 mA movies; background sd 0.006-0.011) inside the aperture, then opened, gives
the area swept after k pulses; divided by the number of rows the wall crosses, it is the mean displacement u_k.
v = Theil-Sen slope of u_k vs k tp over pulses 1..n (robust to single bad frames), n = last pulse before the swept
domain enters the left 15 % of the aperture. Pulse 1 is excluded from the slope: u_1 is often larger than the later
steps (first-pulse transient); it is reported separately. Points with u_n - u_1 < 3 um are flagged as unreliable.
--long adds steady-field runs (one frame every 2 s) to the figure as open symbols.
Outputs next to this script: pulse_velocity.txt and fig_velocity_current.png.
"""
import argparse
import glob
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy.ndimage import binary_opening, gaussian_filter
from scipy.stats import theilslopes

UM = 0.1175
ap = argparse.ArgumentParser()
ap.add_argument("data", nargs="?", default=os.path.expanduser("~/MOKEdata_local/GdFeCo_B1743_2017_02_16"))
ap.add_argument("--thr", type=float, default=0.075)
ap.add_argument("--long", nargs=3, action="append", default=[], metavar=("I_mA", "V_nm_s", "LABEL"))
a = ap.parse_args()
here = os.path.dirname(os.path.abspath(__file__))

rows_out = []
for d in sorted(glob.glob(os.path.join(a.data, "B1743_Bg1iz*mA*s_*"))):
    name = os.path.basename(d)
    I = int(name.split("iz")[1].split("mA")[0]); tp = float(name.split("mA")[1].split("s_")[0])
    fs = sorted(glob.glob(os.path.join(d, "Pos0", "*.tif")))
    if len(fs) > 100:                       # the long steady-field run, not a pulse movie
        continue
    im = [np.asarray(Image.open(f), float) for f in fs]
    inner = gaussian_filter(im[0], 5) > 0.6 * np.percentile(im[0], 90)
    cols = np.where(inner.any(0))[0]; strip = inner & (np.arange(im[0].shape[1]) < cols[0] + 0.15 * (cols[-1] - cols[0]))
    u, reach = [0.0], []
    for k in range(1, len(im)):
        x = gaussian_filter(np.log(np.maximum(im[k], 1) / np.maximum(im[0], 1)), 2)
        sw = binary_opening((x > a.thr) & inner, iterations=2)
        r = np.where(sw.sum(1) > 3)[0]
        u.append(sw[r].sum() / max(len(r), 1) * UM if len(r) else 0.0)
        reach.append(sw[strip].mean() > 0.01)
    u = np.array(u); t = np.arange(len(u)) * tp
    n = len(u) - 1 if not any(reach) else reach.index(True)   # last pulse before the edge
    ok = n >= 3 and np.median(u[n - 1:n + 1]) - u[1] >= 3
    v = theilslopes(u[1:n + 1], t[1:n + 1])[0] if n >= 2 else np.nan
    rows_out.append((I, tp, len(fs), u[-1], u[1], n - 1, v, ok))
    print(f"{name:26s} I={I:5d} mA tp={tp:4.0f} s u_k={np.round(u, 1)} | pulses 1..{n}: v={v*1e3:9.3g} nm/s"
          f"{'' if ok else '  (unreliable)'}")

R = np.array(rows_out, float)
np.savetxt(os.path.join(here, "pulse_velocity.txt"), R,
           header="I_mA tp_s frames u_total_um u_1_um steps_fitted v_um_s reliable", fmt="%.6g")
ok = R[:, 7] > 0
fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for m, lab, st in [(ok, "pulses, 16 Feb 2017", "o"), (~ok, "pulses, unreliable", "x")]:
    ax[0].semilogy(R[m, 0], R[m, 6] * 1e3, st, color="C0", label=lab)
    ax[1].plot(R[m, 0] ** -0.25, np.log(R[m, 6] * 1e3), st, color="C0", label=lab)
for i, (I, vv, lab) in enumerate(a.long):
    I, vv = float(I), float(vv)
    ax[0].semilogy(I, vv, "s", mfc="none", color=f"C{i+1}", ms=9, label=lab)
    ax[1].plot(I ** -0.25, np.log(vv), "s", mfc="none", color=f"C{i+1}", ms=9, label=lab)
if ok.sum() >= 3:
    p = np.polyfit(R[ok, 0] ** -0.25, np.log(R[ok, 6] * 1e3), 1)
    xx = np.linspace(R[:, 0].min() ** -0.25, R[:, 0].max() ** -0.25, 50)
    ax[1].plot(xx, np.polyval(p, xx), "k--", lw=1, label=f"creep fit: slope {p[0]:.0f}")
    print(f"ln v = {p[1]:.1f} {p[0]:+.1f} I^-1/4 (reliable pulse points, I in mA, v in nm/s)")
ax[0].set_xlabel("coil current I [mA]"); ax[0].set_ylabel("v [nm/s]"); ax[0].set_title("velocity vs current (B1743)")
ax[1].set_xlabel("I$^{-1/4}$ [mA$^{-1/4}$]"); ax[1].set_ylabel("ln v [nm/s]"); ax[1].set_title("creep plot (H ∝ I assumed)")
for x in ax:
    x.legend(fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(here, "fig_velocity_current.png"), dpi=90)
