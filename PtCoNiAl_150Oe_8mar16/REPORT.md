# Study 13: `PtCoNiAl_150Oe_8mar16` (PtCoNiAl, 150 Oe, nucleating bubbles, Δt = 3 ms)

First report, 2026-10-16. Tools: `arrival_map.py`, `contour_analysis.py` (see the [README](../README.md)). Independent of the other studies.

## Data
- **Source:** group movie `PtCoNiAl_8Marzo2016_BCI450mA_H150Oe_5X_650umX486um_Dt3mSeg_5fps_0_015vecesRealTime BIS.avi` (not in the repository).
  - 35 frames, 1392×1040 px, compressed video.
  - From the file name: PtCoNiAl, 450 mA = **150 Oe**, 5X, 650 × 486 µm, so **0.467 µm/px**; **Δt = 3 ms**.
- **Content:** about ten bubbles nucleate (some in clusters at the bottom) and merge within ~15 frames. Switching makes the pixels brighter.

```bash
python arrival_map.py ".../PtCoNiAl_8Marzo2016_BCI450mA_H150Oe_5X_650umX486um_Dt3mSeg_5fps_0_015vecesRealTime BIS.avi" \
       --outdir PtCoNiAl_150Oe_8mar16 --qc 3,5,7,9,11,13,16,20,26
python contour_analysis.py PtCoNiAl_150Oe_8mar16/arrival.npz PtCoNiAl_150Oe_8mar16 --frames 4 12 --um-per-px 0.467 \
       --centre 314 875 --roi-erode 30
```

## Arrival map and checks
- Otsu threshold 13.0 (median amplitude 24.4 switched, 0 not); 91 % of the image switched.
- **Segmentation check:** `fig_qc_overlay.png` was approved by the user. The nucleation points are "never switched" (already switched in frame 0); noise on the left rim of the aperture (ROI eroded by 30 px).
- **Bubble analysed:** the isolated one at (y, x) = (314, 875) px, frames 4–12, before it merges with the bubble above. Mean radius 57 → 175 px (27 → 82 µm), 780 px of wall per frame.
- **Velocity:** 14.8 px per frame = 6.9 µm / 3 ms = **2.3 mm/s** (flow or depinning regime rather than creep).

## Results
| quantity | value |
|---|---|
| w at E = 10 / 20 / 50 px | 0.38 / 0.66 / 1.07 px |
| local width ζ, raw | 0.68 (8–80 px), 0.51 (20–80 px); windows only up to 53 px (the smallest radius) |
| local width ζ, calibrated (Study 12 table, w ≈ 1.07 px at 50 px) | ≈ 0.52 (8–80 px), < 0.5 (20–80 px) |
| P_ov | 0 |
| local Λ (L = 3.7, 7.5 µm) | 1.07 [0.68, 1.55], 1.00 [0.78, 1.60] |
| macroscopic anisotropy | n = 1: 0.029 (drift, cancels in Λ); n = 2–4 ≤ noise; Λ(φ) = 0.83–1.15 |

## Reading
A smooth, single-valued bubble wall (w ≈ 1 px at 23 µm) growing isotropically, Λ ≈ 1. Calibrated ζ ≈ 0.5 over 4–25 µm, lower than the Study 12 bubble. Only 9 frames and windows ≤ 25 µm: not well constrained.

## Caveats
- Scale and field from the file name; 9 frames only; one bubble; merging not analysed.
