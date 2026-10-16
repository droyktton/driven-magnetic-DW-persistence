# Study 14: `PtCoNiAl_10_4Oe_13mar16` (PtCoNiAl, 10.4 Oe, two creep walls, Δt = 20 s)

First report, 2026-10-16. Tools: [`arrival_map.py`](../README.md#arrival_mappy-arrival-map-of-a-movie-bubbles-several-walls) and [`contour_analysis.py`](../README.md#contour_analysispy-roughness-overhangs-and-kpz-term-from-an-arrival-map). Independent of the other studies.

## Data
- **Source:** group movie `PtCoNiAl__13Marzo2016__BGI180mA_H10_4Oe_20X_163umX122um_Dt20Seg_50fps_1000vecesRealTime.avi` (not in the repository).
  - 250 frames, 1392×1040 px, 8-bit compressed video.
  - From the file name: PtCoNiAl, 180 mA in the coil = **10.4 Oe**, 20X, field 163 × 122 µm, so **0.117 µm/px**; **Δt = 20 s** between frames (the video plays at 50 fps, 1000× real time), 83 min in total.
- **Geometry:** circular aperture. Two walls enter, from the lower left and from the top, sweep diagonally across the field and annihilate (frames ~20–220). Switching makes the pixels brighter.

```bash
python arrival_map.py ".../PtCoNiAl__13Marzo2016__BGI180mA_H10_4Oe_20X_163umX122um_Dt20Seg_50fps_1000vecesRealTime.avi" \
       --outdir PtCoNiAl_10_4Oe_13mar16 --detrend --qc 20,40,60,90,110,125,140,160,187
python contour_analysis.py PtCoNiAl_10_4Oe_13mar16/arrival.npz PtCoNiAl_10_4Oe_13mar16 --frames 20 200 --every 2 \
       --um-per-px 0.117 --roi-erode 60
```

## Arrival map and checks
- **Global brightness jump:** at frame 50 the whole image darkens by 39 grey levels (camera or lamp), followed by a slow upward drift. Without correction the step detector finds that jump instead of the switching. `--detrend` removes it (cumulative median of the frame-to-frame change).
- Otsu threshold on the step amplitude 24.1 (median 27.7 for switched pixels, 21.8 for the others: the switching contrast is weak); 76.6 % of the image switched.
- **Segmentation check:** `fig_qc_overlay.png` was approved by the user before the analysis. The walls are followed well in the interior; the darker rim of the aperture gives noisy detections, hence the ROI eroded by 60 px (71.6 % of the image).
- **Defects:** never-switched islands (14 500 px: an elongated island at the upper left and bright spots), contour points within 8 px of them dropped.
- 91 frames (every 2nd of 20–200), 1660 px of wall per frame in ~6 pieces.
- **Velocity:** 2.9 px per frame = 0.34 µm / 20 s = **17 nm/s**.

## Roughness and overhangs
| quantity | value |
|---|---|
| w at E = 10 / 20 / 50 px | 0.60 / 1.13 / 2.65 px |
| local width ζ, 8–80 px (0.9–9 µm), raw | **0.91** (same over 20–80 px) |
| pixel bias (synthetic bubbles, w ≈ 2.2 px at 50 px, Study 12) | 0.91 ↔ ζ ≈ 0.72 (8–80 px) or ≈ 0.84 (20–80 px) |
| P_ov (tol 3 px) at 2.3 / 5.9 / 12 µm | 0.04 / 0.11 / 0.08 |
| l₀^geom (P_ov = 1/2) | > 122 px (14 µm) |

`fig_contour_analysis.png`: w²(E) is a single power law from 0.5 to 14 µm, no break. The wall is mostly single-valued at these scales (P_ov ≤ 0.15).

## KPZ term
| L | dt = 1 frame | dt = 2 frames |
|---|---|---|
| 8 px (0.9 µm) | Λ = 1.08 [1.03, 1.14] | 0.94 [0.88, 1.01] |
| 16 px (1.9 µm) | 0.87 [0.75, 1.00] | 0.79 [0.69, 0.89] |
| 32 px (3.7 µm) | 1.03 [0.84, 1.29] | 1.04 [0.90, 1.19] |

R = |Λ|⟨m²⟩/2 ≤ 0.10. **Λ ≈ 1:** the wall moves along its normal with a speed that does not depend on its orientation (isotropic medium, geometric KPZ term), as in Studies 1–11. No sign of qKPZ (Λ ≫ 1).

## Reading
A creep wall at 10.4 Oe in PtCoNiAl, at the same resolution as Studies 2–11: ζ ≈ 0.9 raw at 1–9 µm (local width), ≈ 0.7–0.85 after the pixel-bias correction of Study 12 (synthetic bubbles of similar amplitude; the wall's w²(E) is straighter than theirs, so the correction is approximate), between qEW (1.25) and the creep equilibrium value (2/3), no visible l₀, Λ ≈ 1. Same picture as the Pt/Co/Pt and P7925 walls.

## Caveats
- Scale and field from the file name; compressed 8-bit video.
- Only 83 min of motion and two walls cut by the aperture into pieces: the local width is limited to 14 µm, so l₀ (if above) is not reached.
- Persistence, ξ(τ), χ4 not computed (steps 2–4 need h(x, t) of a single wall).
