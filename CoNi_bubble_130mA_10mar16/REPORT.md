# Study 12: `CoNi_bubble_130mA_10mar16` (Pt/(CoNi)×3/Ta, one bubble, 10 ms pulses)

First report, 2026-10-16. Tools: `arrival_map.py`, `contour_analysis.py`, `synthetic_bubble.py` (see the [README](../README.md)). Independent of the other studies.

## Data
- **Source:** group movie `Pt6(Co02Ni06)x3Ta5_10ms(t_total1s)130mABobinaChica 5XSinBKG 1080.mov` (background-subtracted version of `… 5X(convertido).mov`; not in the repository).
  - 99 frames, 1392×1040 px, compressed video.
  - From the file name: Pt(6)/(Co0.2Ni0.6)×3/Ta(5), field pulses of 10 ms, 1 s in total (so presumably **one frame per pulse**), 130 mA in the small coil (field not known), 5X.
  - **Scale 0.467 µm/px assumed** (5X of the same microscope: 650 µm / 1392 px in the PtCoNiAl movies).
- **Content:** one bubble nucleates at a point near the centre and grows from r ≈ 20 px to fill the aperture. A **planar defect**, a straight seam from the nucleation point towards the right, holds the wall back and cuts a notch into the bubble.

```bash
python arrival_map.py ".../Pt6(Co02Ni06)x3Ta5_10ms(t_total1s)130mABobinaChica 5XSinBKG 1080.mov" \
       --outdir CoNi_bubble_130mA_10mar16 --qc 10,20,30,40,49,55,60,70,85
python contour_analysis.py CoNi_bubble_130mA_10mar16/arrival.npz CoNi_bubble_130mA_10mar16 --frames 10 54 \
       --um-per-px 0.467 --centre 459.1 625.3 --exclude-deg -30 40
```

## Arrival map and checks
- Otsu threshold 89.9 (median step amplitude 145 switched, 37 not): clean contrast.
- **Segmentation check:** `fig_qc_overlay.png` was approved by the user. The contour follows the wall closely up to frame ~60; after that the wall enters the vignetted right side of the aperture and the detection becomes noisy. **Frames 10–54** used.
- **Centre:** (459.1, 625.3) px, centroid of the bubble at frame 9 (r_eq = 72 px).
- **Notch:** the defect slows the wall in the sector θ ≈ −15° … +25° (0° = +x, angles towards +y); the sector **−30° … +40° is excluded**. With the notch included the polar S(q) gives ζ = 0.89 and the height distribution is strongly skewed: both are the defect.
- Mean radius 79 → 408 px (37 → 190 µm), 1300 px of wall per frame.
- **Velocity:** 7.5 px per frame = **3.5 µm per pulse** (0.35 mm/s during the pulse, if the wall moves only during the pulse).

## Roughness and overhangs
| quantity | value |
|---|---|
| w at E = 10 / 20 / 50 px | 0.33 / 0.60 / 1.16 px (0.15 / 0.28 / 0.54 µm) |
| local width ζ, raw | 0.75 (8–80 px = 3.7–37 µm), 0.69 (20–80 px) |
| **local width ζ, calibrated** | **≈ 0.64 (8–80 px), ≈ 0.67 (20–80 px)** |
| P_ov | 0 up to 75 px (35 µm) |

The wall is only ~1 px rough, so the pixel staircase and the 1 px smoothing bias the local-width exponent. **Calibration** with synthetic bubbles of known ζ, same radii (79 → 408 px), amplitude w ≈ 1.05 px at 50 px, analysed with the same `contour_analysis.py` (2 seeds each, `synthetic_bubble.py Z SEED --out F.npz --w50 1.16 --r0 79 --r1 408`, `contour_analysis.py F.npz OUT --frames 1 320 --every 6 --centre 519.5 519.5 --skip-kpz`):

| ζ input | 0.50 | 0.63 | 0.80 | 1.00 | 1.25 |
|---|---|---|---|---|---|
| measured, 8–80 px | 0.65–0.67 | 0.74–0.75 | 0.87 | 1.00 | 1.16 |
| measured, 20–80 px | 0.54 | 0.65 | 0.79–0.80 | 0.96–0.97 | 1.19 |

Seed scatter ±0.01. The measured 0.75 and 0.69 map back to ζ ≈ 0.64 and ≈ 0.67. A qEW wall (1.25) would have given 1.16–1.19. With a twice larger amplitude (w ≈ 2.2 px at 50 px) the same estimator gives 0.77–0.79 for ζ = 0.5, 0.85 for 0.63, 0.96 for 0.8, 1.07 for 1.0 and 1.22 for 1.25: the bias depends on the amplitude, so each wall needs its own calibration.

## Isotropy and KPZ term
- **Macroscopic:** normal velocity vs orientation (radial advance between frames, allowed angles): |δv_n|/v = 0.058 for n = 1 (a drift of the centre, which cancels in λ = v + v″), 0.012, 0.004, 0.006 for n = 2–4, below the noise (0.025–0.04). Krug–Spohn **Λ(φ) = (v + v″)/v = 0.86–1.19** (n ≤ 4): isotropic growth.
- **Local (step10 B):** Λ = 0.43 [0.35, 0.52] at L = 3.7 µm, 0.46 [0.35, 0.59] at 7.5 µm, 1.36 [0.95, 1.91] at 15 µm; R ≤ 0.011. Λ = O(1), not ≫ 1.
- Height skewness (r(θ) minus modes n ≤ 4): +0.14, −0.55, +0.23 by thirds, frame-to-frame std 0.59: not significant with one bubble.

## Reading
An expanding bubble in a Pt/(CoNi) multilayer, wall smooth and single-valued, **ζ ≈ 0.64–0.67 over 4–37 µm**, isotropic growth and Λ = O(1). The exponent is close to the qKPZ value (0.63) and to the creep equilibrium value (2/3), but the KPZ term is not anomalous. A bubble at this resolution (≈ 0.5 µm/px, fits over tens of µm) gives ζ ≈ 0.6 without any sign of qKPZ dynamics. This matters for bubble experiments that read ζ ≈ 0.63 as qKPZ (e.g. Durin et al., PRB 110, L020405 (2024), CoFeB at 0.13–0.16 mT, 0.4 µm/px). This is not a resolution or scale effect: the walls of Studies 2, 4, 5, 8 and 10, coarse-grained to 0.47 µm/px and fitted over the same 3.7–37 µm with the same code, give ζ ≈ 0.9–1.35 and are two to three times rougher (control in [`UNIVERSALITY.md`](../UNIVERSALITY.md)). The bubble wall is genuinely smoother, with a smaller ζ; whether because of the expanding geometry, the fast pulsed drive or the sample is open.

## Caveats
- Scale assumed; field unknown; frame = pulse assumed from the file name.
- One bubble: the large-scale modes are one realisation (S(q) and the height distribution are not usable).
- Sector of the defect excluded by hand; the defect may also affect neighbouring angles.
