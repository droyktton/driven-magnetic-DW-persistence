# Study 6: `0_38A_20s_20x_51C_1` (0.38 A = 21.9 Oe, 51 °C, Δt = 20 s, 10 h)

First report, 2026-10-13. Pipeline and definitions: see the [README](../README.md). The question of the universality class below l₀ is summarised for all studies in [UNIVERSALITY.md](../UNIVERSALITY.md).

## Data

- Folder `…/MOBILEDATAMOKE/DataDWMoviesLargas/2Nov15/0_38A_20s_20x_51C_1/Pos0/` (not in the repository; analysed from a local copy, since the external disk resets often). Micro-Manager 1.4, recorded 2015-11-02 18:10 → 2015-11-03 04:09 (10.0 h).
- 1800 frames, 1392×1040 px, 1.2 s exposure, **Δt = 20 s** (`Interval_ms = 20000`), 20x objective.
- **0.38 A = 21.9 Oe at 51 °C** (thesis `Data/Info.txt`: coil 57.6 Oe/A).
- Sample, microscope and field of view: same as [Study 2](../0_8A_15s_20x_RT_1/REPORT.md), whose run was at 0.8 A, room temperature. The ring of defects and the large bubbles are the same ones.
- Scale 0.1175 µm/px, as in Study 2.
- Part of the field series of M. Grassi's thesis, with [Study 5](../0_36A_20s_20x_51C_2_1/REPORT.md) (0.36 A), [Study 6](../0_38A_20s_20x_51C_1/REPORT.md) (0.38 A) and [Study 7](../0_42A_20s_20x_52C_1/REPORT.md) (0.42 A). Each run starts from a new wall.

```bash
python step1_extract_wall.py ~/MOKEdata_local/2Nov15/0_38A_20s_20x_51C_1/Pos0 --um-per-px 0.1175 --flatfield 2 --rows 0 880
python step2_persistence.py 0_38A_20s_20x_51C_1 --row-mean --eps 0.5 1 1.5 2 3 5
python step7_roughness.py 0_38A_20s_20x_51C_1 --every 5 --rmin 25
python step8_crossover.py 0_36A_20s_20x_51C_2_1 0_38A_20s_20x_51C_1 0_42A_20s_20x_52C_1
python step10_universality.py 0_42A_20s_20x_52C_1 0_38A_20s_20x_51C_1 0_36A_20s_20x_51C_2_1 0_8A_15s_20x_RT_1 \
       PtCoPt_0_50A_2sep16 PtCoPt_0_48A_1sep16 magnetic_fliped
```

## Wall extraction and checks

- **Drift.** The sample drifts much more than in Study 2: dy −41.7…0 px, dx −11.7…+11.3 px (Study 2: up to 9 px). The registration follows it, with a scatter of 0.13–0.17 px around the smoothed curve.
- **Threshold of the swept area.** With the default Otsu threshold the normalisation region was empty: the halos of the moved defects dominate the contrast histogram, and the threshold fell at −0.32. Now, when Otsu marks more than half of the field, it is recomputed on the clipped contrast (README, step 1). The threshold is 0.034.
- `--flatfield 2`: without it, the whole initial domain showed switching contrast from the illumination change, and the swept area covered the left half of the field.
- `--rows 0 880`: noisy stripes on the bottom rows. Rows 12–879 are kept (868 rows, 102 µm).
- **Segmentation check:** `fig_qc_overlay.png`, `qc/overlay.mp4` (all frames), `qc/overlay_NNNN.png`. The line follows the boundary in all frames. Approved by the user before the analysis.
- **Motion:** total advance 38 px = 4.5 µm, **v = 0.02 px/frame = 0.12 nm/s**. No retreats > 1 px.
- **Noise:** σ_Δh = 0.34 px, so ε = 3σ_Δh ≈ 1 px is the reference threshold.

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted (`step2 --row-mean`). τ in frames (1 frame = 20 s). Figures: `fig_xi_chi4_sub_rm_eps<ε>.png`, `fig_Cn_tau_sub_rm_eps<ε>.png`. ξ_max is taken over τ < 0.8 τ_max.

| ε [px] (µm) | Π(τ=1) | τ at Π = ½ | **τ\*** (max χ4) | χ4(τ\*) | **ξ(τ\*)** [px] (µm) | ξ_max [px] (µm) at τ |
|---|---|---|---|---|---|---|
| 0.5 (0.06) | 0.89 | 96 | 203 (68 min) | 32.8 | 30.7 ± 1.3 (3.60) | 52.4 (6.16) at 450 |
| **1 (0.12)** | 0.97 | 141 | **228 (76 min)** | 31.5 | **30.8 ± 1.1 (3.62)** | 52.9 (6.22) at 460 |
| 1.5 (0.18) | 1.00 | 188 | 240 (80 min) | 34.2 | 34.4 ± 1.1 (4.04) | 53.6 (6.30) at 475 |
| 2 (0.23) | 1.00 | 210 | 261 (87 min) | 31.3 | 34.5 ± 1.0 (4.06) | 50.6 (5.94) at 473 |
| 3 (0.35) | 1.00 | 255 | 292 (97 min) | 27.2 | 34.0 ± 0.9 (3.99) | 44.8 (5.26) at 487 |
| 5 (0.59) | 1.00 | 330 | 427 (142 min) | 23.7 | 36.3 ± 1.2 (4.26) | 40.6 (4.77) at 490 |

## Roughness and l₀

`step7_roughness.py`, `step8_crossover.py` (`fig_roughness.png`, `fig_local_width.png`, `fig_local_sq.png`, `fig_crossover.png`).

- Base plane θ0 = −17.8°; 895 rows (105 µm); overhangs 0.9 %, almost no folded segments.
- **The method changes ζ:** local width with rotation over 3–26 µm: contour 0.71, folds excluded 0.80, u_area 0.79. Global S(q): 0.87 ± 0.08 (thirds 1.08, 0.74, 0.46).
- **l₀ is small:** 13 µm from the local width, 8.5 µm from the global S(q), 5 µm model free.
- w²(ℓ) rises steeply up to ~4 µm, then with ζ ≈ 0.9–1.0 up to ~13 µm, then ≈ 0.5.
- Height distribution: width 2.8 µm, skewness −0.06, excess kurtosis −1.4 (flat-topped: the wall barely moves and keeps its large-scale shape).

## qEW or qKPZ below l₀

`step10_universality.py` (`fig_universality.png`).

**A.** ζ over [4 µm, l₀], folds excluded: 0.93 (0.91–0.97); u_area 0.95–1.04. The broken power law does not describe w²(ℓ): the χ² is poor for every model, and the free fit puts l₀ at the upper bound with ζ₁ = 0.76. The range between 4 µm and l₀ is less than a factor 3. **This is the weakest case.** The local exponent is still ≈ 0.9–1.0, not 0.63.

**B.** v_n(φ) is noisy: the wall moves only 4.5 µm in 10 h, so few segments switch. Λ scatters between −1.8 and 1.1 (contour median 0.24; h(x) ≈ 0), R ≤ 0.12. No Λ ≫ 1.

## Not run yet

Steps 3–6 and 9 (ε sweep figure, avalanches, pinning map, pinning tests, wavelets).
