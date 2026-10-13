# Study 7: `0_42A_20s_20x_52C_1` (0.42 A = 24.2 Oe, 52 °C, Δt = 20 s, 10 h)

First report, 2026-10-13. Pipeline and definitions: see the [README](../README.md). The question of the universality class below l₀ is summarised for all studies in [UNIVERSALITY.md](../UNIVERSALITY.md).

## Data

- Folder `…/MOBILEDATAMOKE/DataDWMoviesLargas/3Nov15/0_42A_20s_20x_52C_1/Pos0/` (not in the repository; analysed from a local copy, since the external disk resets often). Micro-Manager 1.4, recorded 2015-11-03 16:53 → 2015-11-04 02:53 (10.0 h).
- 1800 frames, 1392×1040 px, 1.2 s exposure, **Δt = 20 s** (`Interval_ms = 20000`), 20x objective.
- **0.42 A = 24.2 Oe at 52 °C** (thesis `Data/Info.txt`: coil 57.6 Oe/A).
- Sample, microscope and field of view: same as [Study 2](../0_8A_15s_20x_RT_1/REPORT.md), whose run was at 0.8 A, room temperature. The ring of defects and the large bubbles are the same ones.
- Scale 0.1175 µm/px, as in Study 2.
- Part of the field series of M. Grassi's thesis, with [Study 5](../0_36A_20s_20x_51C_2_1/REPORT.md) (0.36 A), [Study 6](../0_38A_20s_20x_51C_1/REPORT.md) (0.38 A) and [Study 7](../0_42A_20s_20x_52C_1/REPORT.md) (0.42 A). Each run starts from a new wall.

```bash
python step1_extract_wall.py ~/MOKEdata_local/3Nov15/0_42A_20s_20x_52C_1/Pos0 --um-per-px 0.1175 --rows 0 640
python step2_persistence.py 0_42A_20s_20x_52C_1 --row-mean --eps 0.5 1 1.5 2 3 5
python step7_roughness.py 0_42A_20s_20x_52C_1 --every 5 --rmin 25
python step8_crossover.py 0_36A_20s_20x_51C_2_1 0_38A_20s_20x_51C_1 0_42A_20s_20x_52C_1
python step10_universality.py 0_42A_20s_20x_52C_1 0_38A_20s_20x_51C_1 0_36A_20s_20x_51C_2_1 0_8A_15s_20x_RT_1 \
       PtCoPt_0_50A_2sep16 PtCoPt_0_48A_1sep16 magnetic_fliped
```

## Wall extraction and checks

- **Drift.** The sample drifts much more than in Study 2: dy −32.9…+0.3 px, dx −18.2…+0.2 px (Study 2: up to 9 px). The registration follows it, with a scatter of 0.13–0.17 px around the smoothed curve.
- Default options, as in Study 2.
- `--rows 0 640`: below y ≈ 640 the initial wall reaches the left edge of the field (x ≈ 0), and a spurious patch appears there. Without the restriction, frames 0–374 were dropped. Rows 13–639 are kept (627 rows, 74 µm), and all 1800 frames are used.
- **Segmentation check:** `fig_qc_overlay.png`, `qc/overlay.mp4` (all frames), `qc/overlay_NNNN.png`. The line follows the boundary in all frames. Approved by the user before the analysis.
- **Motion:** total advance 176 px = 21 µm, **v = 0.10 px/frame = 0.57 nm/s**. No retreats > 1 px.
- **Noise:** σ_Δh = 0.35 px, so ε = 3σ_Δh ≈ 1 px is the reference threshold.

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted (`step2 --row-mean`). τ in frames (1 frame = 20 s). Figures: `fig_xi_chi4_sub_rm_eps<ε>.png`, `fig_Cn_tau_sub_rm_eps<ε>.png`. ξ_max is taken over τ < 0.8 τ_max.

| ε [px] (µm) | Π(τ=1) | τ at Π = ½ | **τ\*** (max χ4) | χ4(τ\*) | **ξ(τ\*)** [px] (µm) | ξ_max [px] (µm) at τ |
|---|---|---|---|---|---|---|
| 0.5 (0.06) | 0.84 | 18 | 31 (10 min) | 9.1 | 15.8 ± 0.3 (1.85) | 16.6 (1.95) at 47 |
| **1 (0.12)** | 0.94 | 28 | **47 (16 min)** | 8.8 | **18.2 ± 0.4 (2.14)** | 18.6 (2.19) at 59 |
| 1.5 (0.18) | 0.98 | 37 | 48 (16 min) | 12.7 | 18.7 ± 0.4 (2.20) | 19.4 (2.28) at 67 |
| 2 (0.23) | 0.99 | 43 | 53 (18 min) | 12.3 | 19.0 ± 0.4 (2.23) | 19.6 (2.31) at 71 |
| 3 (0.35) | 0.99 | 52 | 63 (21 min) | 11.9 | 19.3 ± 0.4 (2.27) | 20.1 (2.36) at 82 |
| 5 (0.59) | 1.00 | 70 | 82 (27 min) | 11.3 | 19.5 ± 0.5 (2.30) | 19.8 (2.33) at 97 |

## Roughness and l₀

`step7_roughness.py`, `step8_crossover.py` (`fig_roughness.png`, `fig_local_width.png`, `fig_local_sq.png`, `fig_crossover.png`).

- Base plane θ0 = −7.2°; the tilt changes from −12.8° to −3.3° during the run. 626 rows (74 µm); overhangs 0.5 %.
- **ζ:** local width with rotation over 3–18 µm: contour 1.04, folds excluded 1.08, u_area 1.07. Global S(q): **1.16 ± 0.04** (thirds 1.21, 1.12, 1.16, stable).
- **l₀:** 14.5 µm from the local width, 15 µm from the global S(q), 8 µm model free. l_tilt = 27 µm.
- Height distribution: width 3.2 µm, skewness −0.38.

## qEW or qKPZ below l₀

`step10_universality.py` (`fig_universality.png`).

**A.** ζ over [4 µm, l₀], folds excluded: 1.04 (1.02–1.06). Broken power law: ζ₁ = 1.18 (1.14–1.28), or 1.24 on u_area. χ² qEW / qKPZ = 32 / 675 (u_area 11 / 616).

**B.** Λ = 0.5–1.0 (contour), 0.1–0.4 (h(x)); R ≤ 0.12. Λ is between the elastic and the isotropic values.

## Not run yet

Steps 3–6 and 9 (ε sweep figure, avalanches, pinning map, pinning tests, wavelets).
