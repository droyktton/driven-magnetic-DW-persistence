# Study 5: `0_36A_20s_20x_51C_2_1` (0.36 A = 20.7 Oe, 51 °C, Δt = 20 s, 10 h)

First report, 2026-10-13. Pipeline and definitions: see the [README](../README.md). The question of the universality class below l₀ is summarised for all studies in [UNIVERSALITY.md](../UNIVERSALITY.md).

## Data

- Folder `…/MOBILEDATAMOKE/DataDWMoviesLargas/16Nov15/0_36A_20s_20x_51C_2_1/Pos0/` (not in the repository; analysed from a local copy, since the external disk resets often). Micro-Manager 1.4, recorded 2015-11-16 16:56 → 2015-11-17 02:55 (10.0 h).
- 1800 frames, 1392×1040 px, 1.2 s exposure, **Δt = 20 s** (`Interval_ms = 20000`), 20x objective.
- **0.36 A = 20.7 Oe at 51 °C** (thesis `Data/Info.txt`: coil 57.6 Oe/A).
- Sample, microscope and field of view: same as [Study 2](../0_8A_15s_20x_RT_1/REPORT.md), whose run was at 0.8 A, room temperature. The ring of defects and the large bubbles are the same ones.
- Scale 0.1175 µm/px, as in Study 2.
- Part of the field series of M. Grassi's thesis, with [Study 5](../0_36A_20s_20x_51C_2_1/REPORT.md) (0.36 A), [Study 6](../0_38A_20s_20x_51C_1/REPORT.md) (0.38 A) and [Study 7](../0_42A_20s_20x_52C_1/REPORT.md) (0.42 A). Each run starts from a new wall.

```bash
python step1_extract_wall.py ~/MOKEdata_local/16Nov15/0_36A_20s_20x_51C_2_1/Pos0 --um-per-px 0.1175 --flatfield 2 --swept-thr 0.02 --rows 0 880
python step2_persistence.py 0_36A_20s_20x_51C_2_1 --row-mean --eps 0.5 1 1.5 2 3 5
python step7_roughness.py 0_36A_20s_20x_51C_2_1 --every 5 --rmin 25
python step8_crossover.py 0_36A_20s_20x_51C_2_1 0_38A_20s_20x_51C_1 0_42A_20s_20x_52C_1
python step10_universality.py 0_42A_20s_20x_52C_1 0_38A_20s_20x_51C_1 0_36A_20s_20x_51C_2_1 0_8A_15s_20x_RT_1 \
       PtCoPt_0_50A_2sep16 PtCoPt_0_48A_1sep16 magnetic_fliped
```

## Wall extraction and checks

- **Drift.** The sample drifts much more than in Study 2: dy −39.8…+0.3 px, dx −5.4…+24.1 px (Study 2: up to 9 px). The registration follows it, with a scatter of 0.13–0.17 px around the smoothed curve.
- **The swept band is very faint.** Its contrast is (B−D)/B ≈ 0.037 against a background of ±0.01. Otsu's threshold, even on the clipped contrast, falls on the halos of the defects (0.08–0.16), so the swept area was a single halo. The threshold was set by hand: `--swept-thr 0.02`.
- `--flatfield 2`: the illumination changes non-uniformly during the run (it is not removed by the per-frame normalisation).
- `--rows 0 880`: noisy stripes on the bottom rows of the camera. Rows 13–879 are kept (867 rows, 102 µm).
- **Segmentation check:** `fig_qc_overlay.png`, `qc/overlay.mp4` (all frames), `qc/overlay_NNNN.png`. The line follows the boundary in all frames. Approved by the user before the analysis.
- **Motion:** total advance 40 px = 4.7 µm, **v = 0.02 px/frame = 0.13 nm/s**. No retreats > 1 px.
- **Noise:** σ_Δh = 0.33 px, so ε = 3σ_Δh ≈ 1 px is the reference threshold.

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted (`step2 --row-mean`). τ in frames (1 frame = 20 s). Figures: `fig_xi_chi4_sub_rm_eps<ε>.png`, `fig_Cn_tau_sub_rm_eps<ε>.png`. ξ_max is taken over τ < 0.8 τ_max.

| ε [px] (µm) | Π(τ=1) | τ at Π = ½ | **τ\*** (max χ4) | χ4(τ\*) | **ξ(τ\*)** [px] (µm) | ξ_max [px] (µm) at τ |
|---|---|---|---|---|---|---|
| 0.5 (0.06) | 0.89 | 76 | 167 (56 min) | 35.3 | 14.1 ± 0.5 (1.66) | 15.1 (1.78) at 243 |
| **1 (0.12)** | 0.98 | 120 | **216 (72 min)** | 32.7 | **19.2 ± 0.4 (2.26)** | 19.3 (2.27) at 233 |
| 1.5 (0.18) | 1.00 | 162 | 238 (79 min) | 34.0 | 20.4 ± 0.4 (2.40) | 20.5 (2.41) at 245 |
| 2 (0.23) | 1.00 | 183 | 242 (81 min) | 31.1 | 21.1 ± 0.4 (2.48) | 21.4 (2.51) at 253 |
| 3 (0.35) | 1.00 | 227 | 254 (85 min) | 26.6 | 22.5 ± 0.5 (2.64) | 26.6 (3.12) at 373 |
| 5 (0.59) | 1.00 | 312 | 268 (89 min) | 20.8 | 21.1 ± 0.6 (2.47) | 28.1 (3.30) at 496 |

## Roughness and l₀

`step7_roughness.py`, `step8_crossover.py` (`fig_roughness.png`, `fig_local_width.png`, `fig_local_sq.png`, `fig_crossover.png`).

- Base plane θ0 = −16.9°; 882 rows (104 µm) in all frames; overhangs on 0.3 % of the (row, frame) pairs, and almost no folded segments.
- **ζ ≈ 1.06–1.08 over 3–26 µm** (local width with rotation: contour 1.06, folds excluded 1.08, u_area 1.06). Global S(q): 0.95 ± 0.08 (thirds 1.18, 0.71, 0.89).
- **l₀ (step8):** 80 µm from the local width (64 µm model free). The global S(q) does not bracket it (lower bound). The tilt spread halves at l_tilt = 60 µm.
- **The largest l₀ of all studies.** Below it, the wall is self-affine with ζ = 1.15–1.25 over more than a decade.
- Height distribution on the base plane: width 5.1 µm, skewness +0.85. The wall has a large step (a nearly horizontal "ledge" at y ≈ 500 px), which skews the global distribution.

## qEW or qKPZ below l₀

`step10_universality.py` (`fig_universality.png`).

**A.** ζ over [4 µm, l₀], folds excluded: 1.15 (1.13–1.17). Broken power law: **ζ₁ = 1.26 (1.24–1.28)**. χ² qEW / qKPZ = 1108 / 10651; both are poor fits, because w²(ℓ) is not exactly a broken power law, but qEW is ten times better.

**B.** Λ ≈ 1.0–1.5 (contour median 1.25, h(x) median 1.21), R ≤ 0.23. A few points at L = 16 µm have Λ up to 3.5, with ±2–8 errors. Λ is about twice its value at 0.42 A (Study 7), while v is 4.4× smaller: a fixed λ (qKPZ) would make it 4.4× larger. The value is compatible with the isotropic one (Λ = 1).

## Not run yet

Steps 3–6 and 9 (ε sweep figure, avalanches, pinning map, pinning tests, wavelets).
