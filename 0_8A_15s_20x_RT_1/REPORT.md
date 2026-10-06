# Study 2: `0_8A_15s_20x_RT_1` (MOKE TIFF sequence, Δt = 20 s, 10 h)

First report, 2026-10-06. Pipeline and definitions: see the [README](../README.md). This study is independent of [Study 1](../magnetic_fliped/REPORT.md).

## Data

- Micro-Manager 1.4 acquisition, recorded 2015-11-19 (17:15:57 → 03:15:37), folder `…/ResumenMaestría_Mati/0_8A_15s_20x_RT_1/Pos0/`, not in the repository.
- 1800 frames, 1392×1040 px, 14-bit (uint16), 1 s exposure, 20x objective, room temperature; the name suggests a drive of 0.8 A.
- **Δt = 20 s, not 15 s as the folder name says.** The metadata `Interval_ms = 20000`, the per-frame elapsed times (1791 of 1799 intervals are 20.0 s) and the wall-clock stamps (35 980 s over 1799 intervals) agree. Total duration 10.0 h.
- **Scale δ = 0.1175 µm/px** (field of view ≈ 164 × 122 µm). The metadata has no pixel size (`PixelSize_um = 0`); the value comes from the thesis material of this same measurement (M. Grassi, master's thesis, `…/ResumenMaestría_Mati/`), by two independent routes:
  1. The thesis avalanche analysis of this measurement discards events smaller than 20 px and states that avalanches below S ≈ 0.3 µm² cannot be resolved; its area histograms (`Data/Histogramas_lineal`, `Data/Histogramas_log`) start at 0.27–0.28 µm². So 20 px ≈ 0.28 µm², i.e. 0.014 µm²/px → 0.118 µm/px.
  2. The scale bars of the full 1392×1040 frames in the defense slides (slides 11 and 13, 100 µm = 213 px) give 0.47 µm/px. Those images belong to the velocity measurements, done with the 5x objective; scaled to the 20x, 0.47/4 = 0.1175 µm/px.

  This is **not** the scale of Study 1 (0.17 µm/px). The uncertainty is a few %; route 1 rests on the histogram binning, route 2 on the nominal 4x ratio of the objectives.
- Geometry: one wall, roughly vertical but tilted (θ ≈ −10° to −14° from vertical), with the dark domain on the left growing to the right. Many static round defects (bubbles/dust), several of them on the wall path.

```bash
python step1_extract_wall.py /media/…/0_8A_15s_20x_RT_1/Pos0 --um-per-px 0.1175
python step2_persistence.py 0_8A_15s_20x_RT_1                       # eps = 0.97 (3σ_Δh), 1, 2
python step2_persistence.py 0_8A_15s_20x_RT_1 --eps 0.5 1.5 3 5
python step3_eps_sweep.py 0_8A_15s_20x_RT_1
python step2_persistence.py 0_8A_15s_20x_RT_1 --row-mean --eps 0.5 0.97 1.5 2 3 5   # results below
python step3_eps_sweep.py 0_8A_15s_20x_RT_1 --row-mean
python step2_persistence.py 0_8A_15s_20x_RT_1 --int                  # cross-check
```

## Wall extraction and checks

Details in the README (step 1, *TIFF folder input*). Points specific to this data set:

- **Sample drift.** The sample moves rigidly by up to dy = −9.3 px and dx = −4.4…+3.5 px over the 10 h (`fig_drift_arrival.png`, left). Uncorrected, the dx drift alone would be a slow fake displacement equal to ~16 frames of real wall motion. Every frame is registered to frame 0; the measured drift scatters by only 0.03–0.04 px around its smoothed curve.
- **Crop.** Field-of-view edge on the right: columns x ≥ 1197 discarded. Rows 12–1017 kept (the vertical drift leaves the edges without data). The wall stays between x ≈ 150 and 950, inside the crop, in all frames: **all 1800 frames are used**.
- **Contrast.** The domain contrast is only ~5 % of the raw intensity, below the vignetting. Per-pixel bright/dark references cancel it, and also cancel the static defects (8841 px of defects inside the swept area get the arrival time of their neighbours).
- **Overhangs.** The wall folds around defects: on average 54 of 1006 rows cross the wall more than once. The position used for the analysis is the switched area per row, x_eff(y,t). It equals the front where the wall is single-valued (91.6 % of row-frames). Using the front instead gave 877 one-frame jumps > 20 px, almost all of them pinch-offs of thin unswitched channels between defects with no real switching (checked on ratio images); with x_eff there are 116, and the largest ones checked are real avalanches (compact patches of 30–50 px switching within one frame).
- **Segmentation check.** The overlays (`fig_qc_overlay.png`, `qc/overlay.mp4` with all 1800 frames, `qc/overlay_NNNN.png` every 50 frames) follow the wall, including next to the defects. The method was approved on a 1-in-10 subset before the full run.
- **Motion.** Total advance 509 px = 60 µm, mean velocity **0.28 px/frame = 1.65 nm/s**, no retreats > 1 px (`fig_h_mean.png`). The motion is intermittent, with avalanches of up to 44 px (5 µm) in one row and one frame (`fig_kymograph.png`).
- **Noise.** σ_Δh = 0.32 px = 38 nm (from backward one-frame steps), so the default ε = 3σ_Δh = 0.97 px = 0.11 µm.
- **Tilt.** θ goes from −12.8° to −9.8° (`fig_tilt_overhang.png`). h is measured along x and the distance n along y; the corrections to a normal displacement (cos θ) and to a distance along the wall (1/cos θ) are ≈ 2 % and are not applied below.

## Results

C(n,τ) is computed with the mean persistence of each row subtracted (`step2 --row-mean`, files `*_sub_rm_*`), C(n) = ⟨(p_i − m_i)(p_{i+n} − m_{i+n})⟩ with m_i = ⟨p_i⟩_t. The reason is explained under *Row-mean subtraction* below. χ4 and Π do not change; ξ does, slightly.

| ε [px] (µm) | Π(τ=1) | τ at Π = ½ | **τ\*** (max χ4) | χ4(τ\*) | **ξ(τ\*)** [px] (µm) | max ξ [px] (µm) at τ |
|---|---|---|---|---|---|---|
| 0.5 (0.06) | 0.78 | 6 | 5 (100 s) | 4.6 | 12.2 ± 0.2 (1.43) | 21.4 (2.51) at 80 |
| **0.97 (0.11)** | 0.89 | 9 | **7 (2.3 min)** | 5.1 | **13.5 ± 0.1 (1.59)** | 22.1 (2.60) at 89 |
| 1.5 (0.18) | 0.94 | 13 | 11 (3.7 min) | 6.0 | 15.2 ± 0.2 (1.79) | 23.3 (2.74) at 89 |
| 2 (0.24) | 0.96 | 14 | 14 (4.7 min) | 6.0 | 15.8 ± 0.2 (1.86) | 23.2 (2.73) at 94 |
| 3 (0.35) | 0.98 | 18 | 17 (5.7 min) | 6.0 | 16.3 ± 0.3 (1.92) | 23.2 (2.73) at 99 |
| 5 (0.59) | 0.99 | 25 | 24 (8 min) | 6.4 | 17.2 ± 0.3 (2.02) | 23.6 (2.77) at 112 |

τ in frames (1 frame = 20 s). Errors are statistical (fit). Figures: `fig_xi_chi4_sub_rm_eps<ε>.png`, `fig_Cn_tau_sub_rm_eps<ε>.png`, `fig_eps_sweep_sub_rm.png`; the same without row-mean subtraction: `*_sub_eps<ε>*`, `fig_eps_sweep_sub.png`.

### Main findings

1. **τ\* is resolved.** Unlike Study 1, χ4(τ) rises, peaks and decays. At ε = 3σ_Δh = 0.97 px, **τ\* = 7 frames ≈ 2.3 min**, with **ξ(τ\*) = 13.5 ± 0.1 (stat.) ± 0.3 (fit) px = 1.59 ± 0.04 µm**. Over ε = 0.5–5 px, ξ(τ\*) = 1.4–2.0 µm.
2. **τ\* is set by Π(τ\*) ≈ ½.** For every ε, τ\* is at or 1–3 frames before the lag at which half of the rows have moved by more than ε. So τ\* grows with ε (5 → 24 frames for ε = 0.5 → 5 px) and mostly measures the time to move a distance ε, as expected when χ4 ≈ Π(1−Π) × (correlated length). The peak χ4 itself (4.6–6.4) and ξ(τ\*) (12–17 px) depend more weakly on ε.
3. **The correlation length keeps growing well beyond τ\*.** ξ(τ) increases from ~8 px at τ = 1 to a maximum of **ξ_max ≈ 21–24 px = 2.5–2.8 µm at τ ≈ 80–110 frames (27–37 min)**, then decays as Π → 0. ξ_max changes by only ±1 px over a factor 10 in ε, so **ξ_max and its time scale of ~30 min are the most robust results** of this study, together with τ\*. The raw χ4 peaks much earlier only because its amplitude Π(1−Π) falls.
4. **χ4 and ξ are consistent at short τ.** C(n,τ) is exponential: (C−B)/A is a straight line on the semilog plot down to ~2 % of its amplitude. For τ ≲ 30 frames, the local normalized susceptibility (χ4 − B·L)/[Π(1−Π)] follows 2ξ, and the sum of the fitted model A·coth(1/2ξ) + B·L reproduces the direct χ4 within 4–20 % at τ\* (4 % at ε = 0.97 px, 20 % at ε = 5 px). For τ ≳ 50 frames the model sum falls below χ4: correlations at long lags are not a single exponential.
5. **Integer h** (`--int`, cross-check, no row-mean subtraction): ε = 1 px gives τ\* = 6 frames, ξ(τ\*) = 15.6 ± 0.3 px; ε = 2 px gives τ\* = 13 frames, ξ(τ\*) = 16.7 ± 0.3 px, consistent with the subpixel results.

### Row-mean subtraction

Without it, C(n) = ⟨p_i p_{i+n}⟩ − Π² is **negative at n ≈ 500–1000 px** (top half against bottom half of the wall). C(n) then contains a static term ⟨m_i m_{i+n}⟩ − Π², from rows that advance at different average rates over the 10 h (in the kymograph, the top overtakes the bottom late in the run; some rows are pinned by defects for a long time). This term sums to exactly zero over all pairs, so it is positive at short n and must be negative at large n. Subtracting m_i removes it, and the identity χ4 = Σ_n (1−|n|/L)·C(n) still holds exactly, because χ4 = (1/L) Σ_ij cov_t(p_i, p_j).

Effect at ε = 0.97 px:

| | without | with row-mean subtraction |
|---|---|---|
| C(n = 500–1000)·L at τ\* | −21.7 | −2.2 |
| ξ(τ\*): fit n ≤ 300 / n ≤ 1000 / B = 0 | 13.3 / 16.9 / 13.9 | 13.5 / 13.5 / 13.2 |
| ξ(τ = 90): fit n ≤ 300 / n ≤ 1000 / B = 0 | 22.7 / 23.0 / 21.3 | 22.0 / 19.9 / 19.9 |
| model sum vs χ4 at τ\* | 8.9 vs 5.1 | 4.9 vs 5.1 |

So the fit-range ambiguity of ξ(τ\*) (±2 px before) disappears. At long lags (τ ≈ 90) a ~±1 px dependence remains, where Π ≈ 0.07 and fewer rows are persistent. The fitted plateau B·L is now small and slightly negative (−1 to −3), as expected from estimating m_i from a finite number of start times.

### Caveats

- **Resolution.** ξ(τ\*) ≈ 1.6 µm is close to the optical resolution of the 20x objective (~1 µm, NA 0.4, thesis) and spans ~14 px; ξ is resolved, but its smallest values (ξ ≈ 0.5–0.9 µm at τ = 1, depending on ε) are at the resolution limit.
- **ξ is measured along y**, not along the wall; the 1/cos θ correction (≈ 2 %) is not applied.
- `fig_eps_sweep_*.png`, panel 3, shows ξ and χ4 at τ = 1 frame, which is meaningful for Study 1 (τ\* ≤ 1 frame) but not here, where τ\* = 5–24 frames; the table above gives the values at τ\*.

## Next steps

- Map the per-row mean persistence m_i (and mean velocity per row) against the defect positions: the static heterogeneity removed by the row-mean subtraction is itself a measure of the pinning landscape.
- Adapt step 3 to show ξ(τ\*) and ξ_max against ε instead of the τ = 1 values.
- Avalanche statistics: the 20 s resolution and the clean x_eff(y,t) make a direct measurement of avalanche sizes and lateral extents feasible (Δx_eff between consecutive frames, clusters of moved rows).
