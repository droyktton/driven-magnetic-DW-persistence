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
python step2_persistence.py 0_8A_15s_20x_RT_1 --eps 0.5 0.97 1.5 2 3 5
python step3_eps_sweep.py 0_8A_15s_20x_RT_1
python step2_persistence.py 0_8A_15s_20x_RT_1 --row-mean --eps 0.5 0.97 1.5 2 3 5   # results below
python step3_eps_sweep.py 0_8A_15s_20x_RT_1 --row-mean
python step2_persistence.py 0_8A_15s_20x_RT_1 --int                  # cross-check
```

## Wall extraction and checks

**Movie of the whole measurement:** [`overlay_movie.gif`](overlay_movie.gif) (10 h, every 5th frame = one image every 100 s, drift-corrected, detected wall in red, 20 µm scale bar; 33 MB). Made with `python make_overlay_gif.py 0_8A_15s_20x_RT_1 /media/…/0_8A_15s_20x_RT_1/Pos0 --every 5`. The full-resolution movie of all 1800 frames is `qc/overlay.mp4`.

Details in the README (step 1, *TIFF folder input*). Points specific to this data set:

- **Sample drift.** The sample moves rigidly by up to dy = −9.3 px and dx = −4.4…+3.5 px over the 10 h (`fig_drift_arrival.png`, left). Uncorrected, the dx drift alone would be a slow fake displacement equal to ~16 frames of real wall motion. Every frame is registered to frame 0; the measured drift scatters by only 0.03–0.04 px around its smoothed curve.
- **Crop.** Field-of-view edge on the right: columns x ≥ 1197 discarded. Rows 12–1017 kept (the vertical drift leaves the edges without data). The wall stays between x ≈ 150 and 950, inside the crop, in all frames: **all 1800 frames are used**.
- **Contrast.** The domain contrast is only ~5 % of the raw intensity, below the vignetting. Per-pixel bright/dark references cancel it, and also cancel the static defects (8841 px of defects inside the swept area get the arrival time of their neighbours).
- **Overhangs.** The wall folds around defects: on average 51 of 1006 rows cross the wall more than once. The position used for the analysis is the switched area per row, x_eff(y,t). It equals the front where the wall is single-valued (91.8 % of row-frames). Using the front instead gave 877 one-frame jumps > 20 px, almost all of them pinch-offs of thin unswitched channels between defects with no real switching (checked on ratio images); with x_eff there are 136, and the largest ones checked are real avalanches (compact patches of 30–50 px switching within one frame).
- **Segmentation check.** The overlays (`fig_qc_overlay.png`, `qc/overlay.mp4` with all 1800 frames, `qc/overlay_NNNN.png` every 50 frames) follow the wall, including next to the defects. The method was approved on a 1-in-10 subset before the full run.
- **Motion.** Total advance 509 px = 60 µm, mean velocity **0.28 px/frame = 1.65 nm/s**, no retreats > 1 px (`fig_h_mean.png`). The motion is intermittent, with avalanches of up to 36 px (4 µm) in one row and one frame (`fig_kymograph.png`).
- **Noise.** σ_Δh = 0.33 px = 39 nm (from backward one-frame steps), so ε = 3σ_Δh ≈ 0.97 px = 0.11 µm is used as the reference threshold.
- **Tilt.** θ goes from −12.8° to −9.8° (`fig_tilt_overhang.png`). h is measured along x and the distance n along y; the corrections to a normal displacement (cos θ) and to a distance along the wall (1/cos θ) are ≈ 2 % and are not applied below.

## Results

C(n,τ) is computed with the mean persistence of each row subtracted (`step2 --row-mean`, files `*_sub_rm_*`), C(n) = ⟨(p_i − m_i)(p_{i+n} − m_{i+n})⟩ with m_i = ⟨p_i⟩_t. The reason is explained under *Row-mean subtraction* below. χ4 and Π do not change; ξ does, slightly.

| ε [px] (µm) | Π(τ=1) | τ at Π = ½ | **τ\*** (max χ4) | χ4(τ\*) | **ξ(τ\*)** [px] (µm) | max ξ [px] (µm) at τ |
|---|---|---|---|---|---|---|
| 0.5 (0.06) | 0.77 | 6 | 6 (2 min) | 4.7 | 12.6 ± 0.2 (1.48) | 21.2 (2.49) at 81 |
| **0.97 (0.11)** | 0.89 | 10 | **7 (2.3 min)** | 5.2 | **13.4 ± 0.1 (1.58)** | 21.9 (2.57) at 88 |
| 1.5 (0.18) | 0.94 | 13 | 11 (3.7 min) | 6.1 | 15.1 ± 0.2 (1.78) | 23.3 (2.73) at 91 |
| 2 (0.23) | 0.96 | 14 | 14 (4.7 min) | 6.0 | 15.7 ± 0.2 (1.85) | 23.2 (2.72) at 93 |
| 3 (0.35) | 0.98 | 18 | 19 (6.3 min) | 6.0 | 16.8 ± 0.3 (1.97) | 23.1 (2.72) at 99 |
| 5 (0.59) | 0.99 | 25 | 25 (8.3 min) | 6.4 | 17.4 ± 0.3 (2.05) | 23.5 (2.77) at 115 |

τ in frames (1 frame = 20 s). Errors are statistical (fit). Figures: `fig_xi_chi4_sub_rm_eps<ε>.png`, `fig_Cn_tau_sub_rm_eps<ε>.png`, `fig_eps_sweep_sub_rm.png`; the same without row-mean subtraction: `*_sub_eps<ε>*`, `fig_eps_sweep_sub.png`.

### ξ(τ\*) and ξ_max: two points of the same curve

The table gives two correlation lengths, both read from the curve ξ(τ):

- **ξ(τ\*)** is ξ at τ\*, the lag of the **maximum of χ4**. Since χ4 ≈ (number of correlated sites) × (variance of one site), and that variance Π(1−Π) is largest at Π = ½, τ\* falls where half of the wall has moved by more than ε. It therefore **depends on ε** (τ\* = 2 min for ε = 0.5 px, 8 min for ε = 5 px). Here ξ(τ\*) ≈ 1.6 µm at ε = 0.97 px (1.5–2.1 µm over ε). It is the standard quantity of the dynamic-heterogeneity literature ("the correlation length at the peak of χ4").
- **ξ_max** is the **maximum of ξ(τ) itself**, regardless of χ4. ξ keeps growing well after τ\*, up to τ ≈ 80–115 frames (~30 min), and only then decreases. Here ξ_max ≈ 2.6 µm, almost independent of ε (2.5–2.8 µm).

**Why they differ.** χ4 peaks earlier than ξ because after τ\* its amplitude Π(1−Π) falls quickly: most of the wall has moved and few persistent sites remain. But those that remain are correlated over longer and longer segments, and ξ measures the size of those segments, not how many there are. Thus ξ(τ\*) answers "at the moment of maximal heterogeneity, how large are the regions that move together?" (it mixes size and number of regions and depends on ε), while ξ_max answers "how large do the correlated regions become, and how long does it take?" (the most ε-robust cooperativity scale, but measured when Π is already small, ~5–10 %, with less statistics).

**In practice:** report ξ(τ\*) to compare with work based on the χ4 peak, and ξ_max as a characteristic scale independent of the threshold. In the strong-pinning tests (below) both behave the same: neither changes far from the defects.

### Main findings

1. **τ\* is resolved.** Unlike Study 1, χ4(τ) rises, peaks and decays. At ε = 3σ_Δh = 0.97 px, **τ\* = 7 frames ≈ 2.3 min**, with **ξ(τ\*) = 13.4 ± 0.1 (stat.) ± 0.3 (fit) px = 1.58 ± 0.04 µm**. Over ε = 0.5–5 px, ξ(τ\*) = 1.5–2.1 µm.
2. **τ\* is set by Π(τ\*) ≈ ½.** For every ε, τ\* is within 3 frames before and 1 frame after the lag at which half of the rows have moved by more than ε. So τ\* grows with ε (6 → 25 frames for ε = 0.5 → 5 px) and mostly measures the time to move a distance ε, as expected when χ4 ≈ Π(1−Π) × (correlated length). The peak χ4 itself (4.7–6.4) and ξ(τ\*) (13–17 px) depend more weakly on ε.
3. **The correlation length keeps growing well beyond τ\*.** ξ(τ) increases from ~8 px at τ = 1 to a maximum of **ξ_max ≈ 21–24 px = 2.5–2.8 µm at τ ≈ 80–115 frames (27–38 min)**, then decays as Π → 0. ξ_max changes by only ±1 px over a factor 10 in ε, so **ξ_max and its time scale of ~30 min are the most robust results** of this study, together with τ\*. The raw χ4 peaks much earlier only because its amplitude Π(1−Π) falls.
4. **χ4 and ξ are consistent at short τ.** C(n,τ) is exponential: (C−B)/A is a straight line on the semilog plot down to ~2 % of its amplitude. For τ ≲ 30 frames, the local normalized susceptibility (χ4 − B·L)/[Π(1−Π)] follows 2ξ, and the sum of the fitted model A·coth(1/2ξ) + B·L reproduces the direct χ4 within 3–21 % at τ\* (3 % at ε = 0.97 px, 21 % at ε = 5 px). For τ ≳ 50 frames the model sum falls below χ4: correlations at long lags are not a single exponential.
5. **Integer h** (`--int`, cross-check, no row-mean subtraction): ε = 1 px gives τ\* = 6 frames, ξ(τ\*) = 15.3 ± 0.3 px; ε = 2 px gives τ\* = 11 frames, ξ(τ\*) = 15.9 ± 0.3 px, consistent with the subpixel results.

### Row-mean subtraction

Without it, C(n) = ⟨p_i p_{i+n}⟩ − Π² is **negative at n ≈ 500–1000 px** (top half against bottom half of the wall). C(n) then contains a static term ⟨m_i m_{i+n}⟩ − Π², from rows that advance at different average rates over the 10 h (in the kymograph, the top overtakes the bottom late in the run; some rows are pinned by defects for a long time). This term sums to exactly zero over all pairs, so it is positive at short n and must be negative at large n. Subtracting m_i removes it, and the identity χ4 = Σ_n (1−|n|/L)·C(n) still holds exactly, because χ4 = (1/L) Σ_ij cov_t(p_i, p_j).

Effect at ε = 0.97 px:

| | without | with row-mean subtraction |
|---|---|---|
| C(n = 500–1000)·L at τ\* | −20.9 | −2.4 |
| ξ(τ\*): fit n ≤ 300 / n ≤ 1000 / B = 0 | 13.3 / 16.7 / 13.9 | 13.4 / 13.5 / 13.2 |
| ξ(τ = 90): fit n ≤ 300 / n ≤ 1000 / B = 0 | 22.5 / 22.7 / 21.0 | 21.9 / 19.7 / 19.7 |
| model sum / χ4 at τ\* | 1.7 | 0.97 |

So the fit-range ambiguity of ξ(τ\*) (±2 px before) disappears. At long lags (τ ≈ 90) a ~±1 px dependence remains, where Π ≈ 0.07 and fewer rows are persistent. The fitted plateau B·L is now small and slightly negative (−1 to −3), as expected from estimating m_i from a finite number of start times.

### Caveats

- **Resolution.** ξ(τ\*) ≈ 1.6 µm is close to the optical resolution of the 20x objective (~1 µm, NA 0.4, thesis) and spans ~14 px; ξ is resolved, but its smallest values (ξ ≈ 0.5–0.9 µm at τ = 1, depending on ε) are at the resolution limit.
- **ξ is measured along y**, not along the wall; the 1/cos θ correction (≈ 2 %) is not applied.

## Avalanches

`step4_avalanches.py`. An avalanche is a connected patch (8-neighbours) of pixels that switch within the same measurement window τ_m, taken from the arrival-time map; patches smaller than 20 px (0.28 µm², the resolution limit used in the thesis of this measurement) are discarded, and the static defects (whose arrival time is inpainted) are excluded. Figures: `fig_avalanches.png`, `fig_avalanche_map.png`; comparison with the thesis in its ROI: `fig_avalanches_thesis_roi.png`.

```bash
python step4_avalanches.py 0_8A_15s_20x_RT_1 --tau-m 1 3 10 30
python step4_avalanches.py 0_8A_15s_20x_RT_1 --tau-m 1 3 --roi 430 790 528 790 --tag _thesis_roi \
       --compare /media/…/ResumenMaestría_Mati/Data/Histogramas_lineal/Avalanchas_15s.txt
```

Arrival times are obtained by fitting a step to each pixel's normalized intensity s(t) (the k that maximises Σ_{t<k}(s_t − 0.5)). The first version counted the frames with s > 0.5; with only ~5 % contrast, isolated noisy frames shifted the count and split single avalanches over neighbouring frames.

| τ_m | avalanches ≥ 20 px | share of swept area (rest: lagunas) | median S [µm²] | max S [µm²] | median ℓ_y [µm] | τ (with cutoff) | S_cut [µm²] |
|---|---|---|---|---|---|---|---|
| 20 s (1 frame) | 6755 | 66 % | 0.54 | 4.0 | 2.2 | 0.50 ± 0.07 | 0.54 |
| 60 s | 5348 | 86 % | 0.79 | 8.0 | 2.6 | 0.56 ± 0.04 | 1.3 |
| 200 s | 2333 | 96 % | 1.71 | 23 | 3.6 | 0.55 ± 0.04 | 4.4 |
| 600 s | 684 | 99 % | 4.43 | 94 | 5.1 | 0.74 ± 0.04 | 26 |

τ and S_cut: maximum-likelihood fit of P(S) ∝ S^(−τ)·exp(−S/S_cut) above 20 px (the form used in the thesis); errors from the inverse Hessian.

### Findings

1. **At 20 s, a third of the advance is below the resolution.** Avalanches ≥ 0.28 µm² account for 66 % of the swept area; the rest switches in smaller steps (grey in `fig_avalanche_map.png`). With longer windows the steps coalesce, and at 600 s 99 % of the area is in resolved events.
2. **Avalanches are thin and elongated along the wall.** At 20 s the median lateral extent is ℓ_y ≈ 2.2 µm (19 px) and the typical width S/ℓ_y only ≈ 2 px (0.24 µm), below the optical resolution (~1 µm): the advance in one frame is a thin sliver of the wall, so S is effectively displacement × length. The lateral extent is of the same order as the correlation length of the persistence analysis (ξ ≈ 0.9 µm at τ = 20 s, 1.6 µm at τ\*, 2.6 µm at ~30 min).
3. **The size distribution is dominated by its cutoff.** At τ_m = 20 s, S_cut = 0.54 µm² is only twice the smallest resolved area, so there is no power-law range and τ is poorly determined. S_cut grows strongly with the window (0.54 → 1.3 → 4.4 → 26 µm² for 20 → 600 s): the measured events are mostly coalescences of smaller ones, as discussed in the thesis. The fitted τ (0.50–0.74) is well below the equilibrium-avalanche value 1.17 quoted in the thesis; a pure power law without cutoff would give 2.3 → 1.4, which only reflects the cutoff. Neither should be taken as an avalanche exponent.
4. **S vs ℓ_y.** The fitted 1+ζ goes from 0.6 (20 s) to 1.3 (600 s). At short windows the events are slivers a few px wide, so S ∝ ℓ_y^(<1) reflects the resolution, not the roughness; this is not a measurement of ζ.

### Lagunas

The thesis calls *lagunas* the parts of the swept area not covered by any detected avalanche: the wall went through them in steps smaller than the minimum area. They are the grey regions of `fig_avalanche_map.png` (defects in black). They cover 34 % of the swept area at τ_m = 20 s, 14 % at 60 s, 4 % at 200 s and 1 % at 600 s.

- **Structure.** They are not uniform: they form filaments aligned with the direction of motion, which persist over many frames and curve around the defects like flow lines. Around some defects there is also a grey ring: the optical halo of the defect, where the contrast is low too (within 10 px of a defect the laguna fraction is 45 %, against 33 % elsewhere). At 20 s, 28 % of laguna pixels are isolated single pixels, 30 % are in patches of 2–5 px and 43 % in patches of 6–19 px; 82 % touch an avalanche.
- **What they are.** Two possibilities were compared on the intensity time series of every pixel in a 230×210 px region without defects (x = 470–700, y = 130–340): (a) ambiguous borders between neighbouring avalanches, which would switch in two abrupt steps of ~0.5; (b) accumulation of steps below the resolution, which would switch gradually.

  | | avalanche pixels | laguna pixels |
  |---|---|---|
  | transition time (s from 0.8 to 0.2) | 18 frames (6 min) | 32 frames (11 min) |
  | largest one-frame drop of s near the arrival | 0.23 | 0.12 |

  (frame noise of s: 0.017). The drops in the lagunas are small but well above the noise and there are no 0.5 steps: **the lagunas are mostly places where the wall advanced by an accumulation of small steps, ≲ 1 px per frame, more slowly** (transition ~1.8× longer). The drop of 0.23 in the avalanche pixels corresponds to an advance of ~2 px over the ~8 px point-spread function (1 µm optical resolution): in one frame the wall advances ~2 px along a long segment. The difference between an avalanche and a laguna is the local advance per frame, not a qualitatively different motion.
- **Interpretation (to be tested).** The filaments look like *slow lanes*: positions along the wall that advance by fine creep while their neighbours jump. That they persist and follow the defects suggests they are set by the local disorder; this can be tested by comparing them with the per-row mean persistence (pinning map).

### Comparison with the thesis (same measurement)

The thesis list `Data/Histogramas_lineal/Avalanchas_15s.txt` has 1151 avalanches in a square ROI (x = 430–790, y = 528–790 px, camera coordinates), found by subtracting consecutive frames, Gaussian filtering and thresholding, without drift correction.

- **Same events.** 1042 of the 1151 (90 %) coincide with one of ours in the same frame (±1) within 15 px. In the ROI we find 1324 avalanches ≥ 20 px at 1 frame.
- **The thesis areas are about 2× larger** (median 81 px = 1.1 µm² vs 38 px = 0.52 µm²). For the largest thesis events (460–485 px), the change visible in the frame ratio covers ~105 px, which is what we measure; one of them (frame 1237) shows no visible change at all. The avalanches are 5–8 px wide, so filtering the difference image before thresholding widens them by a large factor.
- As a result, the thesis 1-frame distribution resembles ours at τ_m = 60 s (both S_cut ≈ 1.2 µm²; `fig_avalanches_thesis_roi.png`).
- **Time base.** The thesis labels its windows 15, 30, 45, 60 s (`Data/Exponente/exponentes_0_8A.txt`), i.e. 1–4 frames at 15 s per frame, while the acquisition interval was 20 s (see *Data*). Times in the thesis for this measurement are probably too short by a factor 3/4.

### Comparison with the thesis at 24.2 Oe (~52 °C)

The thesis also analysed 0_42A_20s_20x_52C_1 (24.2 Oe, ~52 °C); only its processed results are available (`Data/S_medio`, `Data/Histogramas_log`), not the frames. `compare_thesis_fields.py` puts them next to ours (`fig_thesis_fields.png`):

```bash
python 0_8A_15s_20x_RT_1/compare_thesis_fields.py /media/…/ResumenMaestría_Mati/Data
```

1. **Calibration on this measurement (thesis ROI).** Both methods count almost the same avalanches at every window (1151 vs 1324 at 1 frame, 476 vs 513 at 8, 149 vs 169 at 20), but the thesis mean area is larger by an amount that grows with the event size: +0.8 µm² at 1 frame (1.49 vs 0.68 µm²), +2.2 at 8 frames, +3.9 at 20 frames. This is what a halo of roughly constant width around each event does, as expected from filtering the difference image; for large windows (≳ 40 frames) the two agree.
2. **24.2 Oe vs 46.1 Oe, same (thesis) method.** At 24.2 Oe the mean area is smaller and grows more slowly with the window: the ratio to 46.1 Oe is 0.77 at 1 frame, 0.64 at 4, 0.58 at 8 and 0.44 at 16 frames. Fewer events coalesce per window, consistent with a smaller advance per frame at the lower field (the velocity of that measurement is not in the available data).
3. **Same shape.** Rescaled by their mean area, the two thesis distributions at 1 frame collapse onto each other: field and temperature change the scale, not the shape. Ours (46.1 Oe) has more weight at small S/S̄, consistent with the halo, which enlarges small events proportionally more.

Without the frames of the 24.2 Oe measurement, these conclusions rest on the thesis method; the halo correction obtained at 46.1 Oe cannot be transferred reliably because it depends on the event shapes. Running this pipeline on the raw frames of the other measurements (0_36A, 0_38A, 0_42A) is pending.

## Pinning map

`step5_pinning_map.py` (`fig_pinning_map.png`):

```bash
python step5_pinning_map.py 0_8A_15s_20x_RT_1          # eps = 0.97 px, tau = tau* = 7 frames
```

- **Local velocity map.** v(x,y) = 1/|∇t_arrival| (arrival map smoothed over 4 px) shows a network of **slow lanes**: lines about 1 µm wide (the resolution), elongated along the direction of motion, several hundred px long, many of them starting at defects like wakes. They are genuinely slow in the direction of motion, not lateral steps of the wall: the gradient of the arrival time is within 45° of the mean direction of motion in 68 % of the slowest pixels, against 73 % elsewhere.
- **Lagunas lie on the slow lanes.** The local velocity in laguna pixels is 0.53, 0.72 and 0.85 times that in avalanche pixels for smoothing over 2, 4 and 8 px: the difference is concentrated at the scale of the resolution. The laguna fraction per row has a correlation length of 10 px (1.2 µm) along the wall, i.e. the lanes are narrow.
- **Per-row persistence** (ε = 0.97 px, τ\* = 7 frames): m_i = 0.57 ± 0.05. It is correlated with the mean velocity of the row (r = −0.41: slower rows are more persistent) and its correlation length along the wall is 55 px = 6.5 µm, larger than ξ. The largest deviation is the band y ≈ 400–480 px behind the ring of defects, where the row velocity drops to 0.19 px/frame (mean 0.28) and m_i rises to 0.7. This large-scale heterogeneity is what `--row-mean` removes.
- **Lagunas are less persistent, not more.** Rows with more lagunas are slightly slower (r = −0.30) but less persistent (r = −0.20, rows without defects). Where the wall advances by fine creep, ≲ 1 px per frame, it accumulates more than ε within τ\*; where it advances by avalanches it stays still between jumps. The hypothesis of pinned lanes is not supported: the lanes are slow but steady.

**Picture.** The wall advances in two ways that coexist along its length: intermittent avalanches (jumps of ~2 px over segments of ~2 µm, with waiting in between) and slow, steady creep in steps below the resolution along narrow lanes, many of which trail from defects. The persistence and χ4 analysis is dominated by the first; the lagunas and the slow lanes are the second.

## Strong pinning vs collective creep

Question raised by colleagues: could the χ4 peak and ξ come from a few **strong pinning centres** rather than from collective creep (many weak pins, elastic correlation)? With strong pins, the correlated persistent regions should be centred on the pins, their size should be set by the pins and saturate early, and they should disappear far from the pins. With collective creep, correlated regions appear anywhere along the wall and ξ grows slowly with τ. `step6_pinning_tests.py` runs four tests (figures `fig_mask_defects_eps0.97.png`, `fig_pinning_tests.png`):

```bash
# near / far from the 33 static defects (mean spacing 124 px = 14.6 µm, ~9 ξ(τ*)); run at most a few at a time
python step2_persistence.py 0_8A_15s_20x_RT_1 --row-mean --eps 0.97 --tau-max 250 --mask-defects D [--near-defects]
python step6_pinning_tests.py 0_8A_15s_20x_RT_1 --eps 0.97
```

**1. Near vs far from the defects.** `--mask-defects D` keeps only the (row, t) pairs whose wall stays farther than D px from every defect during [t, t+τ]; `--near-defects` keeps the complement. D = 0 reproduces the unmasked results (|Δξ| < 1e-8).

| subset | data kept at τ\* | τ\* [frames] | χ4(τ\*) | ξ(τ\*) [µm] | ξ_max [µm] (τ) |
|---|---|---|---|---|---|
| all | 100 % | 7 | 5.2 | 1.58 ± 0.02 | 2.57 (88) |
| far, D = 10 px | 90 % | 7 | 4.7 | 1.61 ± 0.02 | 2.87 (95) |
| far, D = 20 px | 82 % | 7 | 4.8 | 1.59 ± 0.01 | 2.88 (95) |
| far, D = 40 px | 69 % | 9 | 5.1 | 1.68 ± 0.02 | 3.01 (105) |
| far, D = 80 px | 43 % | 11 | 6.2 | 1.59 ± 0.03 | 2.91 (92) |
| near, D = 20 px | 19 % | 21 | 5.0 | 1.71 ± 0.09 | 1.75 (28) |
| near, D = 40 px | 31 % | 7 | 5.1 | 1.54 ± 0.05 | 2.59 (70) |

**The χ4 peak and ξ survive intact far from the defects**: with the wall more than 80 px (9 µm, ~6 ξ_max) from any defect, the peak is still there with the same height and ξ(τ\*) is unchanged (1.59 µm); ξ_max is even slightly larger far from the defects (2.9–3.0 µm) than on the whole wall (2.6 µm). Near the defects the correlation is not stronger: ξ_max is smaller and the χ4 peak is broader and later (D = 20). Part of the reduction near the defects is geometric: the near subset is made of short segments of the wall (≈ defect size + 2D), which limits the measurable ξ; the far subset has no such limitation and is the conclusive one.

**2. Location of the persistent clusters** (runs of consecutive persistent rows at τ\*). Clusters longer than 2ξ(τ\*) (13 056 of 51 368) are not concentrated at the defects: the fraction of cluster centres within 10 px of a defect equals the surrogate (clusters moved at random along the wall at the same t, 200 replicas): ratio 1.00; within 20 and 40 px there is a small excess, ratio 1.10 and 1.04 (p = 0.005). Clusters planted at the defects, as a control, give ratios 7.4, 4.4 and 2.8, so the test is sensitive. Clusters are also *shorter* next to the defects (15.6 px within 10 px, ~20 px farther away): defects cut the correlated regions rather than create them.

**3. Growth of ξ(τ).** ξ grows from 7.8 to 21.9 px between τ = 1 and 88 frames without saturating. A saturating form a(1−e^{−τ/c})+d fits worst (AIC −105), a power law best (exponent 0.20, AIC −266), a log law in between (AIC −178). No early saturation at a pin-set scale.

**4. Hidden strong pins** (not visible as defects). A point-like strong pin holds the wall for a long time in a small spot. In the waiting map 1/v the slowest 1 % of the pixels hold 9 % of the total waiting time (slowest 5 %: 24 %): some concentration, but no dominant point-like spots. The slow regions are extended along the wall: correlation length of ln(1/v) is 11 px (1.3 µm) along the wall, of the order of ξ(τ\*), and 6 px across. (Locating "hidden pins" as the slowest points and testing whether persistent clusters sit there would be circular: persistence and slowness are the same thing.)

**Conclusion.** The visible strong defects do not produce the χ4 peak or the correlation length: both are unchanged, or ξ_max even larger, far from them, and persistent clusters are not attached to them. The slow regions are extended at the scale of ξ and ξ(τ) grows without saturation, as expected for collective creep. These tests cannot exclude pinning centres that are strong but dense and invisible (below the resolution), with spacing ≲ ξ; such a landscape is, in practice, the collective-pinning picture. The decisive remaining tests are (a) the field dependence of ξ and τ\* (collective creep: ξ grows as H decreases; pin-dominated: fixed), which needs the measurements at 20.7–24.2 Oe, and (b) the same analysis on simulations of an elastic line with weak disorder, with and without added sparse strong pins: step2–step6 run on any `h_xt_sub.npy` + `meta.json`.

## Roughness and structure factor

`step7_roughness.py` (`fig_roughness.png`, `roughness.npz`):

```bash
python step7_roughness.py 0_8A_15s_20x_RT_1 --every 5 --rmin 25
```

**Base plane.** The wall is tilted by θ0 = −11.8° on average, which would dominate any roughness measurement. Every 5th frame (360 frames), the domain is rotated by θ0 as a continuous image (linear interpolation, threshold 0.5) and the residual slope of each frame (≤ 3.9°) is removed with a line fit. On the base plane, 915 rows (108 µm) are inside the field in every frame. The rms width of the detrended wall is 5.5 µm.

**Overhangs.** On the base plane they are 4.7 % of the (row, frame) pairs, narrow along the wall (median 4 px, 90 % below 28 px) but deep (front − back: median 63 px, 90 % below 164 px): mostly tongues around the defects. Three heights are measured in every row: u_area (area-conserving column height, the reference), u_front and u_back. Their S(q) agree for q ≲ 1–2 µm⁻¹ and separate at larger q, where the jumps of u_front and u_back at the overhangs add a q⁻² tail. The fits therefore use ℓ = 2π/q from 2.9 to 27 µm (25–229 px), above both the optical resolution (~1 µm) and most of the overhang effect.

**Validation** (`--selftest`): synthetic walls with known ζ, tilted by 12°, rotated and analysed in the same way, give the same exponents as without tilt; S(q) recovers ζ without bias (0.48, 0.68 and 1.01 for 0.5, 0.66 and 1.0), whereas B(r) underestimates ζ ≥ 0.66 (0.58 and 0.78) and w(ℓ) underestimates ζ ≈ 1. S(q) is the reference estimator; its error is a block bootstrap over time (frames are strongly correlated).

| estimator | ζ from S(q), 2.9–27 µm | thirds of the run | ζ from B(r) | ζ from w(ℓ) |
|---|---|---|---|---|
| u_area | 0.82 ± 0.09 | 0.74, 0.53, 1.01 | 0.54 (far from defects 0.61) | 0.84 |
| u_front | 0.76 ± 0.05 | 0.75, 0.55, 0.85 | 0.52 | 0.72 |
| u_back | 0.69 ± 0.04 | 0.74, 0.53, 0.70 | 0.50 | 0.64 |

- **ζ ≈ 0.75 ± 0.1 over 2.9–27 µm** (but this range mixes two regimes; see *Local width with local rotation* below). The three heights bracket ζ_S between 0.69 and 0.82; with fit ranges from 1.9 or 5.9 µm the values stay within 0.68–0.84. B(r) gives 0.50–0.61, which by the self-test calibration corresponds to a true ζ ≈ 0.6–0.7. The roughness is compatible with the equilibrium exponent ζ_eq = 2/3 expected for creep at short scales, and clearly below the depinning value ζ_dep ≈ 1.25.
- **The exponent fluctuates in time** (0.53–1.01 between thirds of the run): the large scales of a single wall evolve slowly, so each third contains few independent configurations. The global value has a correspondingly large error.
- **Scale dependence.** ζ_eff(r) = ½ d ln B/d ln r decreases from ~0.9 at r ≈ 0.1–1 µm (below the resolution, not physical) to ~0.6 at 5–10 µm and then saturates. No clear crossover appears near ξ(τ\*) or ξ_max (1.6–2.6 µm), which are too close to the resolution to separate regimes.

### Local width with local rotation

Also from `step7_roughness.py` (`fig_local_width.png`, `local_width.npz`). Each frame's wall (1 in 5 frames) is cut into segments of projected length ℓ along the base plane (windows sliding by ℓ/2), and **each segment is rotated by its own tilt**: w² is the smallest eigenvalue of the covariance of its points, i.e. the variance normal to the segment's own best (orthogonal least-squares) line. ⟨w²(ℓ)⟩ is averaged over all segments and frames. The points are the full subpixel contour of the wall, resampled at uniform arc length, so overhangs are part of the shape and need no choice of height. Controls: the same on u_area; the line-removed width without rotation (w(ℓ) above); and the result without the segments whose own tilt differs from the base plane by more than 45° (`--max-tilt`). Self-test on synthetic walls: unbiased for ζ = 0.5 and 0.66 (0.53, 0.66), underestimates ζ = 1 (0.77).

**The wall is far from a small-slope line below ~30 µm.** The local tilt of the segments has a standard deviation of 36–38° for ℓ ≲ 10 µm, 29° at 10–30 µm and 12° at 30–108 µm; 16–21 % of the segments up to ~12 µm are tilted by more than 45° (steps, tongues around defects, flanks of bulges; panel 1 shows a 22 µm segment whose best line follows a horizontal tongue). In this regime the measured exponent depends on how the local tilt is handled, so the local rotation matters, as anticipated.

| range | rotation per segment (contour) | same, segments ≤ 45° | rotation, u_area | line removed, no rotation | sd of local tilt |
|---|---|---|---|---|---|
| 1–3 µm | 1.12 | 1.40 | 1.30 | 1.09 | 38° |
| 3–10 µm | 1.05 | 1.22 | 1.17 | 0.99 | 36° |
| 10–30 µm | 0.93 | 1.13 | 0.89 | 0.69 | 29° |
| 30–108 µm | 0.57 | 0.60 | 0.38 | 0.57 | 12° |
| 3–27 µm (default fit) | 0.98 ± 0.03 | 1.15 | 1.03 | 0.84 | 33° |

- **Two regimes.** At small scales (≲ 10 µm) the wall is very rough, ζ_eff ≈ 1.0–1.2 with every method; at large scales (≳ 30 µm), where the wall is locally flat, ζ_eff ≈ 0.5–0.6 and the methods agree (u_area is lower, 0.38). The crossover is at ~10–30 µm. The ζ ≈ 0.75 from S(q) over 2.9–27 µm averages over both regimes.
- **Possible interpretation (to be confirmed).** The theory of creep predicts several regimes: ζ_eq = 2/3 below the optimal (thermal-nucleus) length, ζ_dep ≈ 1.25 between that length and the size of the avalanches, and ζ_th = 1/2 at larger scales. A depinning-like roughness at a few µm and ζ ≈ 1/2 above ~30 µm would match the last two, with the crossover at the scale of the largest avalanches (their lateral extent reaches 10–55 µm in windows of 20–600 s). The resolution (~1 µm) hides the equilibrium regime if L_opt is sub-micron.
- **Caveats.** Below ~3 µm the resolution and the contour noise contribute; the large-scale regime spans only a factor 3.5 in ℓ and relatively few independent configurations of one wall; and with local slopes of 30–40° the description of the wall as a single-valued elastic line is itself only approximate at small scales.

### Height distribution

Also from `step7_roughness.py` (`fig_height_distribution.png`, `height_distribution.npz`): the distribution of the deviations δu(s,t) = u(s,t) − ⟨u⟩_s(t) of every point of the wall from its mean position, **on the base plane of each frame** (rotation by θ0 plus the residual line of the frame, equivalent to rotating each frame by its own θ(t)), normalised by the width σ(t) of each frame. u increases in the direction of motion, so δu < 0 means a part of the wall lagging behind. Errors: block bootstrap over time.

| | width σ | skewness | excess kurtosis |
|---|---|---|---|
| whole wall, base plane (u_area) | 5.5 µm | −0.38 ± 0.19 | −0.18 ± 0.32 |
| whole wall, mean only (no per-frame line) | 5.6 µm | −0.56 ± 0.23 | +0.23 ± 0.36 |
| windows of 20 µm, line removed | 1.7 µm | −0.28 ± 0.05 | +0.86 ± 0.26 |
| windows of 5 µm, line removed | 0.5 µm | +0.10 ± 0.12 | +7.3 ± 1.3 |

- **Whole wall:** not Gaussian, with a sharp cut at +2σ, a tail behind and a dip at the centre. At this scale the distribution is dominated by the few largest modes of a single wall (the large bulge), so it does not self-average; the skewness is barely significant. Without the per-frame line (mean only) the residual tilt adds to the asymmetry: the rotation matters here.
- **Windows of 20 µm:** close to Gaussian, with a small but significant **negative skewness** (−0.28 ± 0.05): more parts of the wall lag behind the local mean than run ahead of it, as expected for a wall held back locally by pinning while the rest advances.
- **Windows of 5 µm:** fat tails (excess kurtosis +7). They grow strongly with u_front and u_back (+16 and +26), so they come mostly from the overhangs and from steps at the scale of the resolution (σ = 0.5 µm ≈ 4 px here); they should not be read as a property of the elastic line.

## Next steps

- Map the per-row mean persistence m_i (and mean velocity per row) against the defect positions: the static heterogeneity removed by the row-mean subtraction is itself a measure of the pinning landscape.
- Avalanches: separate the coalescence effect, e.g. by following how S_cut(τ_m) and ℓ_y(τ_m) grow, and relate ℓ_y(τ_m) to ξ(τ) from the persistence analysis.
