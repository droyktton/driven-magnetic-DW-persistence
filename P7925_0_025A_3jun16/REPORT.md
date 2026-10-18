# Study 15 (exploratory): P7925, 0.025 A, 3 Jun 2016

**Status: limited extraction; the standard pipeline (steps 1–11) does not apply.** Only the shape of the wall and its mean velocity were measured, with a dedicated method. Read the caveats before using any number.

## Data

- `~/MOKEdata_local/3jun16` (copy of `MOBILEDATAMOKE/MOKE 2016/3jun16`): `P7925_0,025A_1,0s_N.tiff`, N = 1–3000, 1392×1040, camera offset ~33000.
- **Drive:** 0.025 A from the file names; `P7925__info.txt` says **I = 0.03 A**. Not resolved.
- Exposure 1.2 s; Δt = 6.25 s from the file times (median); 3000 frames = 5.2 h.
- **Scale 0.1175 µm/px assumed**, as in Studies 8, 9 and 11; the magnification of this run is not recorded.
- **Same sample as Studies 8, 9 and 11 (0.04 A)**, measured three weeks earlier. That makes it the low-field end of a P7925 field series: 0.02 A (27jun16, not analysed), 0.025 A (this study), 0.035 A (22jun16, not analysed), 0.04 A (Studies 8, 9, 11).

## What is in the field of view

- **One wall** crosses the whole field diagonally, from x ≈ 100 at the top to x ≈ 450 at the bottom.
  - The bright domain is on the left and grows to the right. Across the wall, log(left/right) = +0.023, i.e. **~2.3% domain contrast**.
  - In 5 h the wall sweeps only a **band ~10–15 px wide**, invisible in single raw frames.
- **A closed domain** of ~300 px sits on the right (x ≈ 900–1250, y ≈ 200–600) and also grows slightly. It is excluded here (`x < 760`) and not analysed.
- **Many point-like spots** (defects, or focus/halo changes) appear as dipoles in the ratio images and grow in contrast with time.
- **Drift:** dx −5.8 px, dy +0.5 px over the run (`drift_raw.npy`, from step1).

## Why step 1 fails

Run with `--offset auto --flatfield 2 --init-wall --bright-domain --crop-right 800`. The initial wall is found in only 89 of ~1000 rows (496–584). The ~2% domain contrast is below the gradient criterion, and the advance comes out NaN. Its outputs (`h_xt*.npy`, `qc/`, kymograph) are not committed. Only `meta.json` and `fig_drift_arrival.png` (the drift, reused below) are kept.

## Method used (`tiff_arrival_map.py`, new)

1. **Arrival map.** `python tiff_arrival_map.py ~/MOKEdata_local/3jun16 --outdir P7925_0_025A_3jun16/arrival --drift P7925_0_025A_3jun16/drift_raw.npy --x0 20 --x1 760 --y0 20 --y1 1020 --block 10`
   - Each frame is registered and divided by its median.
   - Frames are averaged in blocks of 10 (62.5 s).
   - Per pixel, the block of the largest step of log I is kept (`ta_all`), with its amplitude (`amp`).
   - The Otsu threshold is useless here: aperture edges and defects dominate it, and only 0.3% of the pixels pass. The swept band is found as follows.
2. **Band (wall) tracking** (`arrival/ridge_dp.npy`): a dynamic-programming ridge of the amplitude, smoothed σ = (3, 1.5) px, clipped at 0.08 to suppress defects, |Δx| ≤ 2 px per row. It follows the band over rows 60–990 (`arrival/fig_ridge.png`). Median amplitude profile across the ridge: FWHM ≈ 10 px.
3. **Wall shape** (`arrival/centre_line.npy`): subpixel centroid of the amplitude within ±12 px of the ridge, row by row. This is the **centre of the swept band**, i.e. the wall averaged over the run, not an instantaneous wall.
4. **Advance vs time** (`arrival/band_signal.npy`): mean log I over the band pixels (amplitude > 0.035 within ±15 px of the ridge), relative to strips 25–45 px on each side. Computed every 100 frames, each point a 10-frame average, directly from the frames. This avoids the per-pixel step estimator, which cannot place steps near the start or end of the run.

## Results

`arrival/fig_summary.png`: centre line, advance vs time, local width.

- **Total advance ≈ 10–15 px (1.2–1.8 µm)**, depending on the amplitude threshold: 15 px at 0.035, 10 px at 0.05, 19 px at 0.025.
- **Mean velocity ≈ 0.07–0.10 nm/s** (about 0.004 px/frame).
  - Study 8 (0.04 A): 0.30 nm/s. Study 11 (0.04 A, 23jun16): 0.15 nm/s.
  - The velocity is **not constant**: it roughly doubles between the first and second halves (0.13 → 0.27 of the total per hour). Possible causes: the transient after the field was switched on, or slow drifts of temperature or illumination.
- **Roughness of the centre line** (local width, linear detrending; quadratic above 100 px):
  - ζ ≈ 1.0–1.2 at l = 4–10 px (0.5–1.2 µm). This is the range affected by the smoothing (σ = 1.5 px) and the ridge constraint.
  - ζ ≈ 0.4 at 20–100 px (2.3–12 µm).
  - ζ(8–40 px) = 0.70, ζ(10–80 px) = 0.52.
  - Above ~150 px, w² grows with the macroscopic curvature of the wall.
  - w(50 px) ≈ 3.7 px (0.44 µm).

## Caveats

1. **The drive is uncertain:** 0.025 A (file names) vs 0.03 A (info file).
2. **The scale is assumed** (0.1175 µm/px). Lengths and v scale with it; exponents do not.
3. **The roughness comes from a different estimator than step7.** It is the centre of a ~10 px swept band, tracked with a constrained ridge and a centroid, so its ζ is not comparable with the step7 values of Studies 8, 9 and 11 (ζ ≈ 0.7–1.0) without a calibration of this estimator on the same synthetic walls. The low ζ ≈ 0.4 at 2–12 µm could come from the method. The ridge cannot follow slopes > 2 px/row, and the band centre averages the initial and final walls, which smooths features that changed during the run. **Do not use it as evidence for a field dependence of ζ.**
4. **There are no dynamics:** with a total advance of ~10–15 px, the persistence, avalanches and Λ steps cannot be run. The arrival times of single pixels are noisy, with step amplitudes ~1–7% against a background σ ≈ 1%.
5. **The velocity ratio with Studies 8 and 11** is only a factor 2–4, for a field ratio of 0.63–0.75. In a creep regime with U/kT ~ 30 a much larger factor would be expected. The three runs are on different days and probably at different positions, and Studies 8 and 11 themselves differ by a factor 2. **The v(H) comparison is not reliable from these runs.**
6. **The non-constant velocity** (×2 over the run) means v is an average over a non-stationary period.

## Also tried on 18 Oct 2026: 31ago16 (Pt/Co/Pt 0.44 A)

`PtCoPt_0_44A_31ago16/` holds a failed step1 run with the Study 3 options. Like 30ago16 (0.42 A), the wall barely moves and the domain contrast fades, and step1 loses the wall by frame ~600. `tiff_arrival_map.py` with the ridge tracking above may rescue its shape and mean velocity. Not done yet.
