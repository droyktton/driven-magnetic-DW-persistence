# Driven magnetic domain wall: persistence and dynamic correlation length

Analysis pipeline for movies or TIFF frame sequences of a driven magnetic domain wall, for example from polar MOKE microscopy. The pipeline:

1. Extracts the wall position h(x,t).
2. Computes the persistence field between times t and t+τ.
3. Measures the dynamic correlation length ξ(τ) and the four-point susceptibility χ4(τ).

## Studies

Each data set is a separate study with its own folder and report:

| Study | Data | Report |
|---|---|---|
| 1 | `magnetic_fliped/`: video, 78 frames at 25 fps (Δt = 40 ms), 0.17 µm/px, wall moving upward | [`magnetic_fliped/REPORT.md`](magnetic_fliped/REPORT.md) |
| 2 | `0_8A_15s_20x_RT_1/`: 1800 Micro-Manager TIFF frames, Δt = 20 s (10 h), 20x, 0.1175 µm/px, wall moving to the right | [`0_8A_15s_20x_RT_1/REPORT.md`](0_8A_15s_20x_RT_1/REPORT.md) |
| 3 | `PtCoPt_0_48A_1sep16/`: Pt/Co/Pt, drive 0.48 A, 1800 TIFF frames, Δt = 22.05 s (11 h), 20x, 0.1175 µm/px, wall moving to the right | [`PtCoPt_0_48A_1sep16/REPORT.md`](PtCoPt_0_48A_1sep16/REPORT.md) |
| 4 | `PtCoPt_0_50A_2sep16/`: same sample and protocol as Study 3, drive 0.50 A, 1799 frames, Δt = 22.05 s (11 h) | [`PtCoPt_0_50A_2sep16/REPORT.md`](PtCoPt_0_50A_2sep16/REPORT.md) |
| 5 | `0_36A_20s_20x_51C_2_1/`: same sample and field of view as Study 2 (thesis of M. Grassi), 0.36 A = 20.7 Oe at 51 °C, 1800 frames, Δt = 20 s | [`0_36A_20s_20x_51C_2_1/REPORT.md`](0_36A_20s_20x_51C_2_1/REPORT.md) |
| 6 | `0_38A_20s_20x_51C_1/`: as Study 5, 0.38 A = 21.9 Oe at 51 °C | [`0_38A_20s_20x_51C_1/REPORT.md`](0_38A_20s_20x_51C_1/REPORT.md) |
| 7 | `0_42A_20s_20x_52C_1/`: as Study 5, 0.42 A = 24.2 Oe at 52 °C | [`0_42A_20s_20x_52C_1/REPORT.md`](0_42A_20s_20x_52C_1/REPORT.md) |
| 8 | `P7925_0_04A_24jun16/`: sample P7925, 0.04 A, 6000 TIFF frames, Δt ≈ 9.5 s (15.9 h), scale 0.1175 µm/px assumed, horizontal wall (`--rotate 3 --bright-domain`) | [`P7925_0_04A_24jun16/REPORT.md`](P7925_0_04A_24jun16/REPORT.md) |
| 9 | `P7925_0_04A_26jun16/`: same sample and drive, two days later, 6000 frames, Δt ≈ 7.8 s (13 h) | [`P7925_0_04A_26jun16/REPORT.md`](P7925_0_04A_26jun16/REPORT.md) |
| 10 | `GdFeCo_B1743_600mA_16may17/`: ferrimagnetic GdFeCo (B1743, 10 nm), 600 mA, 10800 Micro-Manager TIFF frames, Δt = 2 s (6 h), scale 0.1175 µm/px assumed, wall moving to the right (`--bright-domain --crop-right 1200`) | [`GdFeCo_B1743_600mA_16may17/REPORT.md`](GdFeCo_B1743_600mA_16may17/REPORT.md) |
| 11 | `P7925_0_04A_23jun16/`: sample P7925, 0.04 A (as Studies 8–9, one day earlier than 8), 3000 TIFF frames, Δt = 8.25 s (6.9 h), scale assumed; manual swept threshold and rows 130–979 (`--swept-thr 0.018 --rows 130 980`) | [`P7925_0_04A_23jun16/REPORT.md`](P7925_0_04A_23jun16/REPORT.md) |
| 12 | `CoNi_bubble_130mA_10mar16/`: Pt/(Co0.2Ni0.6)×3/Ta, one expanding bubble, 10 ms field pulses (130 mA, small coil), 99 video frames (one per pulse), 5X, 0.467 µm/px assumed; a planar defect cuts a notch (sector excluded) | [`CoNi_bubble_130mA_10mar16/REPORT.md`](CoNi_bubble_130mA_10mar16/REPORT.md) |
| 13 | `PtCoNiAl_150Oe_8mar16/`: PtCoNiAl, 150 Oe, many bubbles nucleating and merging, 35 video frames, Δt = 3 ms, 5X (0.467 µm/px); one isolated bubble analysed | [`PtCoNiAl_150Oe_8mar16/REPORT.md`](PtCoNiAl_150Oe_8mar16/REPORT.md) |
| 14 | `PtCoNiAl_10_4Oe_13mar16/`: PtCoNiAl, 10.4 Oe, two creep walls crossing the field, 250 video frames, Δt = 20 s, 20X (0.117 µm/px) | [`PtCoNiAl_10_4Oe_13mar16/REPORT.md`](PtCoNiAl_10_4Oe_13mar16/REPORT.md) |
| 15 | `P7925_0_025A_3jun16/` (exploratory): sample P7925, 0.025 A (0.03 A in the info file), 3000 TIFF frames, Δt = 6.25 s (5.2 h), scale assumed; ~2 % contrast and ~10–15 px total advance, so step1 fails: wall shape and mean velocity only, from `tiff_arrival_map.py` and ridge tracking; see its caveats | [`P7925_0_025A_3jun16/REPORT.md`](P7925_0_025A_3jun16/REPORT.md) |
| 16 | `GdFeCo_B1743_550mA_15may17/`: same sample as Study 10, one day earlier, 550 mA, 10000 TIFF frames, Δt = 2 s (5.6 h), scale assumed; the wall rotates (tilt 26° → 43°), so step 7 runs with `--from-h` | [`GdFeCo_B1743_550mA_15may17/REPORT.md`](GdFeCo_B1743_550mA_15may17/REPORT.md) |
| 17 | `GdFeCo_B1743_500mA_16feb17/`: same sample, 16 Feb 2017, 500 mA, 10800 TIFF frames, Δt = 2 s (6 h), exposure 20 ms, rotated 180° (`--rotate 2`), 108 µm wall; non-stationary motion, step 2 also on the moving window; also v(I) from the pulse movies of the same day (`pulse_velocity.py`) | [`GdFeCo_B1743_500mA_16feb17/REPORT.md`](GdFeCo_B1743_500mA_16feb17/REPORT.md) |

Cross-study question, qEW or qKPZ below the crossover length l₀: [`UNIVERSALITY.md`](UNIVERSALITY.md). Studies 12–14 (bubbles, walls crossing a circular aperture) use the arrival-map tools `arrival_map.py` and `contour_analysis.py` instead of steps 1–11.

Requires numpy, scipy, matplotlib and scikit-image; for videos also imageio-ffmpeg (or an `ffmpeg` binary on the PATH); for TIFF folders also tifffile and Pillow, and ffmpeg for the QC movie; `arrival_map.py` needs OpenCV (cv2).

## Usage

Each movie gets its own output folder. By default the folder is named after the video file without its extension. For a TIFF folder it is named after the folder, or after its parent when the folder is a Micro-Manager `PosN`.

```bash
python step1_extract_wall.py MOVIE.mp4 --um-per-px 0.17        # -> MOVIE/ with h(x,t) and meta.json
python step1_extract_wall.py /path/RUN/Pos0 --um-per-px 0      # TIFF frames -> RUN/
python step2_persistence.py MOVIE --eps 0.25 0.5 0.75 1 1.5 2 3  # subpixel h
python step2_persistence.py MOVIE --int                          # integer h, ε = 1, 2
python step3_eps_sweep.py MOVIE                                  # sensitivity to ε
python step4_avalanches.py MOVIE --tau-m 1 3 10                  # avalanches (TIFF studies)
python step5_pinning_map.py MOVIE                                # pinning map (TIFF studies)
python step7_roughness.py MOVIE --rmin 25                        # roughness on the base plane (TIFF)
```

Step 1 writes the frame rate and the spatial scale to `MOVIE/meta.json`, and steps 2 and 3 read them from there. By default, the frame rate is the one stored in the video (for TIFF: from the timestamps in Micro-Manager's `metadata.txt`) and the scale is 0.17 µm/px. You can override either with `--fps` or `--um-per-px` in any step. `--um-per-px 0` reports everything in pixels only.

The original video of the example is not in the repository. Without it you can start from step 2, because `magnetic_fliped/h_xt*.npy` and `meta.json` are included. To rerun step 1, copy `magnetic_fliped.mp4` into the repository root.

## Coordinates: base plane, x and h

The correlation function C(n,τ), and therefore ξ, is defined **along the wall**. Before computing it, a coordinate x along the wall must be defined, and for that the **base plane** of the wall has to be found: the straight line (a plane, for a 2D interface) that fits the wall on average, i.e. its mean orientation. Then

- **x** is the position along the base plane, sampled at equally spaced points x_i, and n = |i − j| is a distance along it;
- **h(x,t)** is the displacement of the wall normal to the base plane;
- p_i, C(n,τ), ξ and χ4 are all computed in these coordinates, and the same base plane has to be used at t and t+τ, so that a column i follows the same position along the wall.

If x is taken along an image axis that is not the base plane, the geometry is mixed. For a wall tilted by θ with respect to that axis, a distance n along the image axis corresponds to n/cos θ along the wall, and a displacement Δ along the other image axis to a normal displacement Δ·cos θ. With a large tilt, a fixed image column does not even stay on the same part of the wall as it advances.

**What the scripts do.** The base plane is taken as an image axis:

- **Video input:** the horizontal axis of the image after `--rotate`/`--flipud`/`--fliplr`. Rotate the movie so that the wall is roughly horizontal; a residual tilt is not corrected.
- **TIFF input:** the vertical axis of the image (rows y), with h along the image x axis. Step 1 fits the actual base line in every frame, x = a + y·tan θ(t), and reports θ(t) in `tiff_extra.npz` (`theta_deg`) and `fig_tilt_overhang.png`.

So, for a wall tilted by θ, ξ measured along the image axis has to be multiplied by 1/cos θ to be a length along the wall, and ε corresponds to a normal threshold ε·cos θ. In Study 2, θ = −10° to −14°, so the corrections are 1.5–3 % and are not applied. For a strongly tilted or curved wall, the frames should be rotated onto the base plane (or h measured normal to it) before step 2; this is not implemented.

**Effect of the tilt on what the pipeline computes.** In Study 2 the wall is tilted on average by θ ≈ −12° (−14° to −8.5° over the run), and h(y,t) contains that slope (~200 px between the top and the bottom of the field). This is harmless for the persistence analysis:

- the persistence compares h in the *same* row at t and t+τ, so a fixed tilt cancels exactly in Δh. The tilt changes by ~3° in 10 h, which moves the ends of the wall by ~0.01 px within τ\* = 7 frames and ~0.1 px within τ ≈ 90 frames, against ε ≈ 1 px;
- C(n,τ) and χ4 are computed on the binary field p_i, not on h, so the slope of h does not enter;
- what remains are the geometric factors above (1/cos θ for lengths along the wall, cos θ for the normal threshold), 1.5–3 % here.

The tilt does matter for the **roughness** of the wall (the correlation function or structure factor of h itself): there the slope would dominate and give a spurious roughness exponent ζ ≈ 1. `step7_roughness.py` therefore rotates the wall onto its base plane and subtracts the residual slope of each frame before measuring the roughness.

## Input requirements

The method assumes the following about the movie. Items 1, 2 and 4 refer to the video detection; the TIFF detection has its own assumptions, listed in its section below. Step 1 checks the ones marked ✔ and prints warnings, which are also stored in `meta.json` and repeated by step 2.

1. **Orientation.** The bright domain is on top, the dark domain at the bottom, and the wall moves upward, so h decreases with t. For other geometries, use `--rotate K` (K × 90° counter-clockwise), `--flipud`, `--fliplr` or `--invert` (swap bright and dark). ✔ Step 1 warns if the mean wall position does not move upward.
2. **A single wall that spans the full width and is single-valued in x.** There must be no overhangs and no nucleated bubbles or dark spots in the bright domain ahead of the wall. h is defined as the first dark row from the top, so any dark feature above the wall is taken as the wall. Use `--roi X0 X1 Y0 Y1` to crop out edges, scale bars or overlaid text.
3. **The wall must stay inside the frame.** If the wall touches the top edge, or has not yet entered from the bottom, h stays fixed. Those columns would then look persistent and inflate Π and χ4. ✔ Step 1 keeps only frames where the wall is more than `--edge-margin` px (default 2) from the top and bottom edges in every column. It uses the longest contiguous run of such frames and reports which frames it dropped.
4. **Bimodal contrast that is stable in time.** The Otsu threshold is global over the whole stack, so illumination drift or bleaching during the movie shifts the detected wall.
5. **Monotonic advance.** The automatic noise estimate σ_Δh treats backward displacements as pure noise. ✔ Step 1 reports the fraction of one-frame retreats larger than `--retreat-px` (default 1 px) and warns if it exceeds 1 %. In that case, pass `--eps` explicitly to step 2.
6. **Intermittent motion.** Estimating the noise needs columns that stay still for some frames. ✔ Step 2 warns if there are fewer than 500 backward Δh samples.
7. **Enough frames.** τ runs up to n_frames/2, but Π(τ) vanishes once the wall has moved more than ε almost everywhere. The example has 78 frames, which is already on the low side.

## Scripts

### `step1_extract_wall.py`: wall detection

1. Extracts all frames with `ffmpeg -vsync 0` into `MOVIE/frames/`. This folder is not tracked by git.
2. Loads the frames as grayscale and applies `--roi`, `--rotate`, `--flipud`, `--fliplr` and `--invert`, in that order. The ROI is given in the coordinates of the original video.
3. Computes a **global Otsu threshold** over the whole stack, so that the threshold is the same for every frame.
4. Smooths each column vertically with `gaussian_filter1d`, σ = `--sigma-y` (default 2 px).
5. Defines **h(x,t)** as the first row from the top where the smoothed intensity drops below the threshold, that is, the bright-to-dark transition. It also computes a **subpixel version** by linear interpolation of the threshold crossing between rows r−1 and r.
6. Applies a median filter along x (kernel `--median-x`, default 5) to both versions, to remove single-column glitches.
7. Runs the edge and monotonicity checks described above.

Outputs in `MOVIE/`:
- `h_xt.npy` (integer) and `h_xt_sub.npy` (subpixel). Both have shape (frames used, width) and are in px.
- `meta.json`: fps, µm/px, frames used, threshold, transformations and warnings.
- `fig_qc_overlay.png` and `fig_h_mean.png`.

#### TIFF folder input

When the argument is a folder, step 1 reads its `*.tif` frames in numeric order and uses a different detection, built for long, slow acquisitions with weak contrast, uneven illumination, static defects and sample drift. It assumes a **dark domain on the left growing to the right**; the wall may be tilted and rough.

1. **Drift.** Each frame is registered to frame 0 by phase correlation of the high-passed image (defects give the texture). The drift curve is smoothed in time and applied with subpixel shifts. Rows that the vertical drift leaves without data are trimmed, plus `--trim-y` px (default 12) at each end.
2. **Right crop.** Columns from where the illumination of the field of view falls below 90 % of its plateau (minus 20 px) are discarded, in camera coordinates and the same for every frame. `--crop-right N` sets the column by hand.
3. **Per-pixel references.** Every frame is divided by the median of the pixels that never switch (global brightness drift). The bright reference B is the mean of the first `--nref` frames and the dark reference D the mean of the last ones. The swept area is where (B − D)/B exceeds its Otsu threshold. The normalised intensity s = (I − D)/(B − D) is 1 before and 0 after the wall passes, which cancels vignetting and static defects.
4. **Arrival-time map.** For each pixel, the switching frame is the best fit of a step (s = 1 before, 0 after) to its normalised intensity s(t), i.e. the k that maximises Σ_{t<k}(s_t − 0.5); this assumes a monotonic advance and is robust to isolated noisy frames. Defects inside the swept area show no contrast; they get the arrival time of the nearest valid pixel.
5. **Wall.** The domain at frame t is the set of pixels with arrival ≤ t connected to the left edge, with holes filled, so isolated spots ahead of the wall are ignored. The drawn wall (overlays) is the right-most domain pixel of each row y, i.e. the front where the wall folds back; the subpixel position comes from the 0.5 crossing of s. The position used for the analysis is x_eff(y,t), the number of switched pixels in row y: it equals the front where the wall is single-valued, and changes only by the area that really switched where the wall folds around defects (the front can jump when a thin unswitched channel is pinched off). Then h = W − 1 − x_eff, so that h decreases with t as in the video convention, and the columns of h are the image rows y.

**Options for aperture images with a camera offset** (Studies 3 and 4; all off by default, which leaves the steps above unchanged):

- `--offset auto` subtracts the camera offset, taken as the median of the four 40×40 px corners of each frame (dark, outside the illuminated aperture), or a fixed number. Without it, ratios such as (B − D)/B are computed on a signal that sits on a large constant offset. With an offset, only pixels inside the illuminated aperture (eroded by 30 px, holes filled) are used for the thresholds and the normalisation; the aperture is fixed on the camera, so its edge moves in the registered frames.
- `--flatfield N` corrects an illumination profile that changes during the run: in each frame, a 2D polynomial of order N in log(I/F), with F the mean of the first frames, is fitted to the never-switching pixels (robust refit, residuals beyond 3 MAD excluded) and divided out. If the never-switching region is too small, the whole aperture is used with the robust fit.
- `--crop-left N` discards the columns x < N of the registered frame (e.g. a vignetted aperture edge).
- `--init-wall` locates the initial wall directly as the dark-to-bright step: a dynamic-programming path through the rows that maximises the gradient d/dx log I of the first frames (|Δx| ≤ 3 px between rows, gradient clipped so that short defects do not capture it). Everything left of it is the initial domain, and every swept patch that touches it counts as wall motion (hysteresis threshold for the swept area). Needed where parts of the wall never move during the run: there the swept area has no pixels, and the method above cannot place the wall. Only the longest run of rows where the path is inside the aperture, on a real step, and with `--init-room` px (default 150) of aperture to its right is kept.
- `--init-wall` also fixes the initial domain: pixels left of the path (minus 3 px) are never counted as switched, even where the hysteresis threshold would leak into them (Study 9).
- `--rotate K` also applies to TIFF input (K × 90° counter-clockwise, after the offset subtraction), for walls that are horizontal in the camera (Studies 8 and 9). The drift saved for `--reuse-drift` is in the rotated frame.
- `--bright-domain`: the growing domain is the bright one (the sign of the switching contrast is reversed).
- `--rows Y0 Y1` restricts the output to the image rows Y0 ≤ y < Y1 of the registered frame (the y axis of the QC figures). The segmentation itself still uses all rows; use it to drop parts of the wall that run along defects or the aperture edge.
- `--swept-thr X` sets the threshold of the switching contrast (B − D)/B for the swept area by hand. By default it is Otsu's threshold; when Otsu marks more than half of the field or less than 1 % of it (with a large drift, the halos of the moved defects can dominate the histogram), it is recomputed on the contrast clipped to its 1–99.5 percentiles. Use `--swept-thr` when the swept band is fainter than the defect halos even then (Study 5).
- `--reuse-drift` reuses the frame registration saved by an earlier run in `DIR/drift_raw.npy` (the slowest part), to iterate on the other options.

h is measured along x in each row y, i.e. with the image y axis as base plane; see *Coordinates: base plane, x and h* above for the tilt correction.

Extra outputs: `tiff_extra.npz` (drift, arrival map, swept area, defects, x(y,t), tilt, overhang rows per frame), `fig_drift_arrival.png`, `fig_kymograph.png`, `fig_tilt_overhang.png`, and in `qc/` an overlay movie of every frame (`overlay.mp4`) plus `overlay_NNNN.png` every `--qc-every` frames.

### `step2_persistence.py`: persistence, C(n,τ), ξ(τ), χ4(τ)

Reads `MOVIE/h_xt_sub.npy`, or `MOVIE/h_xt.npy` with `--int`. For each ε (`--eps`) and each τ = 1…n_frames/2 it computes the quantities below. All averages run over columns i and over all start times t.

- p_i(t,τ) = 1 if |h(x_i,t+τ) − h(x_i,t)| < ε (the column did not move), and 0 otherwise.
- Π(τ) = ⟨p_i⟩, the fraction of persistent columns.
- C(n,τ) = ⟨p_i p_{i+n}⟩ − Π², computed with an FFT for every n.
- **ξ(τ)** comes from fitting C(n,τ) = A·e^(−n/ξ) + B. The constant B is needed because the global Π² is subtracted, so fluctuations of Π between different t leave a plateau at large n.
- **χ4(τ)** = L·Var_t[Π(t,τ)].
- **Normalized χ4**, χ4/[Π(1−Π)], in px. χ4 is roughly amplitude × correlation length, and the amplitude C(0) = Π(1−Π) falls quickly with τ as Π → 0. Dividing it out leaves an effective correlated length, which is ≈ 2ξ when the plateau term B·L is small. The script also reports the local part, (χ4 − B·L)/[Π(1−Π)], which removes the global frame-to-frame fluctuations of Π.
- As a consistency check, it compares three quantities:
  - χ4, computed directly.
  - Σ_n (1−|n|/L)·C(n). This is an exact identity, so it always matches.
  - The sum of the fitted model, A·coth(1/2ξ) + B·L. This is the non-trivial check.
- τ\* is the τ at which χ4 is maximal.

ξ is fitted only when there are at least 200 persistent events (`--min-pers`); otherwise it is left as NaN.

**`--row-mean`.** Computes C(n,τ) = ⟨(p_i − m_i)(p_{i+n} − m_{i+n})⟩, with m_i = ⟨p_i⟩_t the mean persistence of column i. This removes the static term ⟨m_i m_{i+n}⟩ − Π², which appears when columns advance at different average rates (for example columns pinned by defects for a long time). That term sums to zero over all pairs, so it is positive at short n and negative at large n, and it biases the fit of ξ. The identity χ4 = Σ_n (1−|n|/L)·C(n) still holds exactly. The output files get the tag `rm_` (e.g. `persistence_sub_rm_eps1.npz`), and `step3_eps_sweep.py --row-mean` reads them. Recommended for long acquisitions; see Study 2.

**`--mask-defects D` (TIFF studies).** Uses only the (row, t) pairs whose wall stays farther than D px from every static defect during [t, t+τ]; `--near-defects` keeps the complement. Π, χ4 = L_eff·Var_t Π(t) and C(n) are computed with these weights (tags `md<D>_` and `nd<D>_`). D = 0 reproduces the unmasked results. **`--tau-max T`** limits the lags (the masked computation is ~3× slower; with 1800 frames, τ ≤ 250 is enough for τ\* and ξ_max). Each run is single-threaded and takes several minutes; on a desktop, run at most ~3–4 at a time with `nice -n 19`.

**Noise estimate.** The script reports two estimates:
- σ_noise(x): high-frequency roughness along x.
- **σ_Δh**: the width of the peak of stationary columns in Δh(τ = 1). It is measured from the backward Δh values, which are pure noise when the wall never retreats.

Without `--eps`, the script uses ε = 1, 2 and 3σ_Δh. With `--int` it uses ε = 1 and 2, because integer h cannot resolve noise below 1 px.

When a scale is available, the summary also gives ξ and the mean wall velocity in µm and µm/s, and the ξ panels get a right-hand µm axis.

Outputs in `MOVIE/`, one set per ε:
- `persistence_<tag>eps<ε>.npz`: all curves, plus σ_Δh, fps and µm/px.
- `fig_Cn_tau_<tag>eps<ε>.png`
- `fig_xi_chi4_<tag>eps<ε>.png`

The tag is `sub_` for subpixel h and empty for integer h. The script also prints a table per τ and a summary.

### `step3_eps_sweep.py`: sensitivity to ε

Reads every `MOVIE/persistence_sub_eps*.npz` file (the integer-h files with `--int`, the row-mean ones with `--row-mean`). It overlays ξ(τ) and χ4(τ) for each ε, with τ\* (the maximum of χ4) marked, and plots against ε: ξ(τ\*), ξ_max and ξ(τ=1), with χ4(τ\*) on a second axis; and the time scales τ\*, the lag of ξ_max and the lag at which Π = ½. It also prints these values as a table. σ_Δh, fps and µm/px are read from the `.npz` files. Output: `MOVIE/fig_eps_sweep_sub.png` (`_int`, `_sub_rm`, `_rm` for the other variants).

### `step4_avalanches.py`: avalanche statistics (TIFF studies)

Reads the arrival-time map in `DIR/tiff_extra.npz`. An avalanche is a connected patch (8-neighbours) of pixels that switch within the same measurement window τ_m (`--tau-m`, in frames; several values can be given). Patches smaller than `--smin` px (default 20) are discarded. For each avalanche it stores the area S, the extents ℓ_y (along the wall) and ℓ_x, the centroid and the frame, and flags those touching the image edge or the never-reached region.

- P(S) is fitted by maximum likelihood with P(S) ∝ S^(−τ)·exp(−S/S_cut) above `--smin`, with errors from the inverse Hessian. A pure power law (no cutoff) is also printed, for reference only.
- S ∝ ℓ_y^(1+ζ) is fitted on the median S in log bins of ℓ_y.
- `--roi X0 X1 Y0 Y1` keeps avalanches whose centroid is in a region (camera px); `--compare FILE` overlays another list of areas in px (second column), e.g. the thesis list; `--tag` adds a suffix to the output files.

Outputs: `avalanches<tag>_tm<m>.npz`, `fig_avalanches<tag>.png` (P(S) with fits, P(ℓ_y), S vs ℓ_y) and `fig_avalanche_map<tag>.png` (avalanches of the first window coloured by frame; grey: area that switched in steps smaller than `--smin`).

### `step5_pinning_map.py`: pinning map (TIFF studies)

Local wall velocity v(x,y) = 1/|∇t_arrival| from the arrival-time map smoothed over `--sigma` px (default 4), lagunas (swept pixels that switched in 1-frame patches smaller than `--smin` px), and per-row values: mean persistence m_i at threshold `--eps` and lag `--tau` (default: τ\* for that ε), mean velocity and laguna fraction of each row. Prints the velocity in lagunas vs avalanches, the correlations between the per-row quantities (all rows and rows without defects in the path) and their correlation lengths along the wall. Output: `fig_pinning_map.png`.

### `make_overlay_gif.py`: shareable movie (TIFF studies)

`python make_overlay_gif.py DIR TIFF_DIR [--every 5] [--scale 0.5] [--fps-gif 12]` writes `DIR/overlay_movie.gif`: every N-th frame, drift-corrected and cropped, scaled, with the detected wall in red, the elapsed time and a 20 µm scale bar (when the scale is known). A fixed palette keeps the GIF small while preserving the overlay colours.

### `step6_pinning_tests.py`: strong pinning vs collective creep (TIFF studies)

`python step6_pinning_tests.py DIR [--eps 0.97] [--tau T] [--nsurr 200]`. (1) Overlays ξ(τ), χ4(τ) and χ4/[Π(1−Π)] for the whole wall and for the step2 `--mask-defects` (far) and `--near-defects` (near) results, and tabulates τ\*, ξ(τ\*) and ξ_max (`fig_mask_defects_eps<ε>.png`). (2) Tests whether persistent clusters longer than 2ξ(τ\*) sit closer to the defects than clusters moved at random along the wall (surrogate), with a planted-cluster control. (3) Fits ξ(τ) with power, log and saturating laws (AIC). (4) Measures how concentrated the waiting map 1/v is and its correlation lengths along and across the wall (`fig_pinning_tests.png`).

Simulated walls can be analysed with step2/step3 by writing `DIR/h_xt_sub.npy` (frames × columns, h decreasing as the wall advances) and a `DIR/meta.json` with at least `fps` and `um_per_px` (0 for none).

### `step7_roughness.py`: roughness and structure factor

`python step7_roughness.py DIR [--every 5] [--rmin 16] [--mask-d 20]` (TIFF studies), `--from-h` for video studies or simulations (uses `h_xt_sub.npy`, no rotation), `--selftest` for a synthetic check. The wall is put on its **base plane**: the domain of every N-th frame is rotated by the mean tilt θ0 (as a continuous image, linear interpolation, threshold 0.5) and the residual slope of each frame is removed. Along each row of the rotated frame it measures three heights: **u_area** (area-conserving column height, single-valued, the reference), **u_front** (last crossing) and **u_back** (first crossing); where they differ there is an overhang, and the exponent is trusted only in the range of scales where their S(q) agree. Outputs: S(q) (Hann window, averaged over frames and over thirds of the run), B(r) = ⟨[u(s+r)−u(s)]²⟩ with its local slope ζ_eff(r), the local width w(ℓ), B(r) without pairs near defects, ζ from each (fit range `--rmin` < ℓ < L/4; S(q) log-binned, error by block bootstrap over time), overhang statistics, `roughness.npz` and `fig_roughness.png`. The self-test shows that S(q) is unbiased while B(r) and w(ℓ) underestimate large ζ; use S(q) as the reference. It also computes the **local width with local rotation**: the wall of each frame is cut into segments of projected length ℓ along the base plane, each segment is rotated by its own tilt (w² = smallest eigenvalue of the covariance of its points, i.e. the variance normal to its orthogonal least-squares line), and ⟨w²(ℓ)⟩ is averaged over all segments and frames; for TIFF studies the points are the full subpixel contour (overhangs included), with u_area, the unrotated line-removed width and `--max-tilt` (default 45°: segments more tilted than this from the base plane are excluded in an extra curve) and `--max-fold` (default 1.5: segments whose contour inside the window is longer than this factor times the median at the same ℓ, i.e. folded, are excluded in another curve) as controls, plus the spread of the local tilt vs ℓ (`local_width.npz`, `fig_local_width.png`). It also computes the **local structure factor S(q, ℓ)** on the same rotated segments (4·rmin ≤ ℓ ≤ 0.8 L): the points are projected on the segment's own axis and binned at 1 px (mean normal displacement per bin), the line through the end points is removed (the segments are not periodic; end matching is unbiased, whereas a Hann window overestimates ζ in short segments), S(q, ℓ) = |FFT|²/N is averaged in log bins of q, and ζ(ℓ) is fitted over 2π/ℓ < q < 2π/rmin and over the low band 2π/ℓ < q < `--sq-band`·2π/ℓ (default 8), with all segments and without folded ones, compared with the global S(q) in the same band (`local_sq.npz`, `fig_local_sq.png`). When the local tilt is large, the exponent depends on the method; compare them by range of ℓ. It also computes the **height distribution** P(δu/σ) of the deviations of every point of the wall from its mean position, on the base plane of each frame (rotation by θ0 plus the per-frame line), normalised by the width of each frame, with skewness and excess kurtosis (block-bootstrap errors); variants: mean only (no per-frame line), the three heights, and local distributions in windows of 5 and 20 µm with a line removed in each. Outputs `height_distribution.npz` and `fig_height_distribution.png`.

### `step8_crossover.py`: crossover length l₀

`python step8_crossover.py DIR [DIR2 ...] [--m 2] [--n-boot 200]`, `--selftest` for synthetic walls with a known l₀. Reads `local_width.npz` and `roughness.npz` from step7 (with their per-time-block arrays) and estimates the crossover length between the small-scale (qEW-like, ζ₁ ≈ 1.25) and the large-scale (ζ₂ ≈ 0.5) regimes in three ways: (a) a smooth broken power law w² = A ℓ^{2ζ₁}[1 + (ℓ/l₀)^m]^{(2ζ₂−2ζ₁)/m} fitted to the local width with local rotation; (b) the same crossover in q fitted to the global S(q); (c) the scale where the local slope of w² crosses the midpoint between its small- and large-scale values. Variants (fixed or free exponents, m = 1, 2, 4, lower cutoff 2/3/5 µm, contour / folds excluded / u_area) give the systematic spread; a block bootstrap over time gives the statistical error; when l₀ is not bracketed a lower bound is reported. Also the disorder proxies l_fold and l_tilt (where the folded fraction and the sd of the local tilt fall to half their small-scale values). Outputs `crossover.npz`, `fig_crossover.png`, and with several studies `fig_crossover_compare.png`. The helpers shared by steps 7–9 are in `wall_tools.py`.

### `step9_wavelet.py`: roughness exponent from wavelets in local frames

`python step9_wavelet.py DIR [--min-kept 0.98] [--curv-ratio 0.6] [--n-boot 300]`, `--selftest` for synthetic walls (also with an optical blur). Uses the estimator of `wall_exponents.py` (`local_frame_windows`, `dog_kernel`, as in `wavelet_rms`). The curves are the subpixel contours saved by step7 in `contours.npz` (TIFF studies; contour of the domain smoothed by 1 px, without the pixel staircase; not tracked by git, regenerate with step7) or the curves (column, h) from `h_xt_sub.npy`. Each curve is resampled at unit arc length. For each scale a, windows of ~26a points are put in their own PCA frame; windows where t is not strictly monotonic are discarded; W = L1-normalised DOG-3 / DOG-4 coefficient at the window centre; F(a) = √⟨W²⟩ ∼ a^ζ (ℓ ≈ 4a). A scale is used only if the kept fraction is ≥ `--min-kept`, the window length 26a ≤ `--curv-ratio` × R_c (10th percentile of the radius of curvature of the curve smoothed over 13a points) and 4a ≥ rmin of step7. The fit range is the longest run of usable scales whose local slopes agree within ±0.15. Errors come from a block bootstrap over 6 contiguous segments × 12 time blocks. Comparison with the local width and the global S(q) over the same ℓ range. Outputs `wavelet.npz`, `fig_wavelet.png`. The optical blur biases ζ upwards at a ≲ 2× the PSF width (see the self-test). For very tortuous walls (Study 2) no scale passes the criteria and no ζ is reported.

A self-contained description of all the roughness methods (local ζ, crossover l₀, geometric descriptors), with their validation and pitfalls, is in [`ROUGHNESS_METHODS.md`](ROUGHNESS_METHODS.md).

### `step10_universality.py`: qEW or qKPZ below l₀

`python step10_universality.py DIR [DIR2 ...] [--lmin 4] [--d-def 8] [--n-boot 300]`. Tests on the measured walls, at the scales between the resolution and the crossover length l₀ of step8, where the creep avalanches are depinning-like. **A. Roughness:** free power law of the local width with local rotation (step7; contour, folds excluded, u_area) over [`--lmin` (3/4/5 µm, above the PSF bias found in step9), l₀], and the broken power law of step8 (ζ₂ = 0.5, m = 2, l₀ free) with ζ₁ fixed to 1.25 (qEW), 0.63 (qKPZ) or 2/3 (creep equilibrium), or free; χ² per model and block-bootstrap CIs. **B. KPZ term:** the normal velocity of contour segments of arc length L (contours of step7, switched pixels of the arrival map assigned to the nearest contour point, segments near defects or the border dropped; also without the defect mask, and from h(x,t) for all studies) vs the local orientation, fitted as v_n = a[1 + b δφ + (κ/2) δφ²]. Λ = λ/v = 1 + κ is the coefficient of v(m) = v₀[1 + (Λ/2)m²] for the velocity projected on the mean direction of motion: 0 for an elastic line driven along h (qEW), 1 for an isotropic medium (geometric term, λ = v → 0), ≫ 1 for qKPZ (λ finite when v → 0). R(L) = |Λ|⟨m_L²⟩/2 is the relative velocity change that the tilts of scale L produce. Three lags dt (mean advance 2, 5 and 15 px). **C.** Mean advance S/ℓ_y vs ℓ_y of the step4 avalanches (resolution-limited, not used for the decision). **D.** Skewness of the local heights from step7 (+ = bumps ahead). Outputs `universality.npz`, `fig_universality.png`, and with several DIRs `fig_universality_compare.png` in the first one. Summary of all studies: [`UNIVERSALITY.md`](UNIVERSALITY.md).

### `step11_overhangs.py`: onset of overhangs, l₀^geom

`python step11_overhangs.py DIR [DIR2 ...] [--tol 3] [--n-boot 300]`; `--selftest` for synthetic walls. Needs `contours.npz` of step7 and, for the comparison, `crossover.npz` of step8. ζ ≈ 1/2 is expected both above l₀ (effective, the wall is multivalued) and above the creep avalanche scale ξ_av (genuine, the wall is single-valued), so the break of step8 is min(ξ_av, l₀) and the exponent alone cannot tell which. This step measures directly whether the wall is single-valued at each scale. The ordered contour of each frame is cut into windows of arc length λ (sliding by λ/2), each rotated to its own principal axis. A window has an overhang if its projection on that axis goes back by more than `--tol` px (3 px by default, above the pixel staircase and the 1 px smoothing; also 2 and 5 px). **P_ov(E)** = fraction of windows with an overhang vs their extent E along their own axis (the ℓ of step7/8). **l₀^geom** = E where P_ov = 1/2, with a 16–84 % block bootstrap over time. Control: back steps on the base plane per unit length (mean spacing between overhangs). **P_ov at the step8 crossovers:** ≈ 1/2 means the break is the overhang crossover; ≪ 1/2 means the wall is still single-valued there (candidate ξ_av). Windows longer than about half the wall are few and overlap, so a single large fold can dominate P_ov there. Self-test: no false positives on a single-valued tilted wall (ζ = 1); with overhangs every D px, l₀^geom ≈ 0.4 D for D = 100–400 px. Outputs `overhangs.npz`, `fig_overhangs.png`; summary table on stdout.

### `bubble_event_map.py`: expanding bubble from an avalanche event map

`python bubble_event_map.py EVENTS.csv --outdir OUT [--snapshots 300] [--rmin 8] [--rmin-sq 20] [--sector 1] [--smax-tau 30] [--tol 3]`. For data given not as a movie but as an event map: a CSV image where each pixel holds the time-ordered index of the avalanche in which it switched, or −1. There are no frame times, so the event index is the clock and lengths are in px. The −1 pixels are split into the initial bubble (the component at the centroid of the switched area), the outside (touching the border) and holes, which get the index of the nearest switched pixel. This gives an arrival map, from which `--snapshots` domains equally spaced in switched area are taken. For each one: the polar profile r(θ) of the subpixel contour (1 px Gaussian smoothing) on a periodic grid of `--sector` px of arc, and its S(q), fitted between 2π/(L_min/4) and 2π/`--rmin-sq`. Also the local width on arc windows rotated to their own axis (`wall_tools.overhang_segments`, as step11), with the variance of a perfect arc subtracted, for windows up to the smallest mean radius; and P_ov with the same windows. Avalanches: P(S) with a maximum-likelihood τ truncated to [5, `--smax-tau`] px (the untruncated fit is dominated by the cutoff), and extents ℓ_∥, ℓ_⊥ along and across the wall tangent at each avalanche's barycentre. Outputs `bubble.npz`, `fig_bubble_map.png`, `fig_bubble_roughness.png`, `fig_bubble_avalanches.png`.

### `tiff_arrival_map.py`: arrival map from TIFF frames, for faint slow walls

`python tiff_arrival_map.py TIFF_DIR --outdir DIR [--drift DIR/drift_raw.npy] [--offset 33000] [--block 10] [--x0 --x1 --y0 --y1] [--sign 1] [--sigma 1.5]`. The step-function arrival map of `arrival_map.py`, but read from a TIFF folder, registered with the drift of an earlier step1 run, normalised by the frame median and averaged in blocks of frames to raise the signal. Saves `ta` (Otsu-thresholded), `ta_all` (every pixel) and `amp`. For walls whose domain contrast is too low for step1 (Study 15); there the Otsu threshold is set by defects and the swept band has to be selected by hand (ridge tracking, see the Study 15 report).

### `arrival_map.py`: arrival map of a movie (bubbles, several walls)

`python arrival_map.py MOVIE --outdir DIR [--sign 1] [--sigma 1] [--detrend] [--frames F0 F1] [--qc K1,K2,...]`. For geometries that step1 does not handle (bubbles, several walls, walls crossing a circular aperture). For every pixel, the frame of the largest step of its grey level (step-function cross-correlation, as in Durin et al., PRB 110, L020405 (2024)): maximise the signed difference of the means after and before the step times √(n₁n₂/N). Frames smoothed by `--sigma` px; pixels with a step amplitude below the Otsu threshold never switched (−1). `--detrend` removes a global brightness change (camera gain, lamp): each frame is shifted by the cumulative median of the frame-to-frame change over the aperture. Outputs `arrival.npz` (ta, amplitude, mean frame), `meta.json` and `fig_qc_overlay.png` (switched region {0 ≤ ta ≤ k} over frame k), to be approved before the analysis.

### `contour_analysis.py`: roughness, overhangs and KPZ term from an arrival map

`python contour_analysis.py DIR/arrival.npz DIR --frames K0 K1 [--um-per-px s] [--centre y x] [--exclude-deg a b] [--roi-erode px] [--ddef 8] [--every n] [--skip-kpz]`. Per frame, the subpixel contour of {0 ≤ ta ≤ k} (1 px smoothing); points outside the ROI (bright aperture of the mean frame eroded by `--roi-erode`), within `--ddef` px of a never-switched island (defect) or, for a bubble (`--centre`), inside the angular sector `--exclude-deg` are dropped, and the contiguous pieces ≥ 60 px kept. Then:
- local width with local rotation (`wall_tools.overhang_segments`, as step7/step11), minus the variance of a perfect arc for a bubble; ζ over E = 8–80 px and 20–80 px; P_ov(E) at tol 2, 3, 5 px and l₀^geom;
- KPZ term as step10 B: segments of arc length L = 8, 16, 32 px, the pixels switching in (k, k+dt] assigned to the nearest contour point, v_n = count/(L dt), tilt relative to the wall at scale 8L; v/⟨v⟩ = a[1 + b Δφ + (κ/2)Δφ²], Λ = 1 + κ (0 elastic qEW, 1 isotropic, ≫ 1 qKPZ), 12-block bootstrap;
- bubble only: macroscopic normal velocity vs orientation from the radial advance between frames (Fourier fit, modes n ≤ 2 and ≤ 4, on the allowed angles; noise from random rotations of each interval), Krug–Spohn Λ(φ) = (v + v″)/v (a uniform drift, n = 1, cancels); skewness of r(θ) minus modes n ≤ 4.
Outputs `contour_analysis.npz`, `fig_contour_analysis.png`.

### `synthetic_bubble.py`: synthetic bubbles of known ζ (calibration)

`python synthetic_bubble.py ZETA SEED --out FILE.npz|.csv [--w50 1] [--r0 106] [--r1 400]`. A bubble whose wall has S(q) ∼ q^−(1+2ζ) (or a broken spectrum `Z1:Z2:LC`), fixed in θ and rescaled during the growth (nested domains, amplitude stationary in arc length), with local width `--w50` px over 50 px of arc; written as an arrival map (`contour_analysis.py`) or an event map (`bubble_event_map.py`). Used to calibrate the local-width exponent of walls only ~1 px rough (Studies 12–13).

### `phi4_contour_roughness.py`: roughness of φ⁴ critical walls, bubble or straight

`python phi4_contour_roughness.py DIR [DIR2 ...] --L L [--Ly LY] [--bubble R0] [--block B] [--outdir OUT]`. Reads the critical fields written by [vmc-phi4-depinning](https://github.com/droyktton/vmc-phi4-depinning) with `--dump-field` (`criticafield_*.bin`, raw float32) and takes the φ = 0 contour: the spanning wall (periodic in x) or, with `--bubble R0`, the closed wall of a bubble (points within 3 sites of the initial circle dropped). `--block B` first binarises B × B blocks, as a camera pixel. Local width of rotated arc windows and P_ov(E) as `contour_analysis.py`, minus the variance of an arc of the mean radius for a bubble; ζ over several ranges with a bootstrap over samples. Outputs `phi4_roughness.npz`, `fig_phi4_roughness.png`. Used for the bubble-vs-straight test in `UNIVERSALITY.md` (figure in `phi4_bubble/`).

## Figures

| File | What it shows |
|---|---|
| `fig_qc_overlay.png` | Segmentation check: subpixel h(x,t) in red over 6 frames spread across the frames used. Use it to confirm that the curve follows the real wall contour. |
| `fig_h_mean.png` | Left: mean wall position H − ⟨h⟩_x against frame, to check that the advance is monotonic with no jumps. Right: mean velocity −Δ⟨h⟩ in px/frame. Dropped frames, if any, are shaded in grey. |
| `fig_Cn_tau_<tag>eps<ε>.png` | Left: C(n,τ) against distance n for the first τ values with a valid fit. Dots are data and lines are the fit A·e^(−n/ξ)+B. Right: (C−B)/A on a semilog scale, where an exponential decay appears as a straight line of slope −1/ξ. |
| `fig_xi_chi4_<tag>eps<ε>.png` | Top row linear, bottom row log-log. Columns: (1) ξ(τ) with fit error bars, in px on the left axis and µm on the right; (2) χ4(τ) computed directly, from the C(n) sum and from the model sum, with τ\* marked in red; (3) Π(τ); (4) normalized χ4/[Π(1−Π)] (dots), its local part (χ4 − B·L)/[Π(1−Π)] (squares) and 2ξ (dashed) for comparison. |
| `fig_eps_sweep_sub.png` | (1) ξ(τ) for each ε on a log τ axis, with a µm axis on the right; (2) χ4(τ) for each ε; stars mark τ\* in both; (3) ξ(τ\*), ξ_max and ξ(τ=1) against ε, with χ4(τ\*) in red on the right axis; (4) τ\*, the lag of ξ_max and the lag at which Π = ½ against ε. The grey band extends to ε = 3σ_Δh. |
