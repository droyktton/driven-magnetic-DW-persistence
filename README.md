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

Requires numpy, scipy, matplotlib and scikit-image; for videos also imageio-ffmpeg (or an `ffmpeg` binary on the PATH); for TIFF folders also tifffile and Pillow, and ffmpeg for the QC movie.

## Usage

Each movie gets its own output folder. By default the folder is named after the video file without its extension. For a TIFF folder it is named after the folder, or after its parent when the folder is a Micro-Manager `PosN`.

```bash
python step1_extract_wall.py MOVIE.mp4 --um-per-px 0.17        # -> MOVIE/ with h(x,t) and meta.json
python step1_extract_wall.py /path/RUN/Pos0 --um-per-px 0      # TIFF frames -> RUN/
python step2_persistence.py MOVIE --eps 0.25 0.5 0.75 1 1.5 2 3  # subpixel h
python step2_persistence.py MOVIE --int                          # integer h, ε = 1, 2
python step3_eps_sweep.py MOVIE                                  # sensitivity to ε
python step4_avalanches.py MOVIE --tau-m 1 3 10                  # avalanches (TIFF studies)
```

Step 1 writes the frame rate and the spatial scale to `MOVIE/meta.json`, and steps 2 and 3 read them from there. By default, the frame rate is the one stored in the video (for TIFF: from the timestamps in Micro-Manager's `metadata.txt`) and the scale is 0.17 µm/px. You can override either with `--fps` or `--um-per-px` in any step. `--um-per-px 0` reports everything in pixels only.

The original video of the example is not in the repository. Without it you can start from step 2, because `magnetic_fliped/h_xt*.npy` and `meta.json` are included. To rerun step 1, copy `magnetic_fliped.mp4` into the repository root.

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
4. **Arrival-time map.** For each pixel, the number of frames with s > 0.5 is the frame at which it switches (this assumes a monotonic advance). Defects inside the swept area show no contrast; they get the arrival time of the nearest valid pixel.
5. **Wall.** The domain at frame t is the set of pixels with arrival ≤ t connected to the left edge, with holes filled, so isolated spots ahead of the wall are ignored. The wall position x(y,t) is the right-most domain pixel of each row y. Where the wall folds back (overhangs), this is its front. The subpixel position comes from the 0.5 crossing of s. Then h = W − 1 − x, so that h decreases with t as in the video convention, and the columns of h are the image rows y.

h is measured along x in each row y. For a wall tilted by θ from vertical, a step Δx corresponds to a normal displacement Δx·cos θ, and a distance n along y to n/cos θ along the wall. The tilt θ(t) is fitted in every frame and reported.

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

Reads every `MOVIE/persistence_sub_eps*.npz` file, or the integer-h files with `--int`. It overlays ξ(τ) and χ4(τ) for each ε and plots ξ(τ=1) and χ4(τ=1) against ε. σ_Δh and µm/px are read from the `.npz` files. Output: `MOVIE/fig_eps_sweep_sub.png`, or `fig_eps_sweep_int.png` with `--int`.

### `step4_avalanches.py`: avalanche statistics (TIFF studies)

Reads the arrival-time map in `DIR/tiff_extra.npz`. An avalanche is a connected patch (8-neighbours) of pixels that switch within the same measurement window τ_m (`--tau-m`, in frames; several values can be given). Patches smaller than `--smin` px (default 20) are discarded. For each avalanche it stores the area S, the extents ℓ_y (along the wall) and ℓ_x, the centroid and the frame, and flags those touching the image edge or the never-reached region.

- P(S) is fitted by maximum likelihood with P(S) ∝ S^(−τ)·exp(−S/S_cut) above `--smin`, with errors from the inverse Hessian. A pure power law (no cutoff) is also printed, for reference only.
- S ∝ ℓ_y^(1+ζ) is fitted on the median S in log bins of ℓ_y.
- `--roi X0 X1 Y0 Y1` keeps avalanches whose centroid is in a region (camera px); `--compare FILE` overlays another list of areas in px (second column), e.g. the thesis list; `--tag` adds a suffix to the output files.

Outputs: `avalanches<tag>_tm<m>.npz`, `fig_avalanches<tag>.png` (P(S) with fits, P(ℓ_y), S vs ℓ_y) and `fig_avalanche_map<tag>.png` (avalanches of the first window coloured by frame; grey: area that switched in steps smaller than `--smin`).

## Figures

| File | What it shows |
|---|---|
| `fig_qc_overlay.png` | Segmentation check: subpixel h(x,t) in red over 6 frames spread across the frames used. Use it to confirm that the curve follows the real wall contour. |
| `fig_h_mean.png` | Left: mean wall position H − ⟨h⟩_x against frame, to check that the advance is monotonic with no jumps. Right: mean velocity −Δ⟨h⟩ in px/frame. Dropped frames, if any, are shaded in grey. |
| `fig_Cn_tau_<tag>eps<ε>.png` | Left: C(n,τ) against distance n for the first τ values with a valid fit. Dots are data and lines are the fit A·e^(−n/ξ)+B. Right: (C−B)/A on a semilog scale, where an exponential decay appears as a straight line of slope −1/ξ. |
| `fig_xi_chi4_<tag>eps<ε>.png` | Top row linear, bottom row log-log. Columns: (1) ξ(τ) with fit error bars, in px on the left axis and µm on the right; (2) χ4(τ) computed directly, from the C(n) sum and from the model sum, with τ\* marked in red; (3) Π(τ); (4) normalized χ4/[Π(1−Π)] (dots), its local part (χ4 − B·L)/[Π(1−Π)] (squares) and 2ξ (dashed) for comparison. |
| `fig_eps_sweep_sub.png` | (1) ξ(τ) for each ε, with a µm axis on the right; (2) χ4(τ) for each ε on a semilog scale; (3) ξ(τ=1) in µm (blue) and χ4(τ=1) (red) against ε. The grey band extends to ε = 3σ_Δh. |
