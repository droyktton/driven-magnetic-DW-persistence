# Driven magnetic domain wall: persistence and dynamic correlation length

Analysis pipeline for movies of a driven magnetic domain wall, for example from polar MOKE microscopy. The pipeline:

1. Extracts the wall position h(x,t).
2. Computes the persistence field between times t and t+τ.
3. Measures the dynamic correlation length ξ(τ) and the four-point susceptibility χ4(τ).

The included example is `magnetic_fliped/`: 78 frames of 768×370 px at 25 fps (Δt = 40 ms). The scale is δ ≈ 0.17 µm/px, so the field of view is ≈ 131 × 63 µm. The bright domain is on top, the dark domain at the bottom, and the wall moves upward.

Requires numpy, scipy, matplotlib, scikit-image and imageio-ffmpeg. imageio-ffmpeg ships its own ffmpeg binary.

## Usage

Each movie gets its own output folder. By default the folder is named after the video file without its extension.

```bash
python step1_extract_wall.py MOVIE.mp4 --um-per-px 0.17        # -> MOVIE/ with h(x,t) and meta.json
python step2_persistence.py MOVIE --eps 0.25 0.5 0.75 1 1.5 2 3  # subpixel h
python step2_persistence.py MOVIE --int                          # integer h, ε = 1, 2
python step3_eps_sweep.py MOVIE                                  # sensitivity to ε
```

Step 1 writes the frame rate and the spatial scale to `MOVIE/meta.json`, and steps 2 and 3 read them from there. By default, the frame rate is the one stored in the video and the scale is 0.17 µm/px. You can override either with `--fps` or `--um-per-px` in any step. `--um-per-px 0` reports everything in pixels only.

The original video of the example is not in the repository. Without it you can start from step 2, because `magnetic_fliped/h_xt*.npy` and `meta.json` are included. To rerun step 1, copy `magnetic_fliped.mp4` into the repository root.

## Input requirements

The method assumes the following about the movie. Step 1 checks the ones marked ✔ and prints warnings, which are also stored in `meta.json` and repeated by step 2.

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

## Figures

| File | What it shows |
|---|---|
| `fig_qc_overlay.png` | Segmentation check: subpixel h(x,t) in red over 6 frames spread across the frames used. Use it to confirm that the curve follows the real wall contour. |
| `fig_h_mean.png` | Left: mean wall position H − ⟨h⟩_x against frame, to check that the advance is monotonic with no jumps. Right: mean velocity −Δ⟨h⟩ in px/frame. Dropped frames, if any, are shaded in grey. |
| `fig_Cn_tau_<tag>eps<ε>.png` | Left: C(n,τ) against distance n for the first τ values with a valid fit. Dots are data and lines are the fit A·e^(−n/ξ)+B. Right: (C−B)/A on a semilog scale, where an exponential decay appears as a straight line of slope −1/ξ. |
| `fig_xi_chi4_<tag>eps<ε>.png` | Top row linear, bottom row log-log. Columns: (1) ξ(τ) with fit error bars, in px on the left axis and µm on the right; (2) χ4(τ) computed directly, from the C(n) sum and from the model sum, with τ\* marked in red; (3) Π(τ); (4) normalized χ4/[Π(1−Π)] (dots), its local part (χ4 − B·L)/[Π(1−Π)] (squares) and 2ξ (dashed) for comparison. |
| `fig_eps_sweep_sub.png` | (1) ξ(τ) for each ε, with a µm axis on the right; (2) χ4(τ) for each ε on a semilog scale; (3) ξ(τ=1) in µm (blue) and χ4(τ=1) (red) against ε. The grey band extends to ε = 3σ_Δh. |

## Results for `magnetic_fliped`

- **Noise:** σ_Δh = 0.25 px, so ε = 3σ_Δh = 0.75 px was chosen (subpixel h). Only 0.02 % of the one-frame displacements are retreats larger than 1 px.
- **τ\*:** χ4(τ) decreases monotonically from τ = 1 for every ε, so **τ\* ≤ 1 frame = 40 ms** and is not resolved at this frame rate. The wall advances ~4 px/frame, so Π(τ) is almost zero by τ ≈ 10.
- **ξ(τ\*)** at ε = 0.75 px: **17.7 ± 0.2 px = 3.01 ± 0.04 µm** (statistical error). The systematic error is ~±3 px (~±0.5 µm), because ξ(1) grows from 14 to 22 px (2.4 to 3.8 µm) as ε goes from 0.25 to 3 px.
- **Mean wall velocity:** 3.95 px/frame = 16.8 µm/s.
- **Shape of ξ(τ):** a plateau from τ = 1 to 2–3 frames, followed by a decay to ~5 px at τ ≈ 7–8.
- **Normalized χ4** (ε = 0.75 px): χ4/[Π(1−Π)] ≈ 45–48 px for τ = 1–3, then decays. It has the same shape as ξ(τ). So the monotonic decay of the raw χ4 comes mostly from Π → 0, not from a loss of cooperativity at short τ. The plateau term B·L accounts for about a third of χ4 at τ = 1.
