# TbCo [Co(1.7)/Tb(t)]×5, LPS Paris, Nov 2023: pulsed-field walls (in progress)

Data by Vincent Jeudy and Jean Rodriguez (LPS, Nov 2023). Wafer AV6517, 7 positions with Tb thickness 0.57–1.08 nm, room temperature. **Status:** data organised, 0.83 nm analysed on an earlier segmentation, segmentation of all 7 thicknesses not yet approved. The per-sample outputs on disk come from different code versions and are not committed (see the end).

## Data

- Raw data: `~/MOKEdata_local/TbCo_Paris/` (3.2 GB), a copy of `MOBILEDATAMOKE/MOKE-TbCo_Paris` (the external disk is unreliable).
- Hamamatsu, 16-bit, 2048×2048. **Raw images, not difference images**, one per field pulse. X20: 0.136 µm/px (scale bar 50 µm = 366 px). The dodecagon is the field diaphragm, opened in some 0.74 nm frames.
- Pulses: HP8114, 1–50 µs, 180–950 V, **µ0H = 0.205 mT/V × V**. The 0.175 mT/V column of the logbook is for the other generator (HV_1000), used only in a few early 0.66 nm files. The 1.08 nm sample also has a Kepco bias of 23.4 mT (2.81 V); Kepco pulses of 3–600 s.
- `logbook_transcription.md`: the handwritten logbook (`Bitacora_Paris_04-12-2023 16.37.pdf`) transcribed: file, voltage, pulse length and nucleation pulses per sample.
- `origin_tables/*.csv`: the worksheets of Jean's Origin project (`PtCoTb_DW_Jean_2023_11_VJ.opju`, not readable on Linux), exported to PDF and parsed. One row per image, in logbook order. Columns include the hand-measured displacement `disp_px` with its error, v = 0.136·disp/dur, h = H/H_d and the reduced velocity vred = v/v(H_d)·(T/T_d)^0.15.

H_d and v(H_d) from their creep/depinning fits (read back from the h = H/H_d column):

| Tb (nm) | 0.57 | 0.66 | 0.74 | 0.83 | 0.91 | 1.00 | 1.08 |
|---|---|---|---|---|---|---|---|
| H_d (mT) | 48 | 51 | 69 | 83 | 114 | 155 | 217* |
| v(H_d) (m/s) | 29 | 25 | 22 | 18 | 15 | 12 | 8.4 |

\* including the 23.4 mT Kepco bias. All seven samples collapse on one reduced curve in Origin.

**Origin displacement ("Criterio 3"):** a local, near-maximum forward advance measured by hand, not a wall average (confirmed with Jean's description).

## Method

`tbco_pulses.py` (repo root) finds the walls before and after each pulse:
- Groups frames by diaphragm setting and registers them on sample defects.
- Normalises the illumination, chained between consecutive images.
- Builds per-pixel bright and dark levels, then takes the subpixel wall at half contrast.
- Masks static defects and low-contrast areas.

For each pulse (image k−1 → k) it computes:
- the swept area, both from the masks and from the consecutive difference;
- the wall-mean displacement (area / wall length);
- the local displacement of each wall point (median, 90th percentile, max);
- the rotated local width w²(E);
- the KPZ term Λ, as in step10 B but on the contour.

`tbco_gif.py` makes the overlay GIFs (red: wall after the pulse, cyan: before).

## Results, 0.83 nm (H_d = 83 mT, 42 pulses, 68–193 mT)

**Velocity**
- Two independent segmentations (masks and consecutive differences) agree to about 1 % in swept area (overlap 0.93–1.00) wherever the mean advance is above about 3 px (0.4 µm).
- Below that (files 42–46 at 68–74 mT, and the 600 s Kepco pulse, file 55) we are at the noise level.
- Our wall-mean displacement is **0.93× Origin above H_d** but **0.32× (0.18–0.65) in creep**.
- The 90th percentile of the local advance reproduces Criterio 3 in creep: file 49 gives 148 px against 150; files 50–52 give 103/159/56 against 120/120/40.
- So their creep v(H) is the velocity of the forward-most parts of the wall. It is steeper and larger than the wall-averaged velocity, which matters for their creep fits, for H_d and for the collapse.

**Roughness** (rotated local width of the wall before each pulse)
- Creep walls (68–81 mT, 12 pulses) are about 3× rougher than walls above H_d.
- The apparent ζ drifts with the window:
  - about 1.6 below 2 µm (optical resolution about 1 µm, plus the σ = 2 px smoothing);
  - 1.3 at 3–4 µm;
  - about 1.05 at 5–9 µm;
  - about 0.9 at 11–18 µm.
- **Forward model** (`forward_model/`): synthetic walls passed through the measured edge blur (σ_psf ≈ 3–3.6 px) and noise (0.18 contrast units) of these images:
  - an apparent ζ of about 1.0–1.1 over 2–20 µm corresponds to a true ζ of about **0.7–0.8 between 5 and 20 µm**;
  - a qEW regime (ζ = 1.25), if present, can only lie below about 6 µm.
- The S(q) estimator on arc windows failed validation on synthetic walls and is not used.

**KPZ term Λ** (creep pulses only; above H_d the wall moves 3–13 window lengths per pulse, so Λ there is a geometric artefact)
- Λ = 0.05 (−0.14 to 0.25) at L = 2.2 µm.
- Λ = 0.44 (0.26 to 0.71) at 4.4 µm.
- Λ = −0.16 (−0.48 to 0.16) at 8.7 µm.

This is between elastic (0) and isotropic (1), with no qKPZ signal, as in CoFeB.

**Caveats**
- Only 12 creep pulses, from two series that enter from the left and lower-right edges.
- The field of view is about 150 µm.
- The creep walls branch, and overhangs have not been measured.

Summary figure: `AV6517_1_0_#2_0.83nm/fig_pulses.png`. It shows w²(E) for the three regimes, v_n against Δφ (Λ), and v(H) for Origin, the wall mean and the local 90th percentile. It was made from the earlier segmentation, whose `pulses.json` was later overwritten by the rewrite run.

## Segmentation of all 7 thicknesses (in progress)

The first pass and the rewrite of `tbco_pulses.py` are described above. The quality varies with the Kerr contrast, which is weak near compensation:

| Tb | pulses | state |
|---|---|---|
| 0.57 nm | 8 | contrast 0.09; false speckle in files 12–17; hardly usable |
| 0.66 nm | 46 (rewrite) | grey domains in files 16–29 traced; a static false contour on the right remains; reported shifts up to 37 px are suspect |
| 0.74 nm | 51 | five diaphragm states; state 3 (files 58–67) has contrast 0.007 and is unusable; opened-diaphragm frames 68–88 are messy |
| 0.83 nm | 42 (rewrite) | clean, but run before the registration and illumination fixes |
| 0.91 nm | 37 | good walls; drift 15–25 px; fixed speckle patch top right |
| 1.00 nm | 31 | good walls and bubbles |
| 1.08 nm | 27 (rewrite) | contrast 0.04–0.05; walls in the main disc fine; noisy contours along the dim edge of the open field |

Open problems:
- The code still marks 100 % of pixels as "switched".
- The 0.66 nm test needed about 5 GB of RAM: it needs a half-precision stack, restricted to the field region.

Next steps:
1. Reduce the memory use.
2. Check the 0.66 nm shifts and the switched-pixel criterion.
3. Rerun the 7 thicknesses one at a time.
4. Approve the overlays and GIFs.
5. Then compare velocity (wall mean vs local), roughness and Λ across thicknesses (H_d 48–217 mT).

## Files

**Committed**
- `tbco_pulses.py`, `tbco_gif.py` (repo root).
- `REPORT.md`, `logbook_transcription.md`, `origin_tables/*.csv`.
- `AV6517_1_0_#2_0.83nm/fig_pulses.png`.
- `forward_model/`:
  - `forward.py` and `measure.py`;
  - `measured.json`: edge width, noise and step per file;
  - `fwd_results.json`, `forward.log`, `fig_fwd_tbco.png`.

  The scripts read `pulses.json`/`walls.npz` of the 0.83 nm run and the B/D maps written by the earlier segmentation, which are not kept.

- **Segmentations, one folder per thickness `AV6517_<pos>_<Tb>/`:**
  - `fig_overlay_*.png`: every file with the wall before (cyan) and after (red) the pulse, and the analysed region (yellow);
  - `pulses.json`: per-pulse numbers.

  These come from two code versions and are not final:

  | Tb | run | version | pulses in `pulses.json` |
  |---|---|---|---|
  | 0.57 nm | 12 Nov | first pass | 8 |
  | 0.66 nm | 13 Nov | rewrite (diaphragm edge, per-pixel levels) | 50 |
  | 0.74 nm | 12 Nov | first pass | 53 |
  | 0.83 nm | 13 Nov | rewrite, before the registration and illumination fixes | 42 |
  | 0.91 nm | 12 Nov | first pass | 37 |
  | 1.00 nm | 12 Nov | first pass | 34 |
  | 1.08 nm | 13 Nov | rewrite | 27 |

  The quality of each is assessed in the table above. All seven will be rerun with the final code.

**Not committed**
- `walls.npz` (2–18 MB per sample) and `segmentation.gif` (4–13 MB): large, regenerated by `tbco_pulses.py` and `tbco_gif.py`.
- The raw images.
