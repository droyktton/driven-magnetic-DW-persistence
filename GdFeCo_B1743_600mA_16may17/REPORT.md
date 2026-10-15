# Study 10: `GdFeCo_B1743_600mA_16may17` (GdFeCo B1743, 600 mA, Δt = 2 s, 6 h)

First report, 2026-10-15. Pipeline and definitions: see the [README](../README.md). First study of a ferrimagnetic sample; independent of Studies 1–9.

## Data

- **Source:** `/media/…/MOBILEDATAMOKE/MOKE 2017/Data2017/2017_05_16/B1743_Avalanchas/B1743_20x_IzBG1_600mA_ImagenesCada2s_1/Pos0/`, copied to `~/MOKEdata_local/GdFeCo_B1743_2017_05_16/Pos0` (not in the repository).
  - 10800 Micro-Manager TIFF frames, 1392×1040 px, Retiga Blue camera, one wall.
  - Recorded 2017-05-16.
- **Sample B1743:** Si/SiO₂(100) / Ta 5 / Gd₃₃Fe₆₀.₃Co₆.₇ 10 / Pt 5 (nm), perpendicular anisotropy, coercive field 74 Oe.
  - Grown at Berkeley (C.-H. Lambert, J. Gorchon), 2016.
  - Composition from `GdFeCo/20161102_Samples_Javier_001.docx` on the same drive.
- **Drive:** 600 mA in coil "BG1" (`IzBG1` in the folder name). The field per ampere of this coil is not known, so the drive is given in mA only.
- **Δt = 2.000 s:** from the `metadata.txt` timestamps (`Interval_ms` = 2000; 21 599 s between the first and the last frame).
- **Scale 0.1175 µm/px: assumed**, the 20x value of Studies 2–7 (same setup and camera). All lengths in µm scale with it; exponents and Λ do not.
- **Geometry:** the bright domain grows to the right from the left edge of the aperture. Used `--bright-domain`, no rotation.

```bash
python step1_extract_wall.py ~/MOKEdata_local/GdFeCo_B1743_2017_05_16/Pos0 --outdir GdFeCo_B1743_600mA_16may17 \
       --um-per-px 0.1175 --bright-domain --offset auto --flatfield 2 --init-wall --crop-right 1200 --qc-every 100
python step2_persistence.py GdFeCo_B1743_600mA_16may17 --row-mean --eps 0.5 1 1.5 2 3 5 --tau-max 2000
python step7_roughness.py GdFeCo_B1743_600mA_16may17 --every 36 --rmin 25
python step8_crossover.py GdFeCo_B1743_600mA_16may17
python step10_universality.py GdFeCo_B1743_600mA_16may17
```

## Wall extraction and checks

**Segmentation check.** `fig_qc_overlay.png`, `qc/overlay_*.png` (1 frame in 100) and `qc/overlay_small.mp4` were approved by the user before the analysis.

- **`--crop-right 1200` is needed.** The automatic crop (x = 904, where the illumination falls below 90 %) left only 26 % of the field as never-switching. The flat-field fit and the Otsu threshold then failed: negative threshold, crash in the QC figure. With 1200 the contrast map is clean (threshold 0.082, 97 defect px).
- **Rows:** 444–983 (540 rows = 63 µm). Above row 444 the initial wall runs along the aperture edge; the domain only reaches those rows later, by lateral growth.
- **Drift:** dy = −0.4…+0.1 px, dx = −0.7…0 px.
- **Motion:** all 10800 frames used. Total advance 433 px = 51 µm, **mean velocity 0.04 px/frame ≈ 2.4 nm/s**, the fastest of all studies; no retreats > 1 px.
- **Noise:** σ_Δh = 0.38 px (step2), so 3σ_Δh ≈ 1.1 px.
- **Tilt:** θ = 27° at the start, 7° at the end; θ0 = 21.9° in step 7.
- **Overhangs:** the front differs from the area estimate in 11 % of the (row, frame) pairs (P7925: 0.4 %); on the base plane of step 7, 0.7 %. The wall forms bays around defects (e.g. y ≈ 650–800 at 4.8 h).

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted, τ ≤ 2000 frames (67 min). τ in frames of 2 s.

| ε [px] | τ\* (max χ4) | χ4(τ\*) | ξ(τ\*) [µm] | ξ_max [µm] at τ |
|---|---|---|---|---|
| 0.5 | 163 (5.4 min) | 10.6 | 4.01 ± 0.07 | 4.09 at 7.8 min |
| 1 | 207 (6.9 min) | 16.5 | 4.40 ± 0.06 | 5.06 at 29 min |
| **1.5** | **218 (7.3 min)** | 17.1 | **4.46 ± 0.06** | 5.56 at 44 min |
| 2 | 250 (8.3 min) | 17.3 | 4.44 ± 0.07 | 5.84 at 51 min |
| 3 | 324 (11 min) | 17.6 | 4.36 ± 0.08 | 5.71 at 52 min |
| 5 | 335 (11 min) | 18.2 | 4.41 ± 0.08 | 5.48 at 64 min |

- **ξ(τ\*) ≈ 4.4 µm** for ε ≥ 1 px, almost independent of ε. It is about twice the value of Pt/Co/Pt and P7925 (1.6–2.6 µm), and ξ_max ≈ 5.5 µm (others 2.4–3.3 µm).
- **χ4(τ\*) ≈ 17**, twice that of P7925 (≈ 8): larger correlated regions of persistent sites.
- **τ\* ≈ 7–11 min**, against 35–50 min in P7925; the wall is ~8× faster. In units of advance, ⟨v⟩τ\* ≈ 9–13 px, similar to P7925 (≈ 5–8 px).

## Roughness and l₀

`step7` (`fig_roughness.png`, `fig_local_width.png`, `fig_local_sq.png`): 402 rows on the base plane (47 µm), 300 frames, fit range 2.9–12 µm.

| estimator | ζ |
|---|---|
| global S(q), u_area | 0.97 ± 0.12 (thirds: 1.41, 1.08, 0.74) |
| local width, local rotation (contour) | 1.17 ± 0.03 (folds excluded 1.29, tilts > 45° excluded 1.32) |
| local S(q, l), l = 12–17 µm | 0.9–1.2 (unfolded 1.0–1.3) |
| B(r) | 0.64 (biased low) |

- The global ζ drifts from 1.4 to 0.7 over the run, as the wall straightens (tilt 27° → 7°) and the field it crosses changes.
- **Crossover l₀** (`step8`, `fig_crossover.png`): 35.5 µm from the local width (bootstrap 31–41), 17.7 µm from the global S(q), 22.3 µm from the local slope. l_fold = 17.8 µm, l_tilt = 29.6 µm. **l₀ ≈ 18–35 µm.**
- **Height distribution:** skewness −0.40 ± 0.13 on the base plane and −0.6 ± 0.1 in 5 µm windows (lagging valleys); excess kurtosis +5 in 5 µm windows.

## qEW or qKPZ below l₀

`step10_universality.py` (`fig_universality.png`); comparison of all studies in [UNIVERSALITY.md](../UNIVERSALITY.md).

**A. Roughness below l₀** (ζ over [4 µm, l₀ = 35.5 µm]):

| estimator | ζ |
|---|---|
| contour | 1.04 (1.02–1.07) |
| folds excluded | 1.17 (1.14–1.21) |
| u_area | 1.02 (1.00–1.05) |

- Broken power law, folds excluded: ζ₁ = 1.31 (1.26–1.37), l₀ = 27 µm; with the PCA contour ζ₁ = 1.22 (1.17–1.27).
- **χ² qEW / qKPZ / equilibrium = 13 / 1266 / 1097** (folds excluded); 5 / 418 / 343 on the contour.
- ζ moves by ≤ 0.06 for ℓ_min = 3–5 µm.

**B. KPZ term:**
- Λ from the contour: 0.1–0.75 for L = 2–8 µm; 1.3–2.6 at L = 16 µm.
- Λ from h(x): 0.2–0.4 for L = 2–4 µm, −0.5…−0.2 at 8 µm, −1.5…−1.2 at 16 µm.
- At L = 16 µm the two estimates have opposite signs, and only ~700 segments: not a robust KPZ signal.
- The defect mask changes Λ by ≤ 0.6. **R ≤ 0.24.**

**Conclusion: qEW**, as in Pt/Co/Pt and P7925, with ζ ≈ 1.0–1.3 below l₀ and Λ = O(1).

## Caveats

- **Scale assumed** (0.1175 µm/px). If it is wrong, all µm values, v and l₀ scale by the same factor; ζ, Λ and R do not change.
- **Field unknown:** the calibration of coil BG1 is not in the data folder.
- **Short wall:** 63 µm of wall (47 µm on the base plane), so l₀ ≈ 20–35 µm leaves less than a decade above 4 µm.
- **Not run yet:** steps 3–6 and 9.
