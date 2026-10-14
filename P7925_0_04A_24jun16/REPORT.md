# Study 8: `P7925_0_04A_24jun16` (sample P7925, 0.04 A, Δt ≈ 9.5 s, 15.9 h)

First report, 2026-10-14. Pipeline and definitions: see the [README](../README.md). [Study 9](../P7925_0_04A_26jun16/REPORT.md) is the same sample and drive two days later. Independent of Studies 1–7.

## Data

- **Source:** folder `/media/…/MOBILEDATAMOKE/DataDWMoviesLargas/24jun16/`, copied to `~/MOKEdata_local/24jun16` (not in the repository).
  - 6000 TIFF frames, 1392×1040 px, one wall, sample P7925.
  - `info.txt`: I = 0.04 A, exposure 2.6 s, plus a focus-voltage table. The field per ampere of this coil is not known, so the drive is given in A only.
- **Δt:** from the file modification times. Median 9.25 s, mean 9.52 s; a ~90 s pause every 300 frames (saving).
  - The `t` value in `info.txt` does not match the file times.
  - The mean is used (`--fps 0.10507`).
  - Persistence lags are therefore in frames, ≈ 9.5 s each on average.
- **Scale 0.1175 µm/px: assumed**, the 20x value of Studies 2–7. The user does not know the magnification of this run. All lengths in µm scale with it; exponents and Λ do not.
- **Geometry:** the wall is horizontal; the bright domain grows upward. Used `--rotate 3 --bright-domain`.

```bash
python step1_extract_wall.py ~/MOKEdata_local/24jun16 --outdir P7925_0_04A_24jun16 --fps 0.10507 --um-per-px 0.1175 \
       --rotate 3 --bright-domain --offset auto --flatfield 2 --init-wall
python step2_persistence.py P7925_0_04A_24jun16 --row-mean --eps 0.5 1 1.5 2 3 5
python step7_roughness.py P7925_0_04A_24jun16 --every 20 --rmin 25
python step8_crossover.py P7925_0_04A_24jun16 P7925_0_04A_26jun16
python step10_universality.py P7925_0_04A_24jun16 P7925_0_04A_26jun16
```

## Wall extraction and checks

New step 1 options, introduced for these runs (README):
- `--bright-domain`: the growing domain is the bright one.
- The initial domain (left of the `--init-wall` path) can no longer switch. Tested on this run: it changes 0.7 % of the (frame, row) values by up to 10 px, and the advance from 148.1 to 147.2 px.

**Segmentation check.** `fig_qc_overlay.png` and `qc/overlay_small.mp4` (1 frame in 20) were approved by the user before the analysis.

- **Rows:** image rows 142–1143 are kept (1002 rows = 118 µm), the longest run where the initial wall is inside the aperture.
- **Drift:** dy = −0.2…+2.5 px, dx = −0.4…+0.3 px; scatter around the smoothed curve 0.04 / 0.02 px.
- **Motion:** all 6000 frames used. Total advance 147 px = 17.3 µm, **mean velocity 0.025 px/frame ≈ 0.30 nm/s**, no retreats > 1 px.
- **Noise:** σ_Δh = 0.47 px, so 3σ_Δh ≈ 1.5 px.
- **Tilt:** θ ≈ −14° to −13°; θ0 = −12.2° in step 7.
- **Overhangs:** 0.3 % of (row, frame), the fewest of all studies.

## Persistence, ξ(τ) and χ4(τ)

Row-mean subtracted. τ in frames (≈ 9.5 s).

| ε [px] | τ\* (max χ4) | χ4(τ\*) | ξ(τ\*) [µm] | ξ_max [µm] at τ |
|---|---|---|---|---|
| 0.5 | 134 (21 min) | 5.6 | 1.77 ± 0.04 | 3.35 at 3.3 h |
| 1 | 149 (24 min) | 7.6 | 2.03 ± 0.03 | 3.13 at 3.5 h |
| **1.5** | **217 (34 min)** | 8.1 | **2.34 ± 0.03** | 2.85 at 69 min |
| 2 | 242 (38 min) | 8.2 | 2.36 ± 0.03 | 2.83 at 75 min |
| 3 | 279 (44 min) | 8.5 | 2.42 ± 0.03 | 2.87 at 106 min |
| 5 | 319 (51 min) | 9.2 | 2.51 ± 0.04 | 3.04 at 89 min |

Above 3σ_Δh, ξ(τ\*) ≈ 2.3–2.5 µm and ξ_max ≈ 2.9–3.0 µm, close to Studies 2–4 (1.6–2.6 and 2.4–3.0 µm). τ\* ≈ 35–50 min.

## Roughness and l₀

`step7` (`fig_roughness.png`, `fig_local_width.png`): 994 rows on the base plane (117 µm), 300 frames, fit range 2.9–29 µm.

| estimator | ζ |
|---|---|
| global S(q), u_area | 1.00 ± 0.05 (thirds: 1.07, 1.03, 0.70) |
| local width, local rotation (contour) | 0.98 ± 0.03 (folds excluded 1.04, tilts > 45° excluded 1.04) |
| B(r) | 0.68 (biased low) |

- **Crossover l₀** (`step8`, `fig_crossover.png`): 22.7 µm from the local width (bootstrap 22–23), 14.1 µm from the global S(q), 7.6 µm from the local slope. l_fold = 5.2 µm (few folds), l_tilt = 35 µm. **l₀ ≈ 14–23 µm.**
- **Height distribution** (20 µm windows): skewness −0.16 ± 0.14, excess kurtosis +1.1 ± 0.4.

## qEW or qKPZ below l₀

`step10_universality.py` (`fig_universality.png`); comparison of all studies in [UNIVERSALITY.md](../UNIVERSALITY.md).

**A. Roughness below l₀** (ζ over [4 µm, l₀ = 22.7 µm]):

| estimator | ζ |
|---|---|
| contour | 1.00 (0.97–1.03) |
| folds excluded | 1.06 (1.03–1.08) |
| u_area | 1.01 |

- Broken power law, folds excluded: ζ₁ = 1.04 (0.99–1.11), l₀ = 44 µm.
- **χ² qEW / qKPZ / equilibrium = 412 / 2740 / 2091.**
- The value of ζ barely depends on ℓ_min = 3–5 µm.

**B. KPZ term:**
- Λ from the contour: 1.2–1.7 for L = 2–8 µm, 0.3–0.5 at L = 16 µm.
- Λ from h(x): 0.7–1.0, and −0.5…0 at L = 16 µm.
- The defect mask makes no difference (no static defects detected in this field).
- Λ does not depend on the lag (80–610 frames). **R ≤ 0.25.**
- v_n(φ) is flat within ±15 % for |φ| < 45°.

**Conclusion: qEW**, with at most the isotropic geometric term (Λ ≈ 1).

## Caveats

- **Scale assumed** (0.1175 µm/px). If it is wrong, all µm values, v and l₀ scale by the same factor; ζ, Λ and R do not change.
- **Frame timing:** Δt is not uniform (pauses every 300 frames). Lags are counted in frames, so durations are averages.
- **Not run yet:** steps 3–6 and 9.
