# MAG — Multiverse Alpha Gaze

Analysis code for a **multiverse** study on the relationship between **posterior alpha-band power** and **eye-movement–derived gaze metrics** during working memory. The pipeline operates on data from the Alpha Oculomotor Control (AOC) study ([registered report](https://www.authorea.com/users/819882/articles/1218693-modulations-of-posterior-alpha-power-during-working-memory-co-vary-with-task-dependent-eye-movement-patterns); *Psychophysiology*): simultaneous EEG and eye-tracking in **Sternberg** and **N-back** tasks. This secondary analysis was **not preregistered**; all specifications, code, and outputs are nevertheless kept **open and version-controlled**.

---

## Scientific aim

The goal is to quantify how robust conclusions are when theoretically defensible preprocessing and analysis choices vary jointly. A seven-dimensional decision grid defines **1,440 “universes” per task**. Each universe combines choices over:

| Dimension | Options |
|-----------|---------|
| Electrode region | posterior, occipital |
| Spectral 1/f handling | FOOOF-parameterized (“FOOOFed”), non-parameterized |
| Retention latency window | 0–500, 0–1000, 0–2000, 1000–2000 ms (the shortest window is excluded in primary R analyses) |
| Alpha band | canonical 8–14 Hz, individual alpha frequency (IAF) |
| Gaze measure | scan path length, gaze velocity, microsaccades, BCEA, gaze deviation |
| EEG baseline | raw, dB, percent change (pre-stimulus −0.5 to −0.25 s) |
| Gaze baseline | raw, dB, percent change (same pre-stimulus window) |

Trial-level linear mixed models (LMMs) relate alpha to gaze and task condition, including EEG-only and gaze-only model families; the R implementation uses the [`multiverse`](https://cran.r-project.org/package=multiverse) package (Sarma et al., 2021). Optional quality control joins trial-level FOOOF *R*² and can filter low-*R*² FOOOF fits.

---

## Repository layout

| Path | Role |
|------|------|
| `MAG_multiverse_prep.m` | Trial-level CSV export: one row per trial × universe (Sternberg & N-back). |
| `subject-level/MAG_multiverse_prep_subject.m` | Subject-level export: trial-averaged spectra and gaze summaries (smaller tables, comparable to classical subject-level plots). |
| `MAG_multiverse_*_analysis.R` | Fits multiverse LMMs and writes result CSVs (Stenberg / N-back). |
| `MAG_multiverse_*_visualize.R` | Figures from multiverse outputs. |
| `subject-level/MAG_multiverse_*` | Subject-level analysis and visualization scripts. |
| `MAG_multiverse_results_overview.R` | Aggregates key counts from result CSVs into a summary table. |
| `MAG_multiverse_decision_tree.py` | Matplotlib figure: “traditional” single specification vs full multiverse (trunk layout). |

Language mix reflects the pipeline: **MATLAB** for feature assembly from the existing AOC feature store, **R** for statistics and most figures, **Python** for the decision-tree schematic.

---

## Dependencies

- **MATLAB**: AOC project `startup` / `setup` on the path; precomputed EEG and gaze features under the AOC directory layout (see path blocks in the prep scripts).
- **R**: `lme4`, `lmerTest`, `broom.mixed`, `tidyverse`, `multiverse`, and standard tidyverse dependencies.
- **Python** (decision tree only): `matplotlib`, `numpy`.

---

## Paths and environment variables

Scripts default to institute network paths for input (AOC) and output (MAG). Override as needed:

- `AOC_MULTIVERSE_DIR` — directory containing multiverse CSVs and result tables (R).
- `AOC_R2_DIR` — optional directory with FOOOF *R*² control CSVs for trial-level filtering (R).
- `AOC_VERBOSE_LOGGING` — set to `true` for full `lmer` diagnostics (R).

The Python script contains fixed output paths for GitHub vs local figure export; edit `OUT_GH` / `OUT_SC` before running elsewhere.

---

## Outputs

Prep scripts write large trial-level (and smaller subject-level) tables under `data/multiverse/` (location configurable). Analysis scripts produce CSVs of model terms, simple effects, and condition-wise summaries; filenames are prefixed with `MAG_multiverse_`. The results overview script collates significance patterns for manuscript cross-checks.

---

## Status

Manuscript **in preparation**. The repository may evolve until the final analysis freeze.

---

## Citation

If the code is reused, cite the AOC registered report and this repository. A dedicated MAG citation will be added when the paper is available.
