# CONOP Python Reproduction and Validation

Independent Python reproduction and cross-solution validation of selected CONOP9 penalty functions for quantitative stratigraphic correlation.

> Developed as part of the course "Frontiers and Applications of Big Data in Geology" (《地质学大数据前沿与应用》) at Nanjing University.

## 1. Overview

Quantitative stratigraphic correlation aims to align fossil first/last appearance events and marker horizons across sections into a single composite sequence. CONOP9 (Constrained Optimization) is a well-known implementation of this approach: it searches for the composite sequence that minimizes penalties defined from the observed section data.

This repository contains an independent Python implementation of the CONOP9 data model and penalty functions, built for the Seymour Island (Antarctica) Cretaceous–Paleogene ammonite dataset used in the course (12 sections, 49 taxa, 120 composite events). The implementation is validated against 21 archived CONOP9 solutions by fixing each archived sequence and recalculating the penalties from scratch — a direct, algorithm-level test of the penalty functions rather than an end-to-end optimizer comparison.

The project spans scientific validation (exact reproduction of Ordinal and Level penalties), algorithm engineering (an incremental ordinal cost path with numba JIT, ~15× faster than the full-recompute reference), and a documented remaining gap (Eventual penalty only partially reproduced).

## 2. Validation status

| Penalty  | Test solutions | Exact matches |   MAE | Status     |
| -------- | -------------: | ------------: | ----: | ---------- |
| Ordinal  |             21 |            21 |     0 | Reproduced |
| Level    |             21 |            21 |     0 | Reproduced |
| Eventual |             21 |             0 | 33.76 | Partial    |

For Eventual, the mean signed error is **-33.76** (Python consistently lower than CONOP9) and the maximum absolute error is **83**. Eventual is **not fully reproduced**: the counting rule still differs from CONOP9 in a systematic way, and no archived solution matched exactly.

The full per-run comparison is in `results_py/cross_solution_validation/` (`penalty_validation.csv`), and the analysis is described in the [technical report](#8-technical-report).

## 3. Why cross-solution validation matters

Checking a single solution can produce an accidental single-point agreement: one permutation happens to give the same penalty even if the implementation differs. By validating against 21 archived CONOP9 solutions — spanning different annealing parameter groups and distinct optimized permutations — the penalty functions are tested across a range of sequences. For each solution the sequence is fixed and the penalty is independently recalculated in Python, so the comparison isolates the penalty implementation itself from optimization randomness.

## 4. Architecture

```text
CONOP-run reference data
        │
        ▼
conop_py/io.py
        │
        ▼
cost.py
 ├─ Ordinal
 ├─ Level
 └─ Eventual
        │
        ▼
anneal.py / incremental.py
        │
        ▼
scripts/conop.py
        │
        ▼
validation / results / plots
```

## 5. Python implementation

- **Parser / validation** (`conop_py/io.py`) — reads the CONOP input formats (`loadfile.dat`, `events.txt`, `sections.txt`, `bestsoln.dat`) and validates dataset schema.
- **Penalty functions** (`conop_py/cost.py`) — Ordinal (inversion count), Level (L1 isotonic regression with PAV + box constraints), Eventual (forcing-event weighted Level), coexistence checks. This is the full-recompute reference implementation used by the regression suite.
- **Simulated annealing** (`conop_py/anneal.py`) — Metropolis acceptance with geometric cooling, configurable temperature/ratio/steps, optional early stopping and coexistence penalty.
- **Incremental ordinal optimization** (`conop_py/incremental.py`) — `FastOrdinalState` maintains per-section ordinal caches and computes move deltas with a differential formula (O(n_s) per move), with a numba JIT inner loop; it falls back to pure Python automatically when numba is unavailable.
- **Multistart** (`scripts/run_multistart.py`) — parallel multi-restart search via `multiprocessing.Pool`.
- **Plotting** (`conop_py/plotting.py`) — shared matplotlib setup and trajectory parsing.

## 6. Quick start

```bash
# Dataset schema validation
uv run --with-requirements requirements.txt \
  python scripts/conop.py validate

# Evaluate one solution (ordinal / level / eventual)
uv run --with-requirements requirements.txt \
  python scripts/conop.py eval CONOP-run/bestsoln.dat

# Regression suite (scientific baseline)
uv run --with-requirements requirements.txt \
  python -m pytest tests/test_regression.py -v
```

## 7. Reproducibility

- The 21 **archived CONOP9 solutions** in `results/` (7 parameter groups × 3 runs) are the reference input set. Each `run_N/` directory contains the CONOP9 `bestsoln.dat` and the corresponding `outmain.txt` with the reference penalties.
- Python **independently recalculates** Ordinal / Level / Eventual for each fixed archived sequence (`scripts/validate_conop9_penalties.py`); no optimizer is involved in the comparison.
- Optimization randomness (simulated annealing) is thus separated from penalty validation.
- Figures and summary tables under `results_py/cross_solution_validation/` are generated from the same validation run.

Note on the reference software: the archived solutions were produced with the CONOP9 executable (P. Sadler's freeware program; see Sadler & Cooper 2003). Reference outputs are included in this repository so that the validation can be reproduced without re-running the original binary.

## 8. Technical report

- [conop9-penalty-cross-solution-validation.pdf](docs/paper/conop9-penalty-cross-solution-validation.pdf) — *CONOP9 代价函数的跨解复现与 Eventual 偏差诊断* ("Cross-solution reproduction of CONOP9 penalty functions and diagnosis of the Eventual bias"), by Xinyu Du.

This is a **course paper / technical report** written for the course; it is not a journal publication. The report presents the validation methodology, the 21-solution results above, and the Eventual bias analysis.

## 9. Project structure

```text
conop_py/            Python implementation (io, cost, anneal, incremental, plotting)
scripts/             CLI (conop.py) and analysis / validation scripts
tests/               Regression tests (tests/test_regression.py)
CONOP-run/           CONOP9 input data and reference run outputs
results/             21 archived CONOP9 runs (reference solutions + outmain penalties)
results_py/          Python-generated results (cross-solution validation, diagnostics)
docs/paper/          Course paper (PDF)
references/          Reference literature PDFs
```

## 10. Limitations

- **Eventual is not fully reproduced**: on the 21 archived solutions, Python's Eventual is systematically lower than CONOP9 (mean signed error -33.76, max absolute error 83). This is a known, documented gap.
- The validation solutions all come from the same dataset (Seymour Island / Quiriquina ammonite data, 12 sections).
- The tested solutions mainly cover the CONOP9 optimized / high-fit solution region; the results do not establish universal equivalence of the implementations over all permutations or other datasets.
- Ordinal and Level are reproduced exactly on the tested set only; they are not claimed to be bit-for-bit equivalents of CONOP9 beyond the evaluated configurations.

## 11. Course context / Attribution

Project developed for the graduate course "Frontiers and Applications of Big Data in Geology" (《地质学大数据前沿与应用》), School of Earth Sciences and Engineering, Nanjing University (Spring 2026, Prof. Junxuan Fan). The Seymour Island ammonite dataset and the CONOP workflow were provided as course materials.

Key reference: Sadler, P. M., & Cooper, R. A. (2003). Best-fit intervals and consensus sequences. In *High-Resolution Approaches in Stratigraphic Paleontology* (pp. 49–94). Kluwer Academic Publishers.

The group presentation for the course was completed; the final presentation is not currently distributed in this repository.
