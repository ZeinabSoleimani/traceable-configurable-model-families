# Constructing Traceable Configurable Process Model Families for Variant Comparison

This repository contains the artifact package for the paper **“Constructing Traceable Configurable Process Model Families for Variant Comparison.”**

The package includes the implementation, generated model artifacts, visualization outputs, evaluation summaries, and validation scripts for constructing configurable process model families from variant-labeled event logs. The method produces two linked artifacts:

- a **Configurable Directly-Follows Graph (C-DFG)** for relation-level evidence, availability, and dominance information; and
- a **Configurable Process Tree (CPT/CPST)** for executable variant-specific process model views.

Together, these artifacts support traceable comparison of shared and variant-specific behavior.

## Repository structure

```text
evidence_cpt/                    Core implementation
exp6/                            Generated artifacts for the experiment logs
evaluation/                      Evaluation summaries and metric CSV files
evaluate_evidence_cpt.py         Evaluation runner
evaluate_evidence_cpt_v2.py      Updated evaluation runner used for the current package
paper_alignment_smoke_tests.py   Smoke tests for paper-alignment checks
PATCH_NOTES.md                   Summary of implementation patches and validation notes
ER.pdf                           Paper draft/reference copy, if included in the release package
```

Each experiment folder under `exp6/` contains generated artifacts such as:

```text
cdfg.csv                         Relation-level evidence table
cdfg.dot / cdfg.svg / cdfg.png   Full C-DFG visualizations
cdfg_core.*                      Compact C-DFG visualizations
configurable_tree.*              Configurable process tree artifacts
configured_<variant>.*           Variant-specific configured process-tree views
shared_context.*                 Shared-context tree/visualization, when applicable
construction_summary.txt         Construction mode and residual-fragment summary
construction_metadata.json       Machine-readable construction metadata
```

Each folder under `evaluation/` contains:

```text
summary_metrics.csv              Dataset-level and method-level summary metrics
per_variant_metrics.csv          Per-variant coverage and precision-proxy metrics
cdfg_metrics.csv                 C-DFG relation statistics
README_metrics.md                Local explanation of reported metrics
```

## Data availability

The repository includes the non-sensitive logs/artifacts needed for the reproducible public and synthetic experiments.

The raw input log for `vad_drug_events_category_named` is **not included** because it is clinical identifiable data and cannot be shared. Where present, this repository only includes generated artifacts and aggregate evaluation outputs for `vad_drug_events_category_named`; these files should not be interpreted as releasing the underlying clinical event log.

## Requirements

The implementation is written in Python and uses only the Python standard library for the core construction and evaluation code.

Recommended environment:

```bash
python --version   # Python 3.10 or newer recommended
```

For rendering `.dot` files to `.svg` and `.png`, install Graphviz and make sure the `dot` executable is available:

```bash
dot -V
```

If Graphviz is not available, the code still writes `.dot`, `.csv`, `.json`, and `.txt` outputs. You can also use `--no-render` to skip image generation.

## Quick validation

From the repository root, run:

```bash
python -m compileall -q evidence_cpt
python paper_alignment_smoke_tests.py
python -m evidence_cpt.cli example --out runs/depression_example --coverage-check
```

The smoke tests should report:

```text
All paper-alignment smoke tests passed.
```

The example command creates a small built-in TRD/nTRD motivating example and writes C-DFG, configurable-tree, and configured-variant artifacts to `runs/depression_example/`.

## Running the construction on a trace-level CSV

A trace-level CSV should contain one row per case with the default columns:

```text
case_id,variant,trace
```

The `trace` column should contain a comma-separated activity sequence, for example:

```text
P01,TRD,"A,B,C,D,E"
```

Run:

```bash
python -m evidence_cpt.cli run \
  --input path/to/log.csv \
  --format trace \
  --out runs/my_log \
  --coverage-check
```

Custom column names can be supplied with:

```bash
--case-col CASE_COLUMN \
--variant-col VARIANT_COLUMN \
--trace-col TRACE_COLUMN \
--trace-sep ","
```

## Running the construction on an event-level CSV

An event-level CSV should contain one row per event with the default columns:

```text
case_id,variant,activity,timestamp
```

Run:

```bash
python -m evidence_cpt.cli run \
  --input path/to/events.csv \
  --format event \
  --out runs/my_event_log \
  --case-col case_id \
  --variant-col variant \
  --activity-col activity \
  --timestamp-col timestamp \
  --coverage-check
```

The timestamp column is used to order events within each case. If no timestamp column is available, events are ordered by case order in the file.

## Main construction parameters

Common options for `python -m evidence_cpt.cli run` include:

```text
--theta FLOAT                         Minimum directly-follows support threshold
--alpha FLOAT                         Significance level for pairwise support tests
--n-min INT                           Minimum retained relation count
--min-shared-activities INT           Minimum number of shared activities for shared-context construction
--coverage-check                      Run finite-language observed-trace coverage check
--no-render                           Skip PNG/SVG rendering
--cdfg-core-min-support FLOAT         Minimum max support for compact C-DFG visualization
--cdfg-core-min-support-delta FLOAT   Minimum support difference for compact C-DFG visualization
--cdfg-core-dominance-only            Keep only dominant/differential edges in compact C-DFG view
--cpst-show-tau                       Show trivial tau leaves in default CPST visualizations
```

## Evaluation

The evaluation runner compares the proposed configurable model family against baseline discovery strategies:

- `proposed_cpt`: proposed configurable process-tree family;
- `pooled_im`: one process tree discovered from the pooled log;
- `separate_im`: one process tree discovered independently per variant; and
- `separate_exact`: optional exact-log memorization baseline for small logs or sanity checks.

Example command for the built-in motivating example:

```bash
python evaluate_evidence_cpt_v2.py \
  --example \
  --out evaluation/depression_exact \
  --coverage-mode exact
```

Example command for external CSV logs:

```bash
python evaluate_evidence_cpt_v2.py \
  --inputs "path/to/logs/*.csv" \
  --format auto \
  --out evaluation/new_run \
  --coverage-mode bounded \
  --max-loop-depth 3 \
  --lang-loop-depth 2 \
  --max-lang-traces 5000
```

For large loop-heavy logs, bounded coverage or skipped directly-follows proxy evaluation may be more practical:

```bash
python evaluate_evidence_cpt_v2.py \
  --inputs "path/to/logs/*.csv" \
  --format auto \
  --out evaluation/new_run_bounded \
  --coverage-mode bounded \
  --skip-df-proxy
```

## Reported metrics

The evaluation files report:

- observed-trace coverage on train/test splits;
- bounded directly-follows precision/recall proxies, when enabled;
- model-family size and simplicity metrics;
- compactness gain compared with separate per-variant models;
- shared and variant-specific node ratios for the configurable tree;
- C-DFG shared, variant-specific, dominant, and balanced relation counts; and
- runtime.

See the `README_metrics.md` files inside `evaluation/*/` for metric-level details.

## Notes on generated visualizations

The CPST visualizations encode structural availability through color. Relation-level dominance and support evidence remain available in the C-DFG outputs, especially in `cdfg.csv`, `cdfg.dot`, `cdfg_with_supports.dot`, and compact C-DFG views.

The current patch removes activity-level `ctx-dom` labels from CPST visualizations to keep relation-level statistical evidence in the C-DFG layer and structural availability in the CPST layer.

## Reusing the package

To apply the package to a new variant-labeled event log:

1. prepare the log as either a trace-level or event-level CSV;
2. run `python -m evidence_cpt.cli run` with the correct `--format`;
3. inspect `cdfg.csv` for relation-level evidence;
4. inspect `configurable_tree.*` for the family-level configurable model; and
5. inspect `configured_<variant>.*` for per-variant configured views.

## Citation

A final citation entry should be added after publication. Until then, cite the paper title as:

```text
Constructing Traceable Configurable Process Model Families for Variant Comparison
```

## License

Add the intended license before public release.
