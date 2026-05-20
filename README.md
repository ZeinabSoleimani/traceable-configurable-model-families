# Constructing Traceable Configurable Process Model Families for Variant Comparison

This repository contains the artifact package for the paper **“Constructing Traceable Configurable Process Model Families for Variant Comparison.”**

The package contains the implementation, variant-labeled input logs, generated model artifacts, visualization outputs, and evaluation summaries for constructing configurable process model families from variant-labeled event logs. The method produces two linked artifacts:

- a **Configurable Directly-Follows Graph (C-DFG)** for relation-level evidence, availability, and dominance information; and
- a **Configurable Process Tree (CPT/CPST)** for executable variant-specific process model views.

Together, these artifacts support traceable comparison of shared and variant-specific behavior.

## Repository contents

The repository is organized as a lightweight artifact package. Some experiment outputs are provided as compressed archives to keep the repository uploadable.

```text
README.md
pyproject.toml
src.zip
evaluate_src.py
variant_labeled_event_logs.zip
Depression_example_experiment.zip
Antidepressants_experiment.zip
other_experiments_light_bundle.zip.zip
evaluation_results.zip
```

The files have the following roles:

```text
src.zip                              Source-code archive. Extracts to the Python package folder `src/`.
pyproject.toml                       Python packaging metadata.
evaluate_src.py                      Evaluation runner for the proposed method and baselines.
variant_labeled_event_logs.zip        Variant-labeled input logs used for the reproducible experiments.
Depression_example_experiment.zip     Generated artifacts for the motivating depression example.
Antidepressants_experiment.zip        Generated artifacts for the antidepressants experiment.
other_experiments_light_bundle.zip.zip Lightweight generated artifacts for the remaining public benchmark experiments.
evaluation_results.zip               Evaluation CSV files and metric summaries.
```

Before running the code, extract the source-code archive so that the repository root contains a `src/` folder:

```bash
unzip -q src.zip
```

Optionally extract the input logs and generated results:

```bash
unzip -q variant_labeled_event_logs.zip
unzip -q Depression_example_experiment.zip
unzip -q Antidepressants_experiment.zip
unzip -q other_experiments_light_bundle.zip.zip
unzip -q evaluation_results.zip
```

If possible, rename `other_experiments_light_bundle.zip.zip` to `other_experiments_light_bundle.zip` before public release to avoid confusion.

## Expected source-code layout after extraction

After extracting `src.zip`, the repository should contain:

```text
src/
  __init__.py
  cli.py
  cdfg.py
  data.py
  discovery.py
  enrichment.py
  io.py
  metrics.py
  ptree.py
  stats.py
  viz.py
```

The implementation includes:

- variant-labeled event-log loading;
- C-DFG construction and directly-follows statistics;
- configurable process-tree construction;
- anchor-localization and structure-first construction paths;
- variant-specific configured views; and
- DOT/SVG/PNG rendering utilities.

## Data availability and privacy

The repository contains only shareable artifact material for reproducibility.

The public benchmark inputs are provided as variant-labeled logs or can be regenerated from their public sources. The generated experiment artifacts include C-DFG tables, configurable-tree outputs, configured variant views, metadata, and summary files.

The raw clinical input data used for the antidepressants-related experiment is **not included**. Any included antidepressants-related files should be interpreted only as generated model artifacts or aggregate evaluation outputs, not as a release of the underlying clinical event log.

Before making the repository public or sharing it through an anonymized review service, check that the archives do not contain:

```text
raw clinical event logs
patient-level private data
local machine paths
usernames
institution-specific private paths
execution logs with identifying information
```

## Requirements

The core implementation is written in Python and is intended to run with the Python standard library.

Recommended environment:

```bash
python --version   # Python 3.10 or newer recommended
```

Install the package in editable mode from the repository root after extracting `src.zip`:

```bash
python -m pip install --upgrade pip
python -m pip install -e .
```

For rendering `.dot` files to `.svg` and `.png`, install Graphviz and make sure the `dot` executable is available:

```bash
dot -V
```

If Graphviz is not available, the code still writes `.dot`, `.csv`, `.json`, and `.txt` outputs. Use `--no-render` to skip image rendering.

## Quick validation

From the repository root, after extracting `src.zip`, run:

```bash
python -m compileall -q src
python -m src.cli example --out runs/depression_example --coverage-check --no-render
```

The example command creates a small built-in TRD/nTRD motivating example and writes C-DFG, configurable-tree, and configured-variant artifacts to:

```text
runs/depression_example/
```

If Graphviz is installed, omit `--no-render` to also produce SVG/PNG visualizations:

```bash
python -m src.cli example --out runs/depression_example --coverage-check
```

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
python -m src.cli run \
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
python -m src.cli run \
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

Common options for `python -m src.cli run` include:

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

## Generated construction artifacts

Each generated experiment folder may contain files such as:

```text
cdfg.csv                         Relation-level evidence table
cdfg.dot / cdfg.svg / cdfg.png   Full C-DFG visualizations
cdfg_core.*                      Compact C-DFG visualizations
configurable_tree.*              Configurable process-tree artifacts
configured_<variant>.*           Variant-specific configured process-tree views
shared_context.*                 Shared-context tree/visualization, when applicable
construction_summary.txt         Construction mode and residual-fragment summary
construction_metadata.json       Machine-readable construction metadata
```

The CPST visualizations encode structural availability through color. Relation-level dominance and support evidence remain available in the C-DFG outputs, especially in `cdfg.csv`, `cdfg.dot`, and compact C-DFG views.

## Evaluation

The evaluation runner compares the proposed configurable model family against baseline discovery strategies:

- `proposed_cpt`: proposed configurable process-tree family;
- `pooled_im`: one process tree discovered from the pooled log;
- `separate_im`: one process tree discovered independently per variant; and
- `separate_exact`: optional exact-log memorization baseline for small logs or sanity checks.

Example command for the built-in motivating example:

```bash
python evaluate_src.py \
  --example \
  --out evaluation/depression_example \
  --coverage-mode exact
```

Example command for external CSV logs:

```bash
python evaluate_src.py \
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
python evaluate_src.py \
  --inputs "path/to/logs/*.csv" \
  --format auto \
  --out evaluation/new_run_bounded \
  --coverage-mode bounded \
  --skip-df-proxy
```

The evaluation runner writes:

```text
summary_metrics.csv              Dataset-level and method-level summary metrics
per_variant_metrics.csv          Per-variant coverage and precision-proxy metrics
cdfg_metrics.csv                 C-DFG relation statistics
README_metrics.md                Local explanation of reported metrics
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

See the `README_metrics.md` files inside the extracted evaluation results for metric-level details.

## Reproducing the reported artifact outputs

A typical reproduction workflow is:

```bash
# 1. Extract the implementation and input logs
unzip -q src.zip
unzip -q variant_labeled_event_logs.zip

# 2. Install the package
python -m pip install -e .

# 3. Run a small built-in example
python -m src.cli example --out runs/depression_example --coverage-check --no-render

# 4. Run the method on a prepared trace-level or event-level CSV
python -m src.cli run \
  --input path/to/log.csv \
  --format auto \
  --out runs/my_log \
  --coverage-check \
  --no-render

# 5. Run evaluation
python evaluate_src.py \
  --inputs "path/to/logs/*.csv" \
  --format auto \
  --out evaluation/reproduced \
  --coverage-mode bounded \
  --skip-df-proxy
```

Use the actual paths created after extracting `variant_labeled_event_logs.zip`.

## Important consistency check before release

The repository currently uses `src` as the package folder. Therefore, all commands and imports should refer to `src`, for example:

```bash
python -m src.cli example --out runs/depression_example
```

The evaluation runner should also import from `src`, for example:

```python
from src.cdfg import CDFG, SINK, SOURCE, directly_follows_counts
from src.data import Case, Trace, VariantLog, depression_example_log
```

If `evaluate_src.py` still imports from `evidence_cpt`, update those imports before releasing the artifact.

## Citation

A final citation entry should be added after publication. Until then, cite the paper title as:

```text
Constructing Traceable Configurable Process Model Families for Variant Comparison
```

## License

The repository is released under the license specified in `pyproject.toml`.
