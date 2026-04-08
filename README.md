# Reusable Survival Analysis Pipeline

This repository is a reusable extraction of an older R-based survival-analysis workflow. The target architecture follows the same broad shape as NIDAP-Community/multi-gene-correlations, but the repository is not fully standardized yet.

The current state is mixed:

- the main reusable analysis code lives in src/core/
- provenance guardrails live in src/utilities/
- newer pooled and piecewise drivers live in drivers/
- older HTML-oriented drivers still live in src/drivers/
- the root entrypoint now dispatches maintained workflows
- test_data/ now contains synthetic smoke-test placeholders
- renv/ exists, but renv.lock has not been checked in yet

This README documents the current repository honestly so cleanup and refactoring can proceed from a clear baseline.

## Goals

- keep the analysis layer reusable across expression plus survival datasets
- preserve non-interactive Rscript execution
- preserve run metadata and provenance guardrails
- keep optional reporting separate from core analysis code
- minimize hard-coded cohort, path, and column assumptions

## Current Layout

```text
.github/
container/
docs/
drivers/
entrypoint.R
renv/
results/
scripts/
src/
  core/
  drivers/
  utilities/
test_data/
CHANGELOG.md
LICENSE
README.md
```

## Architecture Status

The intended architecture is only partially in place.

- root entrypoint: present and dispatches maintained workflows
- container/: present but currently empty
- scripts/: present but currently empty
- src/core/: present and contains the main reusable survival-analysis function
- src/drivers/: present, but currently contains older cohort-specific and report-oriented drivers
- src/utilities/: present and currently contains provenance and output guardrails
- test_data/: present and contains small synthetic smoke-test placeholders
- results/: present and suitable for run-scoped outputs
- docs/: present, with architecture documentation started
- renv.lock: currently missing

## Code Layers

### Core analysis layer

The reusable analysis engine is in src/core/GeneExpressionBasedSurvivalAnalys.r. It provides the main analyze_survival_expression() function used by the older driver layer. This is the strongest candidate for the long-term reusable core.

### Utility layer

src/utilities/run_provenance_guardrails.R provides output-path guardrails and run metadata helpers. It is already aligned with the repository goal of preserving provenance and should remain a shared dependency.

### Maintained driver layer

The drivers/ directory contains newer command-line scripts for pooled and piecewise modeling:

- run_piecewise_dataset_interaction.R
- run_piecewise_subtype_specific.R
- run_single_cohort_survival.R
- run_single_cohort_reporting_bundle.R
- run_signature_consistency_summary.R
- run_finalized_pooled_results_bundle.R

These scripts are non-interactive and executable with Rscript. The pooled drivers are still tied to TCGA and METABRIC defaults, while the single-cohort analysis and reporting drivers are intended for one dataset at a time.

They now share a common execution path through entrypoint.R and default to run-scoped output directories under results/.

### Legacy and reporting-oriented layer

The src/drivers/ directory contains older drivers that combine analysis with HTML and markdown reporting. They still rely on the reusable core, but they are closer to manuscript-oriented workflow code than to a generalized maintained driver layer.

Until the driver layout is unified, treat src/drivers/ as legacy-style workflow wrappers and drivers/ as the more active extraction layer.

## Input Contracts

Most workflows expect two inputs.

The single-cohort maintained driver expects one expression matrix plus one metadata table.

### Expression matrix

Expected shape:

- genes in rows
- samples in columns
- first column contains gene identifiers

Typical example:

```text
Gene    Sample1    Sample2    Sample3
SAT1    1.2        0.8        1.5
SRM     2.1        1.9        2.4
```

### Metadata table

Expected shape:

- one row per sample
- one sample identifier column
- event and time columns for each endpoint used in the run

Commonly required fields in the current drivers:

- sample identifier column
- subtype column
- OS and OS_time
- PFI and PFI_time

Commonly used covariates in the current drivers:

- age column
- stage column

## Current Driver Assumptions

The repository is intended to become reusable, but several current scripts still encode historical assumptions.

### Dataset and cohort assumptions

- newer pooled drivers currently default to TCGA plus METABRIC inputs
- older src/drivers/ scripts are BRCA-focused and assume BRCA subtype labels
- several scripts assume a TCGA primary-tumor filter based on sample IDs containing _01

### Column-name assumptions

The reusable core supports configurable column names, but the current drivers frequently default to:

- sample column: sampnames
- gene column: Gene
- subtype column: Subtype
- TCGA stage column: ajcc_pathologic_tumor_stage
- METABRIC stage column: stage_collapsed
- age column: Diagnosis_Age
- endpoints: OS and PFI with paired time columns OS_time and PFI_time

For the maintained single-cohort driver, the subtype column is optional unless subtype filtering or subtype-adjusted modeling is requested.

### Feature assumptions

Several drivers hard-code the historical four-gene signature:

- SAT1
- SRM
- ODC1
- AMD1

Several scripts also hard-code a 60-month early versus late split for piecewise modeling.

## Workflow Selection

Use the maintained workflows according to the analysis shape you need.

- Use `single-cohort-survival` when you have one dataset and want to run the reusable core analysis engine directly.
- Use `single-cohort-reporting-bundle` when you already have a completed single-cohort run directory and want an overall forest summary, plus an optional subgroup heatmap if you supply subgroup rows explicitly.
- Use `piecewise-dataset-interaction` when you have two datasets and want pooled piecewise modeling with dataset interaction terms.
- Use `piecewise-subtype-specific` when you have two datasets and want pooled piecewise modeling within subtype strata. This workflow also emits a standardized subgroup summary table for downstream heatmaps.
- Use `signature-consistency-summary` when upstream pooled and subtype result tables already exist and you want a compact comparison summary.
- Use `finalized-pooled-results-bundle` when upstream pooled two-dataset result tables already exist and you want pooled forest plots, subtype heatmaps, and manuscript-style summary artifacts.

## Execution

All workflows should run non-interactively via Rscript.

The canonical CLI surface is now the root entrypoint.

List workflows or show usage:

```bash
Rscript entrypoint.R --help
```

Run a maintained workflow:

```bash
Rscript entrypoint.R piecewise-dataset-interaction \
  --output-dir results/example_piecewise_dataset_interaction
```

```bash
Rscript entrypoint.R piecewise-subtype-specific \
  --subtype-list BRCA_Basal \
  --adjust-covariates none
```

```bash
Rscript entrypoint.R single-cohort-survival \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata.csv \
  --endpoints OS \
  --features Signature
```

```bash
Rscript entrypoint.R single-cohort-reporting-bundle \
  --input-run-dir results/smoke_test_single_cohort \
  --subgroup-summary-csv test_data/single_cohort_subgroup_summary.csv
```

Downstream summary and reporting workflows require explicit upstream result tables:

The shared subgroup-summary contract for downstream heatmaps uses these columns:

- subgroup
- subgroup_type
- feature
- endpoint
- comparison
- HR
- CI_low
- CI_high
- p_value
- p_value_bh

Additional columns such as dataset, window, model, or covariates may be present and are preserved when available.

```bash
Rscript entrypoint.R signature-consistency-summary \
  --overall-implied-csv <path> \
  --overall-summary-csv <path> \
  --subtype-implied-csv <path> \
  --subtype-summary-csv <path> \
  --subtype-terms-csv <path>
```

```bash
Rscript entrypoint.R finalized-pooled-results-bundle \
  --overall-signature-table <path> \
  --pooled-term-table <path> \
  --subtype-all-features-table <path>
```

The maintained driver scripts are also still directly runnable if needed:

```bash
Rscript drivers/run_piecewise_dataset_interaction.R \
  --output-dir results/example_piecewise_dataset_interaction
```

```bash
Rscript drivers/run_piecewise_subtype_specific.R \
  --output-dir results/example_piecewise_subtype_specific
```

```bash
Rscript drivers/run_signature_consistency_summary.R \
  --output-dir results/example_signature_consistency
```

```bash
Rscript drivers/run_single_cohort_survival.R \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata.csv \
  --endpoints OS
```

```bash
Rscript drivers/run_single_cohort_reporting_bundle.R \
  --input-run-dir results/smoke_test_single_cohort \
  --subgroup-summary-csv test_data/single_cohort_subgroup_summary.csv
```

```bash
Rscript drivers/run_finalized_pooled_results_bundle.R \
  --output-dir results/example_pooled_bundle
```

The src/drivers/ scripts are still runnable, but they should be treated as legacy-style wrappers until the architecture is cleaned up.

## Smoke-Test Data

The repository now includes synthetic smoke-test placeholders under test_data/.

These files are intended for interface validation, not scientific use. They are useful for checking:

- input-column contracts
- entrypoint dispatch
- output directory creation
- provenance metadata generation

See test_data/README.md for example commands.

## Dependencies

The repository currently uses R packages including:

- dplyr
- readr
- survival
- ggplot2
- survminer
- stringr
- glmnet

Package management is not fully standardized yet because renv.lock is currently missing. Restoring a fully reproducible local environment should be treated as an immediate follow-up task.

## Outputs and Provenance

The intended output convention is a run-scoped directory under results/.

Where the provenance utility is used, runs can write:

- RUN_METADATA.md
- logs/sessionInfo.txt
- tables/
- figures/

Current output behavior is still inconsistent:

- drivers/ generally write under results/
- older src/drivers/ scripts still default to outputs/ legacy paths unless explicitly overridden

Normalizing all maintained workflows to results/<topic>/<run_id>/ is part of the near-term cleanup plan.

## Smallest Safe Cleanup Plan

Before large refactors, the smallest safe cleanup is:

1. document src/drivers/ as legacy or reporting-oriented until they are moved or replaced
2. add renv.lock and formalize environment restoration
3. decide whether maintained workflow scripts should remain in drivers/ or move under src/drivers/
4. tighten shared CLI/help behavior across maintained drivers
5. begin second-wave reduction of TCGA and METABRIC assumptions

## Immediate Gaps

The following gaps are known and intentional to fix next:

- scripts/ is empty
- container/ is empty
- renv.lock is missing

## Working Principle

Prefer the implementation that is more reusable, more explicit, and more reproducible, even when it requires removing historical assumptions from the original workflow.