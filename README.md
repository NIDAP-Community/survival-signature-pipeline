# Reusable Survival Analysis Pipeline

This repository is a reusable extraction of an older R-based survival-analysis workflow. The architecture target is a reusable, non-interactive survival-analysis pipeline, but the repository is not fully standardized yet.

The current state is mixed:

- the main reusable analysis code lives in src/core/
- shared provenance utilities live in src/utilities/
- shared reporting utilities also now live in src/utilities/
- shared maintained-driver CLI utilities now also live in src/utilities/
- newer pooled and piecewise drivers live in drivers/
- older HTML-oriented drivers still live in src/drivers/
- the root entrypoint now dispatches maintained workflows
- test_data/ now includes one canonical single-cohort smoke-test dataset plus additional synthetic fixtures
- renv/ is active and renv.lock is checked in

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
- drivers/: present and contains the maintained analysis and reporting workflows
- src/drivers/: present, but currently contains older cohort-specific and report-oriented wrappers
- src/utilities/: present and now contains both provenance helpers and reporting-side helpers
- test_data/: present and contains one canonical smoke-test input pair for the maintained single-cohort workflow plus additional synthetic fixtures
- results/: present and suitable for run-scoped outputs
- docs/: present, with architecture documentation started
- renv.lock: present

## Code Layers

### Core analysis layer

The reusable analysis engine is in src/core/GeneExpressionBasedSurvivalAnalys.r. It provides the main analyze_survival_expression() function used by the older driver layer. This is the strongest candidate for the long-term reusable core.

### Provenance utility layer

src/utilities/run_provenance_guardrails.R provides output-path guardrails, run metadata helpers, and the shared subgroup-summary normalization and validation helpers used by reporting workflows.

### Reporting utility layer

src/utilities/reporting_helpers.R provides shared reporting-side helpers for subgroup heatmap preparation, shared reporting output scaffolding, shared reporting table-formatting helpers, and common plotting defaults. These helpers reduce duplicated downstream reporting code, but they do not merge the reporting workflows or change the subgroup-summary contract.

### Maintained CLI utility layer

src/utilities/cli_helpers.R provides the small shared utility layer used by the maintained drivers before workflow-specific preparation begins. It currently centralizes:

- standardized help rendering and argument parsing
- required-argument enforcement and workflow-aware CLI error messages
- compatibility alias handling for maintained flags
- comma-separated list parsing and repeated enum-like flag validation
- maintained OS/PFI endpoint mapping
- minimal required input-file existence checks for maintained analysis drivers
- maintained analysis-table loading for expression and metadata inputs

It does not centralize dataset-specific filtering, metadata validation, endpoint-column checks inside each workflow, cohort assumptions, or modeling logic.

### Maintained driver layer

The drivers/ directory contains newer command-line scripts for pooled and piecewise modeling:

- run_piecewise_dataset_interaction.R
- run_piecewise_subtype_specific.R
- run_single_cohort_survival.R
- run_single_cohort_reporting_bundle.R
- run_signature_consistency_summary.R
- run_finalized_pooled_results_bundle.R

These scripts are non-interactive and executable with Rscript. The pooled drivers currently assume two harmonized datasets and now expose neutral `dataset1_*` and `dataset2_*` flags, while retaining the historical TCGA/METABRIC-shaped names as compatibility aliases. They also accept optional `dataset1_label` and `dataset2_label` flags so emitted pooled outputs can use dataset names that match the supplied inputs. The maintained drivers expect datasets to be post-filtered upstream rather than applying a primary-only sample filter internally. The single-cohort analysis and reporting drivers are intended for one dataset at a time.

They now share a common execution path through entrypoint.R and default to run-scoped output directories under results/.

They also now share a small maintained-driver utility layer in src/utilities/cli_helpers.R for repeated CLI, endpoint, file-check, and loader mechanics. The workflows themselves still own their cohort-specific semantics.

The maintained reporting workflows still remain separate:

- run_single_cohort_reporting_bundle.R handles single-cohort downstream summaries and optional subgroup heatmaps
- run_finalized_pooled_results_bundle.R handles pooled two-dataset downstream forest plots, subtype heatmaps, and manuscript-style reporting outputs

### Legacy layer

The src/drivers/ directory contains older drivers that combine analysis with HTML and markdown reporting. They still rely on the reusable core, but they are closer to manuscript-oriented workflow code than to the maintained drivers/ layer.

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
GENE_A  1.2        0.8        1.5
GENE_B  2.1        1.9        2.4
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

- newer pooled drivers currently default to two example datasets represented by TCGA-like and METABRIC-like smoke fixtures
- older src/drivers/ scripts are BRCA-focused and assume BRCA subtype labels
- maintained drivers expect post-filtered inputs and do not apply a primary-only sample filter internally
- several legacy scripts still assume a TCGA primary-tumor filter based on sample IDs containing _01

These cohort-aware assumptions now mainly appear in the compatibility alias layer and in legacy wrappers, but the maintained pooled CLI now prefers neutral names such as dataset1_expr_file, dataset2_expr_file, dataset1_stage_column, and dataset2_stage_column. The pooled model outputs can label datasets generically through dataset1_label and dataset2_label, even though the workflows remain semantically two-dataset workflows.

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

Several drivers hard-code a placeholder four-feature signature:

- GENE_A
- GENE_B
- GENE_C
- GENE_D

Several scripts also hard-code a 60-month early versus late split for piecewise modeling.

## Workflow Selection

Use the maintained workflows according to the analysis shape you need.

- Use `single-cohort-survival` when you have one dataset and want to run the reusable core analysis engine directly.
- Use `single-cohort-reporting-bundle` when you already have a completed single-cohort run directory and want an overall forest summary, plus an optional subgroup heatmap if you supply subgroup rows explicitly.
- Use `piecewise-dataset-interaction` when you have two harmonized datasets and want pooled piecewise modeling with dataset interaction terms.
- Use `piecewise-subtype-specific` when you have two harmonized datasets and want pooled piecewise modeling within subtype strata. This workflow also emits a standardized subgroup summary table for downstream heatmaps.
- Use `signature-consistency-summary` when upstream pooled and subtype result tables already exist and you want a compact comparison summary.
- Use `finalized-pooled-results-bundle` when upstream pooled two-dataset result tables already exist and you want pooled forest plots, subtype heatmaps, and manuscript-style summary artifacts. This reporting workflow now consumes the shared subgroup-summary contract.

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
  --subgroup-summary-table test_data/single_cohort_subgroup_summary.csv
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
  --subgroup-summary-table <path>
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
  --overall-implied-csv <path> \
  --overall-summary-csv <path> \
  --subtype-implied-csv <path> \
  --subtype-summary-csv <path> \
  --subtype-terms-csv <path> \
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
  --subgroup-summary-table test_data/single_cohort_subgroup_summary.csv
```

```bash
Rscript drivers/run_finalized_pooled_results_bundle.R \
  --overall-signature-table <path> \
  --pooled-term-table <path> \
  --subgroup-summary-table <path> \
  --output-dir results/example_pooled_bundle
```

The src/drivers/ scripts are still runnable, but they should be treated as legacy-style wrappers until the architecture is cleaned up.

## Smoke-Test Data

The repository now includes a canonical single-cohort smoke-test dataset under test_data/, plus additional synthetic fixtures for other interface checks.

The two current validation anchors in this repository are:

- the canonical single-cohort smoke baseline at results/smoke_test_single_cohort_canonical
- the pooled reporting contract smoke baseline at results/pooled_reporting_contract_smoke_shared_plotting

After the maintained-driver CLI and utility cleanup work, both baselines were rerun successfully through the root entrypoint during these checkpoint paths as well:

- results/smoke_test_single_cohort_canonical_cli_help
- results/smoke_test_single_cohort_centralized_parsing
- results/smoke_test_single_cohort_endpoint_helpers
- results/smoke_test_single_cohort_loader_helpers
- results/pooled_reporting_contract_smoke_cli_help
- results/pooled_reporting_contract_smoke_centralized_parsing
- results/pooled_reporting_contract_smoke_endpoint_helpers
- results/pooled_reporting_contract_smoke_loader_helpers

These reruns are checkpoint validations of the maintained driver surface, not new canonical baselines.

The canonical validated smoke test in this repository is:

```bash
Rscript entrypoint.R single-cohort-survival \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features GENE_D \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir results/smoke_test_single_cohort_canonical
```

That run has already been validated in this repository and writes:

- run_manifest.csv
- RUN_METADATA.md
- logs/sessionInfo.txt
- endpoint outputs such as OS/OS_results.rds, OS/OS_summary.md, and OS/results_univariate.csv

The current pooled reporting contract smoke test is:

```bash
Rscript entrypoint.R finalized-pooled-results-bundle \
  --overall-signature-table results/pooled_reporting_contract_smoke_input/overall_signature_table.csv \
  --pooled-term-table results/pooled_reporting_contract_smoke_input/pooled_term_table.csv \
  --subgroup-summary-table results/smoke_test_subgroup_contract/tables/subgroup_effect_summary.csv \
  --output-dir results/pooled_reporting_contract_smoke_shared_plotting
```

That validation confirms the pooled reporting workflow consumes the same shared subgroup-summary contract used elsewhere.

The test_data/ files are intended for smoke testing and interface validation, not scientific use. They are useful for checking:

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

This repository now uses renv for a project-local R environment.

On a fresh machine, restore the environment with:

```bash
Rscript -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org")'
Rscript -e 'renv::restore(prompt = FALSE)'
```

You can verify the project environment is aligned with:

```bash
Rscript -e 'renv::status()'
```

If package restore falls back to source builds, standard system build tools for your platform may still be required.

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