# Survival Signature Pipeline

An R command-line workflow for survival analysis of gene expression and expression signatures. It supports one prepared cohort, a declarative plan of adjusted/subgroup/interaction Cox models, and pooled analyses of two harmonized independent cohorts. All maintained workflows write run-scoped tables, figures, and provenance metadata.

## Features

- Continuous and grouped survival analysis for individual genes and mean-expression signatures
- Kaplan-Meier plots and Cox models for a single cohort
- Declarative analysis plans for adjusted overall effects, filtered subgroup effects, and feature-by-group interaction tests
- Two-cohort piecewise Cox models for comparable, independently collected cohorts
- Generic forest plots from a shared subgroup-effect table
- Run metadata, input checksums, session information, and run-scoped output directories

## Quick Start

### Local R setup

The project lockfile targets R 4.5.1. Install `renv` once, then restore the project library from the repository root:

```bash
Rscript --vanilla -e 'if (!requireNamespace("renv", quietly = TRUE)) install.packages("renv", repos = "https://cloud.r-project.org"); renv::restore(project = getwd(), prompt = FALSE)'
```

Use the project runner for subsequent command-line runs. It selects the restored library and avoids a known noninteractive `renv` startup problem on some local systems:

```bash
./scripts/run_rscript_project.sh entrypoint.R --help
```

List every maintained workflow and its short description:

```bash
./scripts/run_rscript_project.sh entrypoint.R list
```

### Container setup

The same maintained CLI is available as a portable OCI image. Build the image, then run the containerized smoke test:

```bash
container/build_docker.sh survival-signature-pipeline:local
container/smoke_test.sh survival-signature-pipeline:local
```

Use `TARGET_PLATFORM=linux/arm64` or `TARGET_PLATFORM=linux/amd64` when building for a specific deployment target. The image reads caller-provided files from `/data` and writes run artifacts to `/output`; see [the container contract](docs/container_contract.md) for a Docker invocation and integration boundary.

The repository's container-validation workflow runs a native image smoke test and builds the ARM64 target on pull requests and pushes to `main`; it does not publish an image.

## Input Files

### Expression matrix

Supply a tab-separated expression matrix with genes in rows and samples in columns. The gene-identifier column defaults to `Gene` and can be changed with `--gene-column`.

```text
Gene    Sample_1    Sample_2    Sample_3
GENE_A  1.2         0.8         1.5
GENE_B  2.1         1.9         2.4
```

### Metadata table

Supply a CSV metadata table with one row per sample. The sample identifier defaults to `sampnames` and can be changed with `--sample-column`. For standard endpoints, provide paired event/time columns:

```text
sampnames,OS,OS_time,PFI,PFI_time,Age,Group
Sample_1,1,18,1,12,55,Control
Sample_2,0,42,0,42,61,Treatment
```

`OS`/`OS_time` and `PFI`/`PFI_time` are the default endpoint mappings. Change them with the corresponding CLI flags. The analysis-plan workflow also allows a row to define any explicit `event_column` and `time_column` pair.

Prepare filtering and sample inclusion upstream. Maintained drivers analyze the rows supplied to them; they do not impose a tumor-type or sample-type filter.

## Choose a Workflow

| Workflow | Use when |
| --- | --- |
| `single-cohort-survival` | You need grouped and continuous gene/signature survival analysis, Kaplan-Meier plots, and per-endpoint outputs for one cohort. |
| `single-cohort-analysis-plan` | You need a reproducible collection of adjusted overall effects, filtered subgroup effects, or formal feature-by-group interaction tests for one prepared cohort. |
| `subgroup-forest-bundle` | You already have a shared subgroup-effect summary table and need a forest plot. |
| `piecewise-dataset-interaction` | You have two harmonized, independent comparable cohorts and need pooled piecewise models with dataset interactions. |
| `piecewise-subtype-specific` | You have two harmonized cohorts and need pooled piecewise models within defined subgroup values. |
| `single-cohort-reporting-bundle` | You need reporting figures and tables from a completed `single-cohort-survival` run. |
| `finalized-pooled-results-bundle` | You need reporting figures and summary text from completed pooled result tables. |

Do not use the two-cohort workflows to represent biologically different disease states as though they were independent validation cohorts. Use a single-cohort analysis plan with explicit covariates, filters, and interaction axes instead.

## Run a Single-Cohort Analysis

This example runs the bundled smoke-test cohort through the maintained single-cohort workflow:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-survival \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features GENE_D \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir results/single_cohort_smoke
```

For a real run, make the feature, endpoint, covariate, and output choices explicit:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-survival \
  --expr-file <expression.tsv> \
  --meta-file <metadata.csv> \
  --sample-column <sample_id_column> \
  --gene-column <gene_id_column> \
  --features <gene_or_Signature> \
  --signature-genes <gene1,gene2,gene3> \
  --endpoints OS,PFI \
  --adjust-covariates <metadata_column1,metadata_column2> \
  --output-dir results/<run_name>
```

Use `--features Signature` only when `--signature-genes` supplies the feature members. To restrict an analysis to a metadata-defined group, add `--group-column <column>` and `--group-filter <value1,value2>`.

### Detectable versus undetectable genes

For an individual gene, use an explicit threshold:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-survival \
  --expr-file <expression.tsv> \
  --meta-file <metadata.csv> \
  --features <gene> \
  --endpoints OS \
  --split-method detectable \
  --detection-threshold 0 \
  --include-groups Undetectable,Detectable \
  --adjust-covariates none \
  --output-dir results/<detectability_run>
```

This grouping is only available for individual genes, not `Signature`.

## Run a Declarative Analysis Plan

Use `single-cohort-analysis-plan` for a planned collection of continuous Cox models. Each CSV row describes an independent `overall`, `subgroup`, or `interaction` model; cohort-specific column names, filters, covariates, signature members, and interaction references live in the plan rather than in the code.

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-analysis-plan \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata.csv \
  --plan-file test_data/single_cohort_analysis_plan.csv \
  --output-dir results/analysis_plan_smoke
```

The required plan columns are:

| Column | Description |
| --- | --- |
| `analysis_id` | Unique identifier for the model row. |
| `analysis_type` | `overall`, `subgroup`, or `interaction`. |
| `feature` | Gene identifier, or display label for a signature. |
| `endpoint` | `OS` or `PFI`, or a custom label when explicit endpoint columns are included. |

See [the analysis-plan contract](docs/analysis_plan_contract.md) for the complete optional-column schema and output definitions.

## Run a Two-Cohort Analysis

Use the pooled workflows only with two independently collected, harmonized cohorts representing a comparable disease setting. Label the cohorts explicitly and supply all four input tables:

```bash
./scripts/run_rscript_project.sh entrypoint.R piecewise-dataset-interaction \
  --dataset1-expr-file <cohort_1_expression.tsv> \
  --dataset1-meta-file <cohort_1_metadata.csv> \
  --dataset1-label <cohort_1_name> \
  --dataset2-expr-file <cohort_2_expression.tsv> \
  --dataset2-meta-file <cohort_2_metadata.csv> \
  --dataset2-label <cohort_2_name> \
  --features Signature \
  --signature-genes <gene1,gene2,gene3> \
  --endpoints OS,PFI \
  --output-dir results/<pooled_run>
```

For subtype- or other stratum-specific pooled models, use `piecewise-subtype-specific --help` to review its required dataset-specific column flags and `--subtype-list` interface.

## Make a Forest Plot from Subgroup Results

All subgroup-producing workflows can emit or consume the shared summary columns below:

```text
subgroup, subgroup_type, feature, endpoint, comparison,
HR, CI_low, CI_high, p_value, p_value_bh
```

Create a generic forest plot from such a table:

```bash
./scripts/run_rscript_project.sh entrypoint.R subgroup-forest-bundle \
  --subgroup-summary-table <subgroup_effect_summary.csv> \
  --subgroup-type-filter <subgroup_type> \
  --output-dir results/<forest_run>
```

## Output Files

Each run writes to the specified `results/<run_name>/` directory and includes `RUN_METADATA.md` plus `logs/sessionInfo.txt`.

Typical analysis outputs include:

```text
results/<run_name>/
├── RUN_METADATA.md
├── logs/
│   └── sessionInfo.txt
├── tables/
│   ├── model_coefficients.csv                 # analysis-plan workflow
│   ├── subgroup_effect_summary.csv            # analysis-plan/subgroup workflows
│   └── ...
├── figures/                                   # reporting and forest workflows
└── OS/ or PFI/                                # single-cohort endpoint outputs
    ├── adjusted_model_full_coefficients.csv
    ├── grouped_gene_cox_results.csv
    ├── results_univariate.csv
    └── <endpoint>_results.rds
```

Output names vary by workflow. Run `--help` before a new workflow and retain the corresponding `RUN_METADATA.md` with any result shared outside the repository.

## Help and Smoke Tests

Every maintained workflow provides workflow-specific help:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-survival --help
./scripts/run_rscript_project.sh entrypoint.R single-cohort-analysis-plan --help
./scripts/run_rscript_project.sh entrypoint.R piecewise-dataset-interaction --help
```

The small files in `test_data/` are interface smoke-test fixtures only; they are not suitable for scientific interpretation. See [test_data/README.md](test_data/README.md) for the supported smoke-test commands.

## Requirements

- R 4.5.1 or a compatible R release
- The project library restored from `renv.lock`
- Tab-separated expression matrices and CSV metadata tables as described above

## Support

Open an issue at [NIDAP-Community/survival-signature-pipeline](https://github.com/NIDAP-Community/survival-signature-pipeline/issues) for bugs, questions, or feature requests. Consult [docs/analysis_plan_contract.md](docs/analysis_plan_contract.md) for the detailed analysis-plan contract and [docs/architecture.md](docs/architecture.md) for development-oriented architecture notes.
