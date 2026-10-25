# Survival Signature Pipeline

An R command-line workflow for survival analysis of prepared molecular features and clinical outcomes. It supports gene expression and expression signatures, binary molecular biomarkers such as mutation status, a declarative plan of adjusted/subgroup/interaction Cox models, and pooled analyses of two harmonized independent cohorts. All maintained workflows write run-scoped tables, figures, and provenance metadata.

## Features

- Continuous and grouped survival analysis for individual genes and z-scored expression signatures
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

The same maintained CLI is available as a portable OCI image. Build the image locally, then run the containerized smoke test:

```bash
container/build_docker.sh survival-signature-pipeline:local
container/smoke_test.sh survival-signature-pipeline:local
```

The image is built for `linux/amd64`, matching the multi-gene-correlation container standard. It reads caller-provided files from `/data` and writes run artifacts to `/output`; see [the container contract](docs/container_contract.md) for a Docker invocation and integration boundary.

### Published image

The validated image is published to GitHub Container Registry (GHCR). For a
reproducible run, pull and launch it by immutable digest rather than by a
mutable tag:

```bash
docker pull ghcr.io/nidap-community/survival-signature-pipeline@sha256:a25e77f298eb602bf0a5e088a295cce783dd987baf23b7dada4d1cafea9ed35f
```

The corresponding source-revision tag is
`99cbdf1cfdcdb4ecf9c0ed5cbe4ac6ab996f7092`; it is provided for discovery,
while the digest above is the deployment identity. A Docker invocation should
bind caller-controlled inputs read-only at `/data` and a fresh output directory
at `/output`. The repository helper enforces that layout:

```bash
container/run_docker.sh \
  ghcr.io/nidap-community/survival-signature-pipeline@sha256:a25e77f298eb602bf0a5e088a295cce783dd987baf23b7dada4d1cafea9ed35f \
  /path/to/input_directory \
  /path/to/output_directory \
  single-cohort-survival \
  --expr-file /data/expression.tsv \
  --meta-file /data/metadata.csv \
  --features Signature \
  --signature-genes GENE1,GENE2 \
  --output-dir /output
```

HPC users may convert the same OCI image with Apptainer or Singularity. Test
the converted image in the target environment before using it for a production
job:

```bash
apptainer pull survival-signature-pipeline.sif \
  docker://ghcr.io/nidap-community/survival-signature-pipeline@sha256:a25e77f298eb602bf0a5e088a295cce783dd987baf23b7dada4d1cafea9ed35f
```

See [container/README.md](container/README.md) for the distribution model and
[the container contract](docs/container_contract.md) for the full runtime
interface. The repository's container-validation workflow runs an amd64 image
smoke test on pull requests and pushes to `main`; it validates but does not
publish images. A calling application should follow the [AWS adapter
contract](docs/aws_adapter_contract.md) for validated request handling, data
staging, and artifact collection.

## Example Output

### How to read this demonstration

These assets are one **exploratory demonstration**, not a clinical classifier
or a claim about an individual patient. They were generated from the public
METABRIC breast-cancer cohort with an unweighted 11-gene PAM50 proliferation
module. Each signature member was z-scored within the eligible cohort before
the module score was calculated. The same run produced the two KM plots, the
two forest plots, and the table below.

For the visual display, the eligible samples were divided into within-cohort
low, mid, and high score tertiles (`n = 657` per group; `N = 1,971` before
covariate complete-case filtering). Tertiles make the curves easier to inspect;
they are not the only analysis. The statistics panel also reports continuous
Cox models that retain the full score information rather than discarding it
into groups.

| Label in the KM statistics panel | Analysis shown | Reference / interpretation |
| --- | --- | --- |
| `cont (unadj)` | Cox proportional-hazards model with the signature score as a continuous predictor | HR for a one-unit increase in the calculated signature score; all eligible samples; no clinical covariates. |
| `cont (adj)` | Continuous Cox model adjusted for the listed clinical covariates | Same continuous signature effect after complete-case adjustment. |
| `unadj H vs L`, `unadj M vs L` | Cox model with the three-level tertile group as the predictor | Low tertile is the reference; HRs compare the high or mid tertile with low. |
| `adj H vs L`, `adj M vs L` | Three-level tertile Cox model adjusted for the listed clinical covariates | Same tertile comparisons after complete-case adjustment. |

The plotted curves are **unadjusted Kaplan-Meier estimates**. Short marks on a
curve are censoring times; the vertical dashed lines mark the estimated median
survival time for each tertile. The `p` values printed in the statistics panel
and in the table below are two-sided Wald tests for Cox-model coefficients.
They are not an omnibus log-rank p value. A global log-rank comparison is a
separate, unadjusted test and should not be substituted for an adjusted Cox
comparison.

The demonstration uses the input cohort's declared endpoint fields: `OS` /
`OS_time` for overall survival and `PFI` / `PFI_time` for progression-free
interval. Event definitions and follow-up origin come from the supplied
metadata; users must confirm those cohort-specific definitions before comparing
endpoints or cohorts.

![Overall-survival Kaplan-Meier plot for the PAM50 proliferation-module demonstration](assets/pam50_proliferation_metabric_os_km.png)

**Overall survival (OS).** The three curves are the unadjusted KM estimates
for the score tertiles. The panel beneath them reports both continuous-score
and tertile Cox effects. The adjusted rows use `Diagnosis_Age`, `tumor_size`,
`N`, `grade`, `ER_status_norm`, and `treatment`; their smaller `N` reflects
complete-case filtering for those covariates.

![Progression-free-interval Kaplan-Meier plot for the PAM50 proliferation-module demonstration](assets/pam50_proliferation_metabric_pfi_km.png)

**Progression-free interval (PFI).** This uses the same score construction,
tertile cutpoints, and covariate set as the OS example, but fits a separate
model to the PFI event/time fields. It is therefore a separate endpoint
analysis, not a second test of the OS model.

### Overall Effect Table

This compact extract is the **grouped Cox high-versus-low tertile coefficient**
from the same three-level models shown in the KM panels. It does not reproduce
the continuous-score rows, the mid-versus-low rows, or every covariate term;
those are available in the full run outputs. `N` and events differ between
unadjusted and adjusted models because adjustment uses complete cases for all
listed covariates. P values are unadjusted Cox Wald-test p values; no
multiple-testing correction is applied to this four-row illustration.

| Endpoint | Cox model | Adjustment set | Comparison | N (events) | HR | 95% CI | Wald p value |
| --- | --- | --- | --- | ---: | ---: | --- | ---: |
| OS | Grouped, unadjusted | None | High vs Low tertile | 1,971 (891) | 1.623 | [1.374, 1.916] | 1.21e-08 |
| OS | Grouped, adjusted | Age, tumor size, N, grade, ER status, treatment | High vs Low tertile | 1,843 (828) | 1.259 | [1.030, 1.540] | 0.025 |
| PFI | Grouped, unadjusted | None | High vs Low tertile | 1,971 (506) | 3.097 | [2.412, 3.975] | 7.35e-19 |
| PFI | Grouped, adjusted | Age, tumor size, N, grade, ER status, treatment | High vs Low tertile | 1,843 (485) | 2.058 | [1.531, 2.767] | 1.74e-06 |

![Overall forest plot for the PAM50 proliferation-module demonstration](assets/pam50_proliferation_metabric_overall_forest.png)

**Overall-effect forest plot.** This is the visual counterpart to the table:
each point is a grouped high-versus-low Cox HR, each horizontal line is its
95% confidence interval, and the dashed vertical line is HR = 1 (no estimated
association). Open gray markers are unadjusted models; filled blue markers are
adjusted models. The x-axis is logarithmic, and the right-side text gives HR,
95% CI, and the same raw Wald p value. OS and PFI are displayed in separate
facets; the figure is not a pooled or meta-analytic estimate.

![Full adjusted Cox-model forest for the PAM50 proliferation-module demonstration](assets/pam50_proliferation_metabric_full_adjusted_forest.png)

**Full adjusted-model forest.** This expands the adjusted tertile model rather
than adding another test: it shows the signature's mid-versus-low and
high-versus-low coefficients plus every retained covariate coefficient for OS
and PFI. A continuous covariate HR is per one unit in the supplied metadata;
a categorical covariate HR is a displayed level versus its input reference
level. These figures demonstrate the reporting layout only; they are not an
official PAM50 classifier or a clinical prediction model. Review
proportional-hazards and time-dependent diagnostic outputs before interpreting
a grouped KM contrast or forest plot.

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

`single-cohort-survival` analyzes both `OS` and `PFI` by default (`--endpoints OS,PFI`). Pass `--endpoints OS` or `--endpoints PFI` only when intentionally restricting a run to one endpoint.

Prepare filtering and sample inclusion upstream. Maintained drivers analyze the rows supplied to them; they do not impose a tumor-type or sample-type filter.

### Binary biomarker matrix

Use `single-cohort-binary-biomarker-survival` for a precomputed, feature-by-subject table where `0` is the reference state and `1` is the comparison state. It is appropriate for mutation calls after the calling workflow has made and documented the call; it does not infer mutations from expression values. Missing values are permitted and mean that the feature was not evaluable for that subject.

```text
Feature       Sample_1  Sample_2  Sample_3
TP53_mutated  1         0         1
```

Feature identifiers are arbitrary labels in this workflow, so it can represent a gene-level mutation, a copy-number event after binary calling, or another already-defined binary molecular state. The matrix and metadata must use the same subject identifiers.

### Continuous biomarker matrix

Use `single-cohort-continuous-biomarker-survival` for a precomputed, feature-by-subject numeric table, initially for copy-number values such as log2 ratios. Missing values are permitted; every non-missing value must be finite and each requested feature must vary in the supplied cohort. The workflow reports continuous Cox effects per within-cohort standard deviation and uses the selected tertile, quartile, or median split only for the KM display.

```text
Feature                  Sample_1  Sample_2  Sample_3
MYC_log2_copy_number     0.42      -0.11     0.67
```

It does not normalize copy number or turn a numeric value into a deletion/amplification call. Those choices belong to the upstream data-preparation contract.

### Immune-score and cell-fraction contexts

The continuous workflow also accepts prepared immune measurements. Use `--biomarker-context estimate_scores` for ESTIMATE-derived score rows such as `ESTIMATEScore`, `ImmuneScore`, and `StromalScore`. The container does not calculate ESTIMATE; the staged table must retain the upstream method, version, and cohort provenance.

Use `--biomarker-context cell_fractions` only for a complete cell-fraction composition. It requires every subject to have every staged component, every fraction to be in `[0, 1]`, and each subject column to sum to one within `--fraction-sum-tolerance` (default `0.02`). A CD8 T-cell association in this mode is an exploratory one-feature analysis, not a multivariable compositional survival model. Do not use this context for unbounded immune scores or CIBERSORT absolute-mode values.

### Gene identifiers

The maintained gene-consuming workflows can normalize user-supplied human gene symbols with the pinned `l2psupp::updategenes()` resolver (`--gene-resolution l2p`). This opt-in behavior leaves existing calls unchanged. The expression matrix is never renamed or rewritten. An exact input row is preferred; otherwise, an updated symbol is used only when it matches an existing expression row (in both matrices for a two-cohort workflow). Every run writes `tables/gene_identifier_resolution.csv`, preserving the submitted name, L2P result, package version, matched expression row, and any unresolved gene.

`--gene-resolution none` is the default for backward compatibility; use it when requested identifiers are already the intended expression-row identifiers. The run fails rather than dropping an unresolved requested gene or double-counting two inputs that resolve to the same row.

## Choose a Workflow

| Workflow | Use when |
| --- | --- |
| `single-cohort-survival` | You need grouped and continuous gene/signature survival analysis, Kaplan-Meier plots, and per-endpoint outputs for one cohort. |
| `single-cohort-binary-biomarker-survival` | You have a prepared 0/1 molecular-feature matrix, such as mutation status, and need mutated/positive versus wild-type/negative KM and Cox results. |
| `single-cohort-continuous-biomarker-survival` | You have a prepared numeric molecular-feature matrix, such as copy-number log2 ratios, and need standardized continuous Cox effects plus grouped KM displays. |
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
  --gene-resolution l2p \
  --features <gene_or_Signature> \
  --signature-genes <gene1,gene2,gene3> \
  --endpoints OS,PFI \
  --adjust-covariates <metadata_column1,metadata_column2> \
  --output-dir results/<run_name>
```

Use `--features Signature` only when `--signature-genes` supplies the feature members. To restrict an analysis to a metadata-defined group, add `--group-column <column>` and `--group-filter <value1,value2>`.

Signature genes are z-scored within the eligible cohort. The default score is their sum. The default Kaplan-Meier display uses low/mid/high tertiles, while continuous signature Cox models use **all eligible samples**; this preserves the full score information for the primary estimate even when a grouped plot is shown. Use `--signature-aggregation mean` when the analysis specification calls for an average score, or `--signature-continuous-cohort grouped` only when a deliberately group-restricted continuous analysis is required.

Use `--signature-label <display_label>` to give a signature plot a clear human-readable title without changing the score calculation.

### Complete analysis is the default

A normal run retains result objects, checks proportional-hazards assumptions,
and performs signature time-dependent analysis when indicated. Together with
the default `OS,PFI` endpoint set, these artifacts let users inspect the full
analysis rather than seeing only a Kaplan-Meier figure. Use the reporting
bundle after a completed run to add the overall forest plot.

For a rapid exploratory figure only, users may explicitly opt out with
`--run-ph-checks false --run-time-dependent false --save-results-rds false`.
That mode is not a substitute for the complete diagnostic run.

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

## Run a Binary Molecular-Biomarker Analysis

This workflow keeps molecular calling upstream and models the submitted binary values directly. In this example, `1` is labeled `Mutated` and is compared with `0`, labeled `Wild type`:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-binary-biomarker-survival \
  --biomarker-file test_data/binary_biomarker_matrix.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features TP53_mutated \
  --negative-label "Wild type" \
  --positive-label Mutated \
  --endpoints OS,PFI \
  --adjust-covariates none \
  --output-dir results/binary_biomarker_smoke
```

For each endpoint, the run writes the KM plot, `grouped_biomarker_cox_results.csv`, an endpoint summary, result object, and the same run-level provenance files used by the expression workflow. The KM plot and log-rank test are unadjusted; grouped Cox rows provide the `1` versus `0` hazard ratio and, when requested, its covariate-adjusted estimate. The binary workflow accepts only `0`, `1`, or missing values and does not run L2P because its feature labels need not be gene symbols.

## Run a Continuous Molecular-Biomarker Analysis

This workflow accepts prepared numeric data without treating values as expression. For example, a log2 copy-number feature can be analyzed as follows:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-continuous-biomarker-survival \
  --biomarker-file test_data/continuous_biomarker_matrix.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features MYC_log2_copy_number \
  --split-method tertiles \
  --endpoints OS,PFI \
  --adjust-covariates none \
  --output-dir results/continuous_biomarker_smoke
```

The primary continuous Cox result is per one standard deviation of the submitted feature among eligible cohort members. The selected split affects only grouped KM and grouped Cox results. Use a separately documented binary matrix when a study-specific amplification/deletion threshold is scientifically required.

### ESTIMATE immune scores

Run prepared ESTIMATE scores through the same continuous workflow, identifying the context in provenance and summaries:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-continuous-biomarker-survival \
  --biomarker-file test_data/estimate_score_matrix.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features ImmuneScore \
  --biomarker-context estimate_scores \
  --endpoints OS,PFI \
  --adjust-covariates none \
  --output-dir results/estimate_score_smoke
```

### CIBERSORT-like cell fractions

For complete relative-abundance output, use the fraction context. Include all components from the composition in the staged table even when analyzing one feature:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-continuous-biomarker-survival \
  --biomarker-file test_data/cell_fraction_matrix.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features CD8_T_cells \
  --biomarker-context cell_fractions \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir results/cell_fraction_smoke
```

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
│   ├── gene_identifier_resolution.csv          # submitted-to-expression gene mapping
│   ├── model_coefficients.csv                 # analysis-plan workflow
│   ├── subgroup_effect_summary.csv            # analysis-plan/subgroup workflows
│   └── ...
├── figures/                                   # reporting and forest workflows
└── OS/ or PFI/                                # single-cohort endpoint outputs
    ├── <endpoint>_summary.md                  # GitHub-rendered model tables, including covariates
    ├── adjusted_model_full_coefficients.csv
    ├── grouped_gene_cox_results.csv
    ├── results_univariate.csv
    └── <endpoint>_results.rds
```

Output names vary by workflow. Run `--help` before a new workflow and retain the corresponding `RUN_METADATA.md` with any result shared outside the repository.

### Overall forest plot for a single-cohort run

The single-cohort driver writes its endpoint analyses first. With its default `--save-results-rds true`, create an overall forest plot from that completed run with:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-reporting-bundle \
  --input-run-dir results/<single_cohort_run> \
  --signature-label <display_label> \
  --output-dir results/<single_cohort_reporting_bundle>
```

This writes the primary-effect forest and its CSV/GitHub-rendered Markdown table under `tables/`. When the upstream run has a grouped adjusted signature model, it also writes the full adjusted-model forest and table:

- `figure1_single_cohort_overall_forest`: primary individual-gene and signature effects for each completed endpoint.
- `figure2_single_cohort_full_adjusted_model_forest`: the signature comparison and every retained covariate coefficient from each endpoint's grouped adjusted Cox model. Factor levels are labeled against their model reference level.

Use the normal analysis defaults when these reports are needed, rather than a figure-only run with `--save-results-rds false`.

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

## Distribution and versioning

Use the published OCI digest for a reproducible deployment and retain it with
the run's `RUN_METADATA.md`. A formal semantic-version Git tag and GitHub
Release are intentionally deferred; when one is created, its notes and
versioned image reference will be recorded under [docs/releases/](docs/releases/).

## License

This project is available under the [MIT License](LICENSE).

## Support

Open an issue at [NIDAP-Community/survival-signature-pipeline](https://github.com/NIDAP-Community/survival-signature-pipeline/issues) for bugs, questions, or feature requests. Consult [docs/analysis_plan_contract.md](docs/analysis_plan_contract.md) for the detailed analysis-plan contract and [docs/architecture.md](docs/architecture.md) for development-oriented architecture notes.
