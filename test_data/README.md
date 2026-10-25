# Smoke-Test Data

These files are small synthetic smoke-test inputs for interface validation and workflow checks.

The documented primary interface for these checks is the root entrypoint.

They are intended to help verify that:

- the root entrypoint dispatches correctly
- maintained drivers can read expected columns
- run-scoped outputs are created under results/
- provenance metadata is written

They are not intended for scientific interpretation or for validating the statistical behavior of the models.

## Canonical Maintained Smoke Test

The canonical validated smoke test for this repository currently targets the maintained single-cohort workflow.

Input pair:

- tcga_expression.tsv
- tcga_metadata_no_subtype.csv

Canonical command:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-survival \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features GENE_D \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir results/smoke_test_single_cohort_canonical
```

This command exercises a real maintained workflow through entrypoint.R while also validating that subtype is optional when subtype-dependent options are not requested.

## Files

- tcga_expression.tsv
- binary_biomarker_matrix.tsv
- tcga_metadata.csv
- tcga_metadata_no_subtype.csv
- metabric_expression.tsv
- metabric_metadata.csv
- single_cohort_subgroup_summary.csv
- gbm_lgg_p53beta_adjusted_effects_subgroup_summary.csv
- gbm_lgg_p53beta_interactions_subgroup_summary.csv
- single_cohort_analysis_plan.csv

## Binary Biomarker Smoke Test

The binary-biomarker fixture is a synthetic, feature-by-subject 0/1 table. It
is an interface test for patient-level binary calls such as mutation status;
it is not a scientifically interpretable mutation dataset.

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-binary-biomarker-survival \
  --biomarker-file test_data/binary_biomarker_matrix.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features TP53_mutated \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir results/binary_biomarker_smoke
```

## Example Commands

Two-dataset interface smoke test:

```bash
./scripts/run_rscript_project.sh entrypoint.R piecewise-dataset-interaction \
  --dataset1-expr-file test_data/tcga_expression.tsv \
  --dataset1-meta-file test_data/tcga_metadata.csv \
  --dataset2-expr-file test_data/metabric_expression.tsv \
  --dataset2-meta-file test_data/metabric_metadata.csv \
  --features Signature \
  --endpoints OS
```

Subtype-specific smoke test:

```bash
./scripts/run_rscript_project.sh entrypoint.R piecewise-subtype-specific \
  --dataset1-expr-file test_data/tcga_expression.tsv \
  --dataset1-meta-file test_data/tcga_metadata.csv \
  --dataset2-expr-file test_data/metabric_expression.tsv \
  --dataset2-meta-file test_data/metabric_metadata.csv \
  --subtype-list BRCA_Basal \
  --adjust-covariates none \
  --features Signature \
  --endpoints OS
```

Alternative single-cohort smoke test:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-survival \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata.csv \
  --features Signature \
  --endpoints OS \
  --adjust-covariates none
```

Detectable-versus-undetectable smoke test:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-survival \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features GENE_D \
  --endpoints OS \
  --split-method detectable \
  --detection-threshold 5 \
  --include-groups low,high \
  --adjust-covariates none \
  --output-dir results/validation_detectable_single_cohort
```

This command is a maintained interface check for the detectability mode. It validates that the run writes `Undetectable` and `Detectable` group outputs through the maintained single-cohort path.

The maintained driver also accepts `--include-groups Undetectable,Detectable` as an alias for `--include-groups low,high` in detectability mode.

Single-cohort reporting heatmap smoke test:

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-reporting-bundle \
  --input-run-dir results/smoke_test_single_cohort \
  --subgroup-summary-table test_data/single_cohort_subgroup_summary.csv
```

The subgroup-summary fixture uses the shared downstream heatmap contract:

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

GBM/LGG p53beta adjusted-effects forest fixture:

```bash
./scripts/run_rscript_project.sh entrypoint.R subgroup-forest-bundle \
  --subgroup-summary-table test_data/gbm_lgg_p53beta_adjusted_effects_subgroup_summary.csv \
  --subgroup-type-filter survival_model \
  --output-dir results/validation_gbm_lgg_adjusted_forest
```

GBM/LGG p53beta interaction forest fixture:

```bash
./scripts/run_rscript_project.sh entrypoint.R subgroup-forest-bundle \
  --subgroup-summary-table test_data/gbm_lgg_p53beta_interactions_subgroup_summary.csv \
  --subgroup-type-filter interaction_model \
  --output-dir results/validation_gbm_lgg_interaction_forest
```

These GBM/LGG fixtures are real field-report tables from an external p53beta signature project. They are intended for reporting-interface and generalization tests, not for standalone scientific interpretation inside this repository.

## Declarative Analysis-Plan Smoke Test

`single_cohort_analysis_plan.csv` is a synthetic interface fixture for the maintained `single-cohort-analysis-plan` workflow. It covers a continuous gene effect, a mean-expression signature effect, a filtered-group effect, and a feature-by-group interaction without adding project-specific logic to the driver.

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-analysis-plan \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata.csv \
  --plan-file test_data/single_cohort_analysis_plan.csv \
  --output-dir results/analysis_plan_smoke
```

See `docs/analysis_plan_contract.md` for the general plan schema and output contract.

## Subtype Guardrail Validation Matrix

The single-cohort driver treats `subtype_column` as optional unless subtype filtering or subtype-adjusted modeling is requested.

| Case | Metadata file | Subtype filter | Adjust covariates | Expected behavior |
| --- | --- | --- | --- | --- |
| no subtype column, no subtype filter, no subtype adjustment | tcga_metadata_no_subtype.csv | ALL | none | should run |
| no subtype column, subtype filter requested | tcga_metadata_no_subtype.csv | BRCA_Basal | none | should fail clearly |
| no subtype column, subtype adjustment requested | tcga_metadata_no_subtype.csv | ALL | subtype | should fail clearly |
| subtype column present, no subtype-specific options | tcga_metadata.csv | ALL | none | should run |
| subtype column present, subtype-specific options enabled | tcga_metadata.csv | BRCA_Basal,BRCA_LumA | subtype | should run |
