# Smoke-Test Data

These files are small synthetic smoke-test inputs for interface validation and workflow checks.

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
Rscript entrypoint.R single-cohort-survival \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata_no_subtype.csv \
  --features AMD1 \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir results/smoke_test_single_cohort_canonical
```

This command exercises a real maintained workflow through entrypoint.R while also validating that subtype is optional when subtype-dependent options are not requested.

## Files

- tcga_expression.tsv
- tcga_metadata.csv
- tcga_metadata_no_subtype.csv
- metabric_expression.tsv
- metabric_metadata.csv
- single_cohort_subgroup_summary.csv

## Example Commands

Two-dataset interface smoke test:

```bash
Rscript entrypoint.R piecewise-dataset-interaction \
  --tcga-expr-file test_data/tcga_expression.tsv \
  --tcga-meta-file test_data/tcga_metadata.csv \
  --metabric-expr-file test_data/metabric_expression.tsv \
  --metabric-meta-file test_data/metabric_metadata.csv \
  --features Signature \
  --endpoints OS
```

Subtype-specific smoke test:

```bash
Rscript entrypoint.R piecewise-subtype-specific \
  --tcga-expr-file test_data/tcga_expression.tsv \
  --tcga-meta-file test_data/tcga_metadata.csv \
  --metabric-expr-file test_data/metabric_expression.tsv \
  --metabric-meta-file test_data/metabric_metadata.csv \
  --subtype-list BRCA_Basal \
  --adjust-covariates none \
  --features Signature \
  --endpoints OS
```

Alternative single-cohort smoke test:

```bash
Rscript entrypoint.R single-cohort-survival \
  --expr-file test_data/tcga_expression.tsv \
  --meta-file test_data/tcga_metadata.csv \
  --features Signature \
  --endpoints OS \
  --adjust-covariates none
```

Single-cohort reporting heatmap smoke test:

```bash
Rscript entrypoint.R single-cohort-reporting-bundle \
  --input-run-dir results/smoke_test_single_cohort \
  --subgroup-summary-csv test_data/single_cohort_subgroup_summary.csv
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

## Subtype Guardrail Validation Matrix

The single-cohort driver treats `subtype_column` as optional unless subtype filtering or subtype-adjusted modeling is requested.

| Case | Metadata file | Subtype filter | Adjust covariates | Expected behavior |
| --- | --- | --- | --- | --- |
| no subtype column, no subtype filter, no subtype adjustment | tcga_metadata_no_subtype.csv | ALL | none | should run |
| no subtype column, subtype filter requested | tcga_metadata_no_subtype.csv | BRCA_Basal | none | should fail clearly |
| no subtype column, subtype adjustment requested | tcga_metadata_no_subtype.csv | ALL | subtype | should fail clearly |
| subtype column present, no subtype-specific options | tcga_metadata.csv | ALL | none | should run |
| subtype column present, subtype-specific options enabled | tcga_metadata.csv | BRCA_Basal,BRCA_LumA | subtype | should run |