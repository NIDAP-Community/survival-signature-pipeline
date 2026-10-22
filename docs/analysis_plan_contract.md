# Single-Cohort Analysis-Plan Contract

`single-cohort-analysis-plan` executes a CSV file of independent, continuous-feature Cox-model specifications against one expression table and one metadata table. It is the maintained automation boundary for combined-cohort effects, filtered subgroup effects, and formal feature-by-group interaction tests.

The driver deliberately has no disease, subtype, cohort-name, or project-specific branches. A plan supplies the metadata columns and values appropriate to its own prepared cohort. It does not replace the separate two-dataset pooled workflows for independent comparable cohorts.

## Required Columns

| Column | Meaning |
| --- | --- |
| `analysis_id` | Unique, stable identifier for one model specification. |
| `analysis_type` | `overall`, `subgroup`, or `interaction`. |
| `feature` | Expression-row identifier for `feature_type=gene`, or a human-readable label for `feature_type=signature`. |
| `endpoint` | `OS` or `PFI` when using the CLI endpoint defaults; any other label is allowed when the row supplies both `event_column` and `time_column`. |

## Optional Columns

| Column | Meaning |
| --- | --- |
| `feature_type` | `gene` (default) or `signature`. |
| `signature_genes` | Comma-separated expression-row identifiers; required for a signature. The score is the per-sample arithmetic mean. |
| `covariates` | Comma-separated raw metadata column names. `none` or blank means unadjusted. Missing requested columns fail that plan row rather than silently changing its model. |
| `filter_column`, `filter_values` | Optional exact-match filter applied before fitting; both are required together. Values are comma-separated. Useful for a subgroup or a restricted interaction analysis. |
| `subgroup_type`, `subgroup_label` | Display metadata for an `overall` or `subgroup` primary-effect row. |
| `interaction_column` | Required for `interaction`; treated as a categorical factor. |
| `interaction_reference` | Optional observed factor level to use as the interaction reference. |
| `interaction_label` | Optional display label for the interaction axis. |
| `event_column`, `time_column` | Explicit endpoint mapping. They must be supplied together and take precedence over the `OS`/`PFI` CLI mappings. |

`overall` fits `Surv(time, event) ~ feature + covariates`. `subgroup` first applies the optional filter and then fits the same model. `interaction` first applies the optional filter and fits `Surv(time, event) ~ feature * interaction_factor + covariates`.

## Outputs

The driver writes run-scoped outputs under `results/`:

- `tables/analysis_plan_resolved.csv`: resolved input plan.
- `tables/analysis_plan_status.csv`: one completed or failed status row per plan row.
- `tables/model_coefficients.csv`: every fitted Cox coefficient.
- `tables/subgroup_effect_summary.csv`: primary feature-effect rows and interaction rows in the shared subgroup-summary contract.

For `overall` and `subgroup` rows, `subgroup_effect_summary.csv` contains the continuous feature coefficient. For `interaction` rows, it contains only feature-by-factor interaction coefficients and marks them `subgroup_type=interaction_model`; subgroup effects and interaction tests therefore remain distinguishable. `p_value_bh` is calculated separately within each feature, endpoint, and analysis type among emitted primary rows.

Invalid plan structure fails before execution. A model-specific issue—such as an unavailable feature, missing metadata column, empty filter, or non-estimable interaction—does not discard other independent rows: it is recorded in `analysis_plan_status.csv`.

## Generic Example

```bash
./scripts/run_rscript_project.sh entrypoint.R single-cohort-analysis-plan \
  --expr-file <expression.tsv> \
  --meta-file <metadata.csv> \
  --plan-file <analysis_plan.csv> \
  --output-dir results/<run_name>
```

The synthetic [`test_data/single_cohort_analysis_plan.csv`](../test_data/single_cohort_analysis_plan.csv) is an interface fixture, not a recommended scientific analysis design.
