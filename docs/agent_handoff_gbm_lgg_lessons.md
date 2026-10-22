# Agent Handoff: GBM/LGG p53beta Survival Lessons

Date: 2026-05-21

## Ask For The Pipeline-Optimization Agent

Use the GBM/LGG p53beta field report and fixture tables to generalize the survival-signature-pipeline for single-cohort stratified survival analyses. The goal is not to encode a GBM/LGG special case, but to make the pipeline handle this kind of project cleanly.

## Read First

- `docs/case_studies/gbm_lgg_p53beta_survival.md`
- `docs/decision_log.md`
- `docs/session_notes.md`
- `test_data/gbm_lgg_p53beta_adjusted_effects_subgroup_summary.csv`
- `test_data/gbm_lgg_p53beta_interactions_subgroup_summary.csv`

## Key Lessons To Preserve

- Do not force biologically distinct disease-state comparisons into the two-dataset pooled validation workflow.
- Keep independent-cohort pooling separate from single-cohort disease/subgroup modeling.
- Treat IDH status as central biology in glioma analyses.
- Treat stage proxies explicitly; do not label disease group or histologic grade as true clinical stage.
- Keep subgroup survival effects separate from formal interaction/effect-modification tests.
- Use the shared subgroup-summary contract as the downstream reporting boundary.

## Concrete Pipeline Issues Observed

- Normal `Rscript` from the pipeline repo can hang during `.Rprofile` startup because `.Rprofile` calls `renv/activate.R`.
- `Rscript --vanilla` avoids the hang.
- `RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript ...` also avoids the hang.
- `entrypoint.R` argument forwarding via `system2()` can mishandle long quoted strings containing spaces.
- Direct driver invocation worked for `drivers/run_subgroup_forest_bundle.R`.

## Candidate Work Items

- Add a noninteractive-safe execution path or documented flag for disabling `renv` autoload during reporting workflows.
- Diagnose the `renv` autoload hang and document the expected pipeline startup behavior.
- Harden `entrypoint.R` argument forwarding for quoted values.
- Add a maintained aggregator that converts multiple single-cohort pipeline run outputs into the shared subgroup-summary contract.
- Add maintained support for signature-by-group interaction models in a single-cohort run plan.
- Extend `subgroup-forest-bundle` examples and tests to cover the GBM/LGG adjusted-effect and interaction fixture tables.
- Consider a config-driven run plan for combined, subgroup, and interaction models.

## Validation Targets

At minimum, the optimizing agent should be able to run:

```bash
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript entrypoint.R subgroup-forest-bundle \
  --subgroup-summary-table test_data/gbm_lgg_p53beta_adjusted_effects_subgroup_summary.csv \
  --subgroup-type-filter survival_model \
  --output-dir results/validation_gbm_lgg_adjusted_forest
```

and:

```bash
RENV_CONFIG_AUTOLOADER_ENABLED=FALSE Rscript entrypoint.R subgroup-forest-bundle \
  --subgroup-summary-table test_data/gbm_lgg_p53beta_interactions_subgroup_summary.csv \
  --subgroup-type-filter interaction_model \
  --output-dir results/validation_gbm_lgg_interaction_forest
```

Both should complete without hanging and should write figure, table, and run metadata outputs.
