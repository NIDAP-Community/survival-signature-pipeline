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
- Use `scripts/run_rscript_project.sh` after `renv::restore()` for a noninteractive-safe execution path: it disables the autoloader while explicitly selecting the restored project library.
- `entrypoint.R` now shell-quotes independently parsed workflow arguments before forwarding them through `system2()`, so quoted values containing spaces are preserved.

## Candidate Work Items

- Continue to investigate the underlying `renv/activate.R` hang; the repository now has a documented noninteractive-safe runner as a workaround.
- Extend analysis-plan fixtures beyond the synthetic smoke cohort and validate their contract behavior with larger synthetic data when available.
- Extend `subgroup-forest-bundle` examples and tests to cover the GBM/LGG adjusted-effect and interaction fixture tables.
- The maintained `single-cohort-analysis-plan` workflow now provides the config-driven run-plan boundary for combined, subgroup, and interaction models without encoding GBM/LGG logic.

## Validation Targets

At minimum, the optimizing agent should be able to run:

```bash
./scripts/run_rscript_project.sh entrypoint.R subgroup-forest-bundle \
  --subgroup-summary-table test_data/gbm_lgg_p53beta_adjusted_effects_subgroup_summary.csv \
  --subgroup-type-filter survival_model \
  --output-dir results/validation_gbm_lgg_adjusted_forest
```

and:

```bash
./scripts/run_rscript_project.sh entrypoint.R subgroup-forest-bundle \
  --subgroup-summary-table test_data/gbm_lgg_p53beta_interactions_subgroup_summary.csv \
  --subgroup-type-filter interaction_model \
  --output-dir results/validation_gbm_lgg_interaction_forest
```

Both should complete without hanging and should write figure, table, and run metadata outputs.
