# GBM/LGG p53beta Survival Field Report

Date: 2026-05-21

## Purpose

This field report captures lessons from applying the reusable survival-signature-pipeline to a GBM/LGG p53beta signature project in a separate analysis repository. The goal is to preserve what the project taught us about general pipeline behavior, not to turn GBM/LGG into a pipeline-specific special case.

## External Analysis Context

- External project: `/Users/maggiec/GitHub/Maggie/Harris/GBM`
- WT signature result folder: `results/p53b_wt_signature_gbm_lgg_survival_pipeline`
- R273H signature result folder: `results/p53b_r273h_signature_gbm_lgg_survival_pipeline`
- Forest summary folder: `results/p53b_signature_gbm_lgg_survival_forest_panels`
- Signature inputs: precomputed GSVA scores for `p53b_WT_up_all_no_TP53` and `p53b_R273H_up_all_no_TP53`

## Scientific Design Lessons

- GBM and LGG should not be treated like two independent validation datasets analogous to TCGA-BRCA and METABRIC.
- GBM/LGG are biologically distinct glioma disease states within a related tumor setting.
- IDH mutation status is central biology and should be a major covariate, subgroup axis, and possible interaction axis.
- TP53 mutation status can be added as a covariate and subgroup axis, but it is entangled with IDH/LGG biology and must be interpreted carefully.
- AJCC or clinical stage is not available in this setting; disease group, histologic grade, or glioma subgroup must be labeled as proxies if used.

## Workflow Needs Exposed

- A generic single-cohort survival workflow should support combined models, subgroup-specific models, and signature-by-group interaction models in one auditable run plan.
- Subgroup survival effects and interaction effects should be reported separately because they answer different questions.
- The shared subgroup-summary contract worked well as a handoff layer for downstream forest plots.
- The pipeline needs a maintained aggregator that can convert a directory of pipeline run outputs into one subgroup-summary table.
- Interaction model outputs should be able to enter the same downstream reporting system as subgroup effects, but with clear labels that identify them as effect-modification tests.

## Technical Lessons

- The maintained `subgroup-forest-bundle` plotting helper is useful for this use case.
- Normal `Rscript` execution from the pipeline repo hung during `.Rprofile` startup because `.Rprofile` activates `renv`.
- `scripts/run_rscript_project.sh` is the maintained workaround after `renv::restore()`: it disables the autoloader while selecting the restored project library.
- `entrypoint.R` now shell-quotes independently parsed workflow arguments before forwarding them through `system2()`, preserving quoted labels containing spaces.

## Recommended Pipeline Follow-Up

- Continue investigating why `renv/activate.R` can hang in noninteractive `Rscript` runs on this system; a noninteractive-safe runner is now available.
- Add a maintained single-cohort stratified/interaction workflow that consumes an analysis configuration rather than requiring a bespoke wrapper script.
- Add a maintained result aggregator for adjusted Cox rows and interaction rows.
- Add regression tests using the GBM/LGG p53beta fixture tables added with this field report.

## Fixture Tables

Two fixture tables were added to `test_data/`:

- `test_data/gbm_lgg_p53beta_adjusted_effects_subgroup_summary.csv`
- `test_data/gbm_lgg_p53beta_interactions_subgroup_summary.csv`

These are real downstream subgroup-summary tables produced from the GBM/LGG p53beta analysis. They are intentionally small enough to use as reporting fixtures while still representing a real multi-signature, multi-endpoint, subgroup-plus-interaction case.

## Interpretation Reminder

Subgroup effect rows answer: "Is the signature associated with survival within this model or subgroup?"

Interaction rows answer: "Is the signature effect different in the named group compared with the reference group?"

These should not be collapsed into one biological interpretation without explicit labeling.
