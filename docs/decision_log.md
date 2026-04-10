# Decision Log

## 2026-04-08: Standardize Around One Maintained Root CLI

Decision:

- Use entrypoint.R as the maintained root command surface for repository workflows.

Rationale:

- The repository needs one stable non-interactive Rscript entrypoint rather than multiple competing top-level invocation patterns.
- This matches the intended architecture and makes future workflow additions easier to document and test.

Current effect:

- Maintained workflows are dispatched through entrypoint.R.

## 2026-04-08: Treat drivers/ as the Maintained Workflow Layer

Decision:

- Treat drivers/ as the maintained workflow layer and src/drivers/ as legacy-style wrappers until a later unification step.

Rationale:

- The repository currently has a mixed structure, but the newer extraction and standardization work is happening in drivers/.
- Forcing an immediate directory reorganization would add churn before the workflow contracts are stable.

Current effect:

- New maintained workflows were added under drivers/.
- Legacy HTML-oriented wrappers remain in src/drivers/ and are documented as such.

## 2026-04-08: Default Maintained Outputs to Run-Scoped results/ Directories

Decision:

- Maintained workflows should default to run-scoped outputs under results/ rather than legacy Output_* or dated ad hoc locations.

Rationale:

- Reproducibility and provenance are easier to preserve when runs are isolated.
- This reduces silent overwrites and aligns with the repository’s explicit provenance goal.

Current effect:

- Shared helpers in src/utilities/run_provenance_guardrails.R now generate run-scoped output paths.

## 2026-04-08: Keep Single-Cohort Analysis Separate From Pooled Two-Dataset Workflows

Decision:

- Maintain a separate single-cohort analysis workflow rather than forcing one-dataset use into pooled TCGA plus METABRIC drivers.

Rationale:

- The reusable core can run on one dataset, but the maintained pooled drivers still encode two-dataset assumptions.
- A separate single-cohort path is cleaner than weakening the meaning of the pooled workflow contract.

Current effect:

- drivers/run_single_cohort_survival.R is the maintained one-dataset analysis path.

## 2026-04-08: Make Subtype Optional in Single-Cohort Analysis Unless Explicitly Required

Decision:

- Treat subtype as optional in the maintained single-cohort workflow unless subtype filtering or subtype-adjusted modeling is requested.

Rationale:

- Subtype is not inherently required by the reusable core.
- Requiring subtype unconditionally would block valid one-dataset analyses with simpler metadata.

Current effect:

- The single-cohort driver now fails clearly only when subtype-dependent options are requested without a subtype column.

## 2026-04-08: Keep Pooled Reporting and Single-Cohort Reporting as Separate Modes

Decision:

- Do not overload the finalized pooled reporting bundle to handle one-dataset or no-subtype cases.
- Use a separate maintained single-cohort reporting workflow instead.

Rationale:

- The pooled reporting bundle is semantically tied to two-dataset pooled inputs and subtype-stratified reporting.
- Separate reporting modes keep the downstream contracts clear and avoid ambiguous behavior.

Current effect:

- drivers/run_finalized_pooled_results_bundle.R now validates and rejects inputs that violate its pooled contract.
- drivers/run_single_cohort_reporting_bundle.R handles single-cohort downstream reporting.

## 2026-04-08: Standardize One Shared Subgroup-Summary Contract for Downstream Heatmaps

Decision:

- Use a shared subgroup-summary schema for downstream heatmap-style reporting.

Required core columns:

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

Rationale:

- Downstream reporting should not depend on driver-specific table names or bespoke column layouts.
- A shared subgroup contract creates a stable handoff between subgroup-producing workflows and downstream reporting workflows.

Current effect:

- drivers/run_piecewise_subtype_specific.R emits subgroup_effect_summary.csv in the shared shape.
- drivers/run_single_cohort_reporting_bundle.R consumes and normalizes subgroup tables using that same shape.

## 2026-04-08: Use Synthetic Fixtures for Interface Validation, Not Statistical Validation

Decision:

- Keep small synthetic test_data fixtures for smoke testing interface and contract behavior.

Rationale:

- The immediate goal is validating workflow structure, CLI contracts, and output expectations.
- Tiny synthetic cohorts are acceptable for smoke testing even when Cox models emit convergence warnings.

Current effect:

- Smoke tests validate dispatch, guardrails, output generation, and subgroup-contract wiring.

## 2026-04-08: Standardize Both Maintained Reporting Workflows on One Shared Subgroup-Summary Contract

Decision:

- Use one shared downstream subgroup-summary contract for both maintained reporting workflows.
- Keep the pooled reporting workflow's legacy subgroup-table flag only as a compatibility alias.

Rationale:

- Downstream reporting should not depend on workflow-specific subgroup table names or bespoke column layouts.
- One normalized contract makes subgroup-producing workflows and reporting consumers easier to validate and document.

Current effect:

- drivers/run_single_cohort_reporting_bundle.R and drivers/run_finalized_pooled_results_bundle.R both normalize and validate subgroup inputs through shared helpers.

## 2026-04-08: Add a Small Shared Reporting Utility Layer Instead of Merging Reporting Workflows

Decision:

- Extract only narrowly shared reporting helpers into src/utilities/reporting_helpers.R.
- Keep single-cohort and pooled reporting as separate maintained workflows.

Rationale:

- The reporting workflows overlap in plotting preparation and output scaffolding, but they still have distinct semantics, layouts, and output scopes.
- A small utility layer reduces duplication without collapsing separate reporting modes into one ambiguous driver.

Current effect:

- src/utilities/reporting_helpers.R now holds shared subgroup heatmap preparation, common plotting defaults, and shared reporting output helpers.

## 2026-04-08: Keep Reporting CSV Contracts Separate From Markdown Presentation

Decision:

- Preserve raw machine-readable CSV outputs separately from formatted Markdown table outputs when extracting shared reporting writers.

Rationale:

- Reporting refactors should not silently change downstream table contracts.
- Markdown presentation formatting and CSV data export serve different consumers and should not be conflated.

Current effect:

- Shared reporting table writers accept separate CSV and Markdown data frames where needed so existing output contracts remain intact.

## 2026-04-08: Standardize Project-Local Reproducibility with renv

Decision:

- Check in the renv activation files and lockfile and treat renv as the standard local execution environment.

Rationale:

- The repository should be runnable and reproducible before broader refactoring.
- A project-local environment makes the maintained workflows easier to validate consistently.

Current effect:

- .Rprofile, renv/, and renv.lock are present and documented as part of the maintained execution path.

## 2026-04-08: Keep Shared Maintained-Driver Utilities Small and Pre-Workflow

Decision:

- Centralize only narrowly repeated maintained-driver mechanics in src/utilities/cli_helpers.R.
- Keep dataset-specific preparation, metadata validation, endpoint-column checks inside each workflow, cohort-aware semantics, and modeling logic local.

Rationale:

- The maintained drivers benefit from less duplication in CLI parsing and pre-load mechanics, but a broad input abstraction layer would hide workflow-specific contracts too early.
- Small shared helpers are easier to validate against the canonical smoke baselines.

Current effect:

- src/utilities/cli_helpers.R now holds shared maintained-driver CLI, parsing, endpoint-mapping, required-file-check, and table-loading helpers without changing workflow-specific modeling behavior.

## 2026-04-10: Prefer Neutral Dataset Names in the Maintained Pooled CLI

Decision:

- Use dataset1_* and dataset2_* as the primary maintained flag names for pooled workflows.
- Retain tcga_* and metabric_* only as compatibility aliases.

Rationale:

- The maintained pooled workflows are intended to support any pair of similarly structured datasets, not only historically named cohorts.
- Neutral flag names reduce unnecessary project-specific coupling at the CLI surface without forcing a deeper model rewrite in the same change.

Current effect:

- drivers/run_piecewise_dataset_interaction.R and drivers/run_piecewise_subtype_specific.R now document and parse neutral dataset1/dataset2 flags first.
- Historical tcga/metabric flag names still resolve to the same maintained options for backward compatibility.

## 2026-04-10: Infer Pooled Reporting Dataset Labels From Inputs

Decision:

- Make the maintained pooled reporting bundle derive dataset labels and exploratory dataset defaults from the supplied upstream result tables rather than hardcoding named cohorts in reporting output.

Rationale:

- The downstream reporting layer does not need hardcoded cohort names when the upstream tables already carry dataset labels.
- This keeps reporting reusable while leaving the deeper pooled modeling semantics unchanged until a later refactor.

Current effect:

- drivers/run_finalized_pooled_results_bundle.R now orders reporting panels from observed dataset and window values, resolves exploratory_dataset dynamically when omitted, and writes generic report text instead of TCGA/METABRIC-specific narrative.
