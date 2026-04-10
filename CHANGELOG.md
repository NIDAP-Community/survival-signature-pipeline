## 2026-04-08

### Added

- Added a maintained root CLI in entrypoint.R.
- Added a maintained single-cohort survival workflow.
- Added a maintained single-cohort reporting workflow.
- Added a shared reporting utility layer in src/utilities/reporting_helpers.R.
- Added project-local renv activation files and renv.lock.
- Added canonical synthetic smoke-test fixtures and subgroup-summary fixtures under test_data/.

### Changed

- Standardized maintained workflows on run-scoped outputs under results/.
- Made subtype optional in the single-cohort workflow unless subtype-dependent options are requested.
- Standardized maintained downstream reporting on one shared subgroup-summary contract.
- Updated the pooled reporting workflow to validate its two-dataset subtype-stratified contract explicitly.
- Consolidated shared reporting output directory creation, CSV plus Markdown table writing, and text-block writing into shared helpers.
- Consolidated maintained-driver CLI parsing, help rendering, repeated parser behavior, endpoint mapping, required input-file checks, and maintained analysis-table loading into src/utilities/cli_helpers.R.
- Updated README and architecture documentation to reflect the maintained versus legacy split, the reporting utility split, and the current validation baselines.

### Notes

- The repository is now runnable end to end through the maintained single-cohort smoke path and through pooled reporting contract smoke tests.
- Broader cohort generalization and legacy driver unification remain deferred.

## 2026-04-10

### Changed

- Generalized the maintained pooled reporting bundle so downstream dataset labels, panel ordering, exploratory dataset selection, and summary text are inferred from supplied result tables instead of hardcoding TCGA and METABRIC.
- Consolidated shared reporting table-formatting helpers into src/utilities/reporting_helpers.R and updated maintained reporting drivers to use them.
- Replaced the maintained pooled CLI surface with neutral dataset1_* and dataset2_* flag names while retaining tcga_* and metabric_* names as compatibility aliases.
- Updated public-facing documentation and smoke-test examples to document the neutral two-dataset interface first and to use placeholder feature names GENE_A through GENE_D.

### Notes

- The maintained pooled workflows are now CLI-agnostic and reporting-agnostic with respect to dataset names, but the underlying modeling layer still retains a two-dataset comparison design and some historical term assumptions.
- The neutral pooled CLI surface was revalidated through the maintained dataset-interaction and subtype-specific smoke paths.
