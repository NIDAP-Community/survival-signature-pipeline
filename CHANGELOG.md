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
- Updated README and architecture documentation to reflect the maintained versus legacy split, the reporting utility split, and the current validation baselines.

### Notes

- The repository is now runnable end to end through the maintained single-cohort smoke path and through pooled reporting contract smoke tests.
- Broader cohort generalization and legacy driver unification remain deferred.
