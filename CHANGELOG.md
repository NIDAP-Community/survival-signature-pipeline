# 2026-10-24

- repaired portable container builds behind an HTTPS-inspecting network by accepting an optional BuildKit CA secret without embedding institutional trust material in the image or repository
- added missing explicit Linux build dependencies and retained the renv package cache across interrupted local builds
- restored the locked R environment into a self-contained project library rather than creating links to an ephemeral build cache
- made the Docker runtime helpers portable to macOS and Linux, and prevented noninteractive survival plots from opening implicit `Rplots.pdf` devices
- built the native image and passed the canonical container smoke test; expected tiny-fixture Cox convergence and upstream ggpubr deprecation warnings remain non-fatal
- aligned the container architecture and CI scope with multi-gene-correlation: linux/amd64 build and smoke validation only
- published the validated linux/amd64 OCI image to GHCR by immutable digest; the CI workflow remains validation-only and does not hold registry publication credentials
- added a generic Apptainer/Singularity runtime wrapper that preserves the documented `/data` read-only and `/output` writable mount contract

# 2026-10-23

- changed maintained single-cohort defaults to low/mid/high tertiles for Kaplan-Meier displays and all eligible samples for continuous signature Cox models
- updated the reusable core function to mirror the maintained CLI defaults for programmatic callers
- improved reusable KM styling with accessible group colors, censoring marks, group sizes, and compact continuous/grouped model summaries
- added full adjusted-model forest reporting with the signature contrasts and every retained covariate coefficient across completed endpoints
- refreshed the README demonstration using a public METABRIC PAM50 proliferation-module example in the expected review order: OS KM, PFI KM, endpoint table, overall forest, then full adjusted-model forest
- re-ran the maintained single-cohort and reporting-bundle smoke paths with the project-local runner; expected tiny-fixture Cox convergence and upstream ggpubr deprecation warnings remain non-fatal

# 2026-04-29

- added explicit signature score aggregation (`sum` or `mean`) and continuous-model cohort selection (`grouped` or `all`) to the maintained single-cohort CLI
- fixed signature plotting when adjusted covariates are supplied and made degenerate quartile/tertile gene groupings non-fatal
- revalidated local pooled BRCA/METABRIC reproduction against its historical result table and single-cohort signature modeling against independent prostate Cox fits
- added a portable OCI container build, `/data` to `/output` runtime contract, Docker helpers, and a canonical container smoke test for integration readiness
- added non-publishing CI validation for a native container smoke test and a separate ARM64 image build
- added `single-cohort-analysis-plan`, a generic CSV-driven maintained workflow for overall, filtered-subgroup, and feature-by-group interaction Cox models
- added a documented analysis-plan schema with explicit optional endpoint mappings, signature definitions, covariates, filters, and interaction references
- emit analysis-plan primary effects and interaction terms through the shared subgroup-summary contract, with per-row status and full-coefficient outputs
- added `scripts/run_rscript_project.sh`, a noninteractive-safe local runner that bypasses a hanging renv autoloader while requiring the restored project library rather than silently using global packages
- documented `Rscript --vanilla` restoration and the project runner for reliable command-line execution
- shell-quoted root-entrypoint child arguments so quoted workflow values containing spaces survive `system2()` forwarding
- added a maintained `detectable` split mode to the single-cohort workflow so individual genes can be analyzed as true detectable-versus-undetectable groups using an explicit `detection_threshold`
- updated maintained single-cohort KM legends, comparison labels, and exported sample-group names so `split_method=detectable` renders as `Undetectable` versus `Detectable` rather than inheriting low/high labels
- documented the maintained detectability mode in the README and restricted it to individual gene features rather than `Signature`
- added a maintained `grouped_gene_cox_results.csv` export for single-cohort endpoints so grouped per-gene comparisons can be consumed without scraping plots or RDS objects, including explicit reference/comparison labels and retained covariates
- updated detectability mode so KM statistics rows use the full `Detectable vs Undetectable` wording and the maintained driver accepts `Undetectable,Detectable` as an `include-groups` alias
- added a maintained `adjusted_model_full_coefficients.csv` export for single-cohort endpoints so full adjusted Cox coefficient tables, including retained covariate terms, are available without extracting them from RDS plot objects

# 2026-04-13

- generalized the maintained single-cohort driver to support a generic `group_column` and `group_filter` interface while retaining `subtype_filter` as a compatibility path
- added a maintained `subgroup-forest-bundle` workflow that builds generic forest plots from the shared subgroup-summary contract
- added reusable subgroup forest plotting helpers so downstream reporting no longer needs to assume BRCA subtype semantics for every subgroup-oriented figure
- updated README examples to describe generic subgroup cohorting such as subtype, response class, or Gleason risk bands
- validated the generic subgroup filter path with a BRCA Basal smoke run and validated the new subgroup forest workflow against the shared subgroup-summary fixture

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
- Generalized the maintained pooled modeling internals so dataset interaction term lookup and implied dataset-window contrasts are resolved from the fitted model coefficients instead of assuming a METABRIC comparison level.
- Added optional dataset1_label and dataset2_label flags to the maintained pooled drivers so pooled outputs can carry dataset names that match the supplied inputs.
- Removed the maintained primary-only sample filtering step so maintained workflows now expect already post-filtered inputs from upstream preprocessing.
- Updated public-facing documentation and smoke-test examples to document the neutral two-dataset interface first and to use placeholder feature names GENE_A through GENE_D.

### Notes

- The maintained pooled workflows are now CLI-agnostic and reporting-agnostic with respect to dataset names, but the underlying modeling layer still retains a two-dataset comparison design and some historical term assumptions.
- The neutral pooled CLI surface was revalidated through the maintained dataset-interaction and subtype-specific smoke paths.
- The pooled modeling layer was revalidated with custom labels CohortA and CohortB to confirm that maintained outputs no longer depend on TCGA/METABRIC-specific coefficient names.
- Maintained workflows no longer attempt to infer or enforce a primary-only sample subset; upstream preprocessing is now the expected place for that decision.

## 2026-04-13

### Changed

- Updated the maintained single-cohort workflow defaults so tertile analyses keep `low,mid,high` by default and plot generation is enabled by default for the maintained path.
- Fixed reusable KM plotting so legend labels and colors are derived from fitted strata rather than assumed factor levels.
- Fixed maintained per-gene and early/late adjusted KM summary generation by using explicit complete-case adjusted model data and qualifying constructed formulas with `survival::Surv(...)`.
- Refactored the maintained single-cohort covariate interface so `--adjust-covariates` accepts explicit raw metadata column names, defaults to `none`, and skips missing requested columns with a warning.
- Kept shorthand aliases for `subtype`, `stage`, `age`, and `race` only as convenience resolution on top of the explicit-column interface.
- Updated README guidance to document explicit dataset-specific covariate selection as the preferred pattern for reusable single-cohort runs.

### Notes

- BRCA validation confirmed that the maintained AMD1 OS early-window adjusted KM rows now render with the expected `N = 956` basis.
- Residual differences between the earlier legacy BRCA figure and the maintained BRCA figure were traced to a different adjustment set, specifically whether `race_collapsed` was included.
- Broader generalization work has been identified around endpoint configuration, preprocessing contracts, and shared maintained covariate resolution across all maintained drivers.
