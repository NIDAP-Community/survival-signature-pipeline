# Repository Architecture

This repository is being standardized as a reusable R-based survival-analysis pipeline. The immediate goal is not to finish the full cohort-generalization refactor; the first wave is to make the repository structurally coherent, runnable through one root entrypoint, and explicit about which code belongs to which layer.

## Target Shape

The intended long-term structure is:

- entrypoint.R as the main non-interactive CLI surface
- container/ for containerization assets
- scripts/ for convenience wrappers only
- src/core/ for reusable analysis functions
- drivers/ for maintained workflow and reporting entrypoints
- src/drivers/ for legacy wrappers until later unification or retirement
- src/utilities/ for shared helpers such as provenance, path guardrails, and CLI support
- test_data/ for small smoke-test inputs and contract fixtures
- results/ for run-scoped outputs
- docs/ for architecture notes and decisions
- renv/ plus renv.lock for reproducible environments

## Current State

The repository currently has one strong reusable core, a small but real utility split, and a still-split driver layer.

- src/core/GeneExpressionBasedSurvivalAnalys.r contains the main reusable survival-analysis function.
- src/utilities/run_provenance_guardrails.R contains shared provenance, output, and subgroup-summary interface helpers.
- src/utilities/reporting_helpers.R contains shared reporting-side subgroup heatmap preparation and minimal common plotting defaults.
- drivers/ contains the maintained pooled and piecewise workflows.
- drivers/ now also contains a maintained single-cohort analysis workflow that runs through the reusable core.
- drivers/ also contains maintained downstream reporting workflows for single-cohort and pooled reporting bundles.
- src/drivers/ contains older HTML-oriented and report-oriented wrappers that still depend on the reusable core.
- entrypoint.R is the canonical root CLI surface and dispatches to maintained workflows.
- renv/ plus renv.lock are now present and the repository can be restored into a project-local R environment.
- test_data/ now includes one validated canonical smoke-test input pair for the single-cohort workflow.

This means the repository is now structurally standardized enough to have a single runnable entrypoint, while still preserving the distinction between maintained workflows and older wrappers.

## Layer Boundaries

### Core layer

The core layer should contain reusable analysis logic that does not depend on manuscript formatting or one-off reporting conventions. Today that means src/core/GeneExpressionBasedSurvivalAnalys.r is the main reusable engine and should remain the place where durable statistical functionality lives.

### Provenance utility layer

The provenance utility layer currently lives in src/utilities/run_provenance_guardrails.R. This layer is responsible for:

- run-scoped output-path defaults
- output guardrails
- run metadata generation
- shared subgroup-summary normalization and validation helpers used by downstream reporting workflows

These helpers are shared infrastructure and should remain independent from analysis-specific or manuscript-specific code.

### Reporting utility layer

The reporting utility layer currently lives in src/utilities/reporting_helpers.R. This layer is responsible for narrowly shared downstream reporting helpers, currently limited to subgroup heatmap preparation and a small set of common plotting defaults.

This layer exists to reduce duplicated reporting code without merging separate reporting workflows or pushing reporting behavior into the core analysis engine.

### Maintained driver and reporting layer

The maintained driver layer is currently the drivers/ directory. These scripts should:

- run non-interactively via Rscript
- accept explicit CLI arguments
- write to run-scoped directories under results/
- record metadata through the provenance utility
- stay focused on workflow execution rather than manuscript text

In the first wave, these drivers remain TCGA and METABRIC aware. The second wave will reduce those assumptions further.

The maintained driver layer now has two categories:

- pooled two-dataset workflows for TCGA plus METABRIC comparisons
- a single-cohort workflow for one dataset at a time, using the reusable core directly

Within drivers/, reporting is still a separate downstream concern:

- run_single_cohort_reporting_bundle.R is the maintained downstream reporting path for single-cohort outputs
- run_finalized_pooled_results_bundle.R is the maintained downstream reporting path for pooled two-dataset outputs

These maintained reporting workflows share some utility code, but they still remain separate workflows with separate reporting scopes, layouts, and outputs.

### Legacy layer

The src/drivers/ directory is currently best treated as a legacy wrapper layer. These scripts still provide useful historical behavior, but they mix execution and HTML-style reporting more tightly than the long-term architecture should.

The legacy layer still applies, but it is no longer the same thing as the maintained reporting layer in drivers/.

## First-Wave Decisions

The current cleanup wave follows these rules:

1. keep one clean root entrypoint before larger refactors
2. parameterize maintained drivers enough to remove hard-coded run dates and code-embedded column assumptions
3. preserve TCGA and METABRIC awareness for now instead of forcing early over-generalization
4. add smoke-test placeholders in test_data/ so the repository has an explicit input contract
5. add a canonical smoke-test dataset in test_data/ so one maintained workflow is runnable end to end
6. defer bigger cohort-generalization and driver unification work to the second wave

## Known Gaps

- container/ is still a placeholder
- scripts/ is still empty
- src/drivers/ has not been moved or archived yet
- the maintained pooled summary workflows still require explicit upstream result tables
- broader cohort generalization has not started yet

## Expected Flow

The recommended execution flow is:

1. use entrypoint.R to dispatch a maintained workflow
2. restore the project environment with renv if needed
3. point the workflow at raw expression and metadata inputs, or at the canonical test_data smoke-test pair for interface checks
4. write outputs into a run-scoped directory under results/
5. use the downstream summary and pooled-reporting drivers only after upstream result tables exist

This keeps the repository aligned with reusable workflow execution rather than historical archive behavior.

## Current Validation Points

The repository currently has two concrete validation anchors:

1. The canonical maintained single-cohort smoke test through entrypoint.R, writing to results/smoke_test_single_cohort_canonical.
2. The pooled reporting contract smoke test through entrypoint.R, writing to results/pooled_reporting_contract_smoke_shared_plotting.

Together these validate the maintained single-cohort baseline, run-scoped provenance outputs, and the shared downstream subgroup-summary reporting interface.
