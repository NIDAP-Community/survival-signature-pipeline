# Repository Architecture

This repository is being standardized as a reusable R-based survival-analysis pipeline. The immediate goal is not to finish the full cohort-generalization refactor; the first wave is to make the repository structurally coherent, runnable through one root entrypoint, and explicit about which code belongs to which layer.

## Target Shape

The intended long-term structure is:

- entrypoint.R as the main non-interactive CLI surface
- container/ for containerization assets
- scripts/ for convenience wrappers only
- src/core/ for reusable analysis functions
- src/drivers/ or drivers/ for maintained workflow entrypoints
- src/utilities/ for shared helpers such as provenance, path guardrails, and CLI support
- test_data/ for small synthetic or placeholder inputs
- results/ for run-scoped outputs
- docs/ for architecture notes and decisions
- renv/ plus renv.lock for reproducible environments

## Current State

The repository currently has one strong reusable core and one strong reusable utility, but the driver layer is split.

- src/core/GeneExpressionBasedSurvivalAnalys.r contains the main reusable survival-analysis function.
- src/utilities/run_provenance_guardrails.R contains shared output and metadata helpers.
- drivers/ contains the maintained pooled and piecewise workflows.
- drivers/ now also contains a maintained single-cohort workflow that runs through the reusable core.
- drivers/ now also contains a maintained single-cohort reporting bundle for downstream overall figures and optional subgroup heatmaps.
- src/drivers/ contains older HTML-oriented and report-oriented wrappers that still depend on the reusable core.
- entrypoint.R is the canonical root CLI surface and dispatches to maintained workflows.

This means the repository is now structurally standardized enough to have a single runnable entrypoint, while still preserving the distinction between maintained workflows and older wrappers.

## Layer Boundaries

### Core layer

The core layer should contain reusable analysis logic that does not depend on manuscript formatting or one-off reporting conventions. Today that means src/core/GeneExpressionBasedSurvivalAnalys.r is the main reusable engine and should remain the place where durable statistical functionality lives.

### Maintained driver layer

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

The reporting layer now also has a single-cohort downstream bundle that is intentionally separate from the pooled two-dataset reporting workflow.

### Legacy and reporting layer

The src/drivers/ directory is currently best treated as a legacy wrapper layer. These scripts still provide useful historical behavior, but they mix execution and HTML-style reporting more tightly than the long-term architecture should.

The pooled forest and heatmap assembly in drivers/run_finalized_pooled_results_bundle.R is also a reporting-oriented layer. It remains maintained, but should stay downstream of analysis drivers rather than becoming a hidden dependency of the core engine.

### Utility layer

Shared helpers belong in src/utilities/. Provenance guardrails already live there, and additional shared CLI or run-directory helpers should continue to accumulate there rather than being copied across drivers.

## First-Wave Decisions

The current cleanup wave follows these rules:

1. keep one clean root entrypoint before larger refactors
2. parameterize maintained drivers enough to remove hard-coded run dates and code-embedded column assumptions
3. preserve TCGA and METABRIC awareness for now instead of forcing early over-generalization
4. add smoke-test placeholders in test_data/ so the repository has an explicit input contract
5. defer bigger cohort-generalization and driver unification work to the second wave

## Known Gaps

- renv.lock is still missing
- container/ is still a placeholder
- scripts/ is still empty
- src/drivers/ has not been moved or archived yet
- the maintained pooled summary workflows still require explicit upstream result tables
- test_data/ is intended for smoke testing and interface validation, not scientific inference

## Expected Flow

The recommended execution flow is:

1. use entrypoint.R to dispatch a maintained workflow
2. point the workflow at raw expression and metadata inputs, or at test_data placeholders for interface checks
3. write outputs into a run-scoped directory under results/
4. use the downstream summary and pooled-reporting drivers only after upstream result tables exist

This keeps the repository aligned with reusable workflow execution rather than historical archive behavior.
