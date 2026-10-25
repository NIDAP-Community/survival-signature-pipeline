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
- src/utilities/reporting_helpers.R contains shared reporting-side subgroup heatmap preparation, shared reporting output scaffolding, shared reporting table-formatting helpers, and minimal common plotting defaults.
- src/utilities/cli_helpers.R contains the small shared maintained-driver CLI and pre-load utility layer.
- drivers/ contains the maintained pooled and piecewise workflows.
- drivers/ now also contains maintained single-cohort expression/signature and binary-biomarker survival workflows that run through the reusable core.
- drivers/ contains a maintained declarative single-cohort analysis-plan workflow for generic adjusted, subgroup, and interaction Cox models.
- drivers/ also contains maintained downstream reporting workflows for single-cohort and pooled reporting bundles.
- src/drivers/ contains older HTML-oriented and report-oriented wrappers that still depend on the reusable core.
- entrypoint.R is the canonical root CLI surface and dispatches to maintained workflows.
- renv/ plus renv.lock are now present and the repository can be restored into a project-local R environment.
- test_data/ now includes one validated canonical smoke-test input pair for the single-cohort workflow.
- container/ now contains a portable OCI build, Docker mount helpers, and a canonical container smoke test for the maintained root CLI.

This means the repository is now structurally standardized enough to have a single runnable entrypoint, while still preserving the distinction between maintained workflows and older wrappers.

The container is deliberately a portable batch-analysis boundary, not an AWS application. It has no database or cloud credentials, consumes mounted `/data` inputs, and writes mounted `/output` artifacts. Future AWS/OncoChat work should stage and validate data outside the image rather than placing integration logic in the reusable analysis repository.

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

The reporting utility layer currently lives in src/utilities/reporting_helpers.R. This layer is responsible for narrowly shared downstream reporting helpers, currently limited to subgroup heatmap preparation, shared reporting output scaffolding, shared reporting table-formatting helpers, and a small set of common plotting defaults.

This layer exists to reduce duplicated reporting code without merging separate reporting workflows or pushing reporting behavior into the core analysis engine.

### Maintained CLI utility layer

The maintained CLI utility layer currently lives in src/utilities/cli_helpers.R. This layer is responsible for narrowly shared pre-analysis and CLI mechanics used by the maintained drivers, including:

- standardized help rendering with entrypoint-first usage
- shared `--flag value` parsing and required-argument enforcement
- maintained compatibility alias handling
- repeated CSV-list, numeric, and enum-style flag parsing
- shared maintained OS/PFI endpoint mapping
- minimal required analysis-input file checks
- shared maintained analysis-table loading for expression and metadata inputs

This layer is intentionally small. It does not replace workflow-specific preparation. Dataset-specific filters, metadata validation, endpoint-column checks inside each workflow, cohort-aware assumptions, and all modeling logic remain local to the individual maintained drivers.

### Maintained driver and reporting layer

The maintained driver layer is currently the drivers/ directory. These scripts should:

- run non-interactively via Rscript
- accept explicit CLI arguments
- write to run-scoped directories under results/
- record metadata through the provenance utility
- stay focused on workflow execution rather than manuscript text

They now rely on a small shared CLI utility layer for repeated surface-level mechanics, while still retaining workflow-specific semantics. The maintained pooled CLI now prefers neutral dataset1/dataset2 flag names, with legacy tcga/metabric aliases retained for compatibility, and the pooled maintained drivers can emit dataset labels from dataset1_label/dataset2_label rather than hardcoding cohort names in their maintained outputs. The maintained drivers also assume datasets are already post-filtered upstream rather than reapplying a primary-only sample filter internally.

In the first wave, these drivers still retain some historical two-dataset naming and model assumptions. The second wave will reduce those assumptions further.

The maintained driver layer now has two categories:

- pooled two-dataset workflows for independent or comparable cohort comparisons
- single-cohort expression/signature and binary-biomarker workflows for one prepared dataset at a time, using the reusable core directly
- a declarative single-cohort analysis-plan workflow for a prepared cohort with explicit covariates, filters, subgroup effects, and interaction tests

Biologically distinct disease-state comparisons, such as GBM versus LGG, should not be forced into the pooled two-dataset workflow by default. Those analyses need a single-cohort or pan-disease-style design with explicit covariates, subgroup analyses, and interaction tests for major biological axes such as IDH mutation status.

Within drivers/, reporting is still a separate downstream concern:

- run_single_cohort_reporting_bundle.R is the maintained downstream reporting path for single-cohort outputs
- run_finalized_pooled_results_bundle.R is the maintained downstream reporting path for pooled two-dataset outputs

These maintained reporting workflows share some utility code, but they still remain separate workflows with separate reporting scopes, layouts, and outputs.

An intentionally retained inconsistency is the compatibility alias layer and the fixed two-dataset comparison design in the pooled maintained workflows. Legacy flags such as tcga_expr_file and metabric_meta_file are still accepted for compatibility, even though the primary maintained interface now uses dataset1_expr_file and dataset2_meta_file and the maintained outputs can label datasets generically.

### Legacy layer

The src/drivers/ directory is currently best treated as a legacy wrapper layer. These scripts still provide useful historical behavior, but they mix execution and HTML-style reporting more tightly than the long-term architecture should.

The legacy layer still applies, but it is no longer the same thing as the maintained reporting layer in drivers/.

## First-Wave Decisions

The current cleanup wave follows these rules:

1. keep one clean root entrypoint before larger refactors
2. parameterize maintained drivers enough to remove hard-coded run dates and code-embedded column assumptions
3. preserve the current two-dataset comparison mechanics for now instead of forcing early over-generalization
4. add smoke-test placeholders in test_data/ so the repository has an explicit input contract
5. add a canonical smoke-test dataset in test_data/ so one maintained workflow is runnable end to end
6. defer bigger cohort-generalization and driver unification work to the second wave

## Known Gaps

- container/ is still a placeholder
- scripts/ is still empty
- src/drivers/ has not been moved or archived yet
- the maintained pooled summary workflows still require explicit upstream result tables
- planned analyses currently support continuous gene and mean-expression signature effects; grouped KM-style effects remain the responsibility of `single-cohort-survival`
- broader cohort generalization should avoid overloading two-dataset interaction workflows for biologically distinct disease groups such as GBM versus LGG

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

The maintained-driver utility cleanup was also checkpointed by rerunning both baselines after the CLI/help, parser, endpoint, file-check, and loader refactors. Those reruns used separate output directories and confirmed that the shared utility extractions did not alter the canonical baseline contracts.
