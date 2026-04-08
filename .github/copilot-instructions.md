# Copilot Instructions — Ready To Copy Into The New Repo

Copy this file into the new repository as:

```text
.github/copilot-instructions.md
```

---

# Copilot Instructions — Reusable Survival Analysis Repo

This repository contains reusable R-based survival-analysis code.

Its architecture should remain aligned to the style used in `NIDAP-Community/multi-gene-correlations`, adapted for survival analysis rather than correlation analysis.

## Primary Goal

Favor reusable survival-analysis code over historical compatibility.

The repository should support similar expression + survival datasets without being locked to the original project structure, cohorts, or manuscript outputs.

## Architectural Standard

Prefer a structure like this:

- `.github/`
- `container/`
- `docs/`
- `scripts/`
- `src/core/`
- `src/drivers/`
- `src/utilities/`
- `test_data/`
- `results/`
- `renv/`
- `renv.lock`
- `entrypoint.R`
- `README.md`
- `CHANGELOG.md`

Keep the repo visually and operationally consistent with `NIDAP-Community/multi-gene-correlations` where practical.

## Runtime Model

- All executable workflows must run non-interactively via `Rscript`.
- Prefer one stable root `entrypoint.R` as the main CLI surface.
- Use `scripts/` only as convenience wrappers.
- Support local execution via `renv`.
- Keep the repo container-ready, even before containerization is complete.

## Input Contract Rules

Treat input contracts as explicit and documented.

At minimum, most workflows will require:
- an expression matrix with a gene identifier column
- a metadata table with a sample identifier column
- event and survival time columns for each endpoint

Where practical, these should be configurable:
- sample ID column
- gene ID column
- subtype column
- stage column
- age column
- endpoint column names

## Reusability Rules

- Prefer CLI arguments or config over hard-coded paths.
- Remove hard-coded cohort names unless logically required.
- Isolate cohort-specific filtering rules so they can be disabled or replaced.
- Avoid user-specific paths.
- Keep manuscript/reporting logic separate from the core analysis engine.

## Provenance Rules

- Outputs should write to run-scoped directories under `results/`.
- Do not silently overwrite prior outputs.
- Preserve or improve run metadata generation.
- Record inputs, key parameters, driver path, and session/environment details when feasible.

## Layer Separation

Keep these layers distinct:

### Core analysis layer
- reusable survival functions
- data preparation helpers
- signature scoring and model helpers

### Maintained driver layer
- overall survival entrypoints
- subtype-specific survival entrypoints
- pooled piecewise modeling entrypoints

### Optional reporting layer
- pooled forest plot assembly
- subtype heatmap assembly
- manuscript-oriented summary tables or prose

Reporting code must not become a hidden dependency of the core analysis engine.

## Documentation Rules

When changing interfaces, update `README.md` with:
- required inputs
- expected columns
- example commands
- output layout
- any remaining dataset-specific assumptions

## Preferred First Tasks

When starting work in this repo, begin by:

1. mapping dependencies across the copied scripts
2. identifying hard-coded assumptions
3. checking how each script fits the standardized architecture
4. separating reusable analysis code from optional reporting code
5. updating README documentation before large refactors

## What To Avoid

- Do not default to preserving legacy output conventions.
- Do not assume all datasets behave like TCGA or METABRIC.
- Do not bury important assumptions in helper functions without documenting them.
- Do not mix one-off manuscript writing logic into reusable analysis code.

## Operating Principle

Choose the implementation that is more reusable, more explicit, and more reproducible, even if it is less convenient than a project-specific shortcut.