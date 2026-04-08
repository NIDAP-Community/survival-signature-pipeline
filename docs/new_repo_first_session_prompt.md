# New Repo First Session Prompt

Paste this into the very first Copilot chat in the new repository.

```text
This repository is a reusable extraction of an older R-based survival-analysis workflow.

The repo should follow the same broad architectural style as NIDAP-Community/multi-gene-correlations, but for survival analysis.

Please work from these priorities:
1. keep the repo reusable and standardized
2. preserve Rscript-based execution
3. preserve provenance and run metadata
4. separate core analysis code from optional reporting/manuscript code
5. minimize hard-coded dataset assumptions

Please start by doing the following:
1. inspect the repo structure and map dependencies
2. identify all hard-coded paths, cohort names, and column-name assumptions
3. tell me whether the current files fit the intended architecture:
   - root entrypoint
   - container/
   - scripts/
   - src/core/
   - src/drivers/
   - src/utilities/
   - test_data/
   - results/
   - docs/
   - renv.lock
4. propose the smallest safe cleanup plan
5. update or draft README documentation before large refactors

Assume I want a reusable analysis repo, not a historical archive.
```