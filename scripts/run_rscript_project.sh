#!/usr/bin/env bash
# Run an R script with this repository's restored renv library, without
# invoking the project's noninteractive renv autoloader. This is intended for
# local command-line use when `Rscript` startup through `.Rprofile` hangs.

set -euo pipefail

# Resolve the repository root relative to this script, so callers can invoke
# the runner from any working directory.
script_dir="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(CDPATH= cd -- "${script_dir}/.." && pwd)"

# renv stores the bootstrap package inside the platform-specific project
# library. Its parent is the library path required by R_LIBS_USER.
renv_package_dir="$(find "${project_dir}/renv/library" -maxdepth 4 -type d -name renv -print -quit 2>/dev/null || true)"
project_library="$(dirname "${renv_package_dir}")"

# A restored plotting stack is required by maintained KM workflows. Fail with
# the exact recovery command instead of silently falling back to incompatible
# globally installed packages.
if [[ -z "${renv_package_dir}" || ! -d "${project_library}/ggplot2" ]]; then
  printf '%s\n' "Project renv library is missing or incomplete." >&2
  printf '%s\n' "Restore it from the repository root with:" >&2
  printf '%s\n' "  Rscript --vanilla -e 'if (!requireNamespace(\"renv\", quietly = TRUE)) install.packages(\"renv\", repos = \"https://cloud.r-project.org\"); renv::restore(project = getwd(), prompt = FALSE)'" >&2
  exit 1
fi

# Keep the R process noninteractive-safe while preferring the restored project
# packages. Base R libraries remain available after this project library.
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export R_LIBS_USER="${project_library}${R_LIBS_USER:+:${R_LIBS_USER}}"

exec Rscript "$@"
