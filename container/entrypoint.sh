#!/usr/bin/env bash
# Resolve the project-local renv library inside the image, then dispatch the
# repository's maintained root CLI. This deliberately avoids sourcing .Rprofile
# so CLI behavior does not depend on interactive/local renv activation.

set -euo pipefail

project_dir="/app"
renv_package_dir="$(find "${project_dir}/renv/library" -maxdepth 6 -type d -name renv -print -quit 2>/dev/null || true)"

if [[ -z "${renv_package_dir}" ]]; then
  printf '%s\n' "Project renv library is missing from the container image." >&2
  exit 1
fi

project_library="$(dirname "${renv_package_dir}")"
export RENV_CONFIG_AUTOLOADER_ENABLED=FALSE
export R_PROFILE_USER=/dev/null
export R_LIBS_USER="${project_library}${R_LIBS_USER:+:${R_LIBS_USER}}"

exec Rscript --vanilla "${project_dir}/entrypoint.R" "$@"
