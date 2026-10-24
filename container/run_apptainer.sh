#!/usr/bin/env bash
# Run a converted OCI image with the same explicit /data -> /output contract
# as run_docker.sh. Either Apptainer or legacy Singularity may provide the
# runtime; this wrapper selects the available command without requiring a
# privileged build or any runtime network access.

set -euo pipefail

if [[ "$#" -lt 4 ]]; then
  cat >&2 <<'USAGE'
Usage:
  container/run_apptainer.sh <image.sif> <data_dir> <output_dir> <workflow> [workflow arguments]

The complete input directory is bound read-only at /data. Pass /data/<file>
paths to the workflow and set --output-dir /output (or a subdirectory of it).
USAGE
  exit 2
fi

image="$1"
data_dir="$2"
output_dir="$3"
shift 3

if [[ ! -f "${image}" ]]; then
  printf 'Image file does not exist: %s\n' "${image}" >&2
  exit 2
fi

if [[ ! -d "${data_dir}" ]]; then
  printf 'Input directory does not exist: %s\n' "${data_dir}" >&2
  exit 2
fi

if command -v apptainer >/dev/null 2>&1; then
  runtime="apptainer"
elif command -v singularity >/dev/null 2>&1; then
  runtime="singularity"
else
  printf '%s\n' 'Neither apptainer nor singularity is available on PATH.' >&2
  exit 127
fi

# -p is portable across macOS and Linux shells. Resolve paths before mounting
# so that the runtime receives stable absolute host paths.
mkdir -p "${output_dir}"
data_dir="$(CDPATH= cd -- "${data_dir}" && pwd -P)"
output_dir="$(CDPATH= cd -- "${output_dir}" && pwd -P)"

exec "${runtime}" run \
  --bind "${data_dir}:/data:ro" \
  --bind "${output_dir}:/output" \
  "${image}" "$@"
