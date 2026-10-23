#!/usr/bin/env bash
# Run the container with an explicit, read-only data mount and writable output
# mount. Keeping the mount contract here avoids host paths leaking into saved
# run metadata and matches the future AWS staging boundary.

set -euo pipefail

if [[ "$#" -lt 4 ]]; then
  cat >&2 <<'USAGE'
Usage:
  container/run_docker.sh <image> <data_dir> <output_dir> <workflow> [workflow arguments]

The complete input directory is mounted read-only at /data. Pass /data/<file>
paths to the workflow and set --output-dir /output (or a subdirectory of it).
USAGE
  exit 2
fi

image="$1"
data_dir="$2"
output_dir="$3"
shift 3

if [[ ! -d "${data_dir}" ]]; then
  printf 'Input directory does not exist: %s\n' "${data_dir}" >&2
  exit 2
fi

mkdir --parents "${output_dir}"
data_dir="$(CDPATH= cd -- "${data_dir}" && pwd -P)"
output_dir="$(CDPATH= cd -- "${output_dir}" && pwd -P)"

exec docker run --rm \
  --user "$(id -u):$(id -g)" \
  --mount "type=bind,src=${data_dir},dst=/data,readonly" \
  --mount "type=bind,src=${output_dir},dst=/output" \
  "${image}" "$@"
