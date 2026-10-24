#!/usr/bin/env bash
# Build a locally tagged OCI image. Set TARGET_PLATFORM explicitly for a
# deployable architecture, for example linux/arm64 for the current AgentCore
# path or linux/amd64 for an x86_64 Fargate task.

set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(CDPATH= cd -- "${script_dir}/.." && pwd)"
image_tag="${1:-survival-signature-pipeline:local}"
vcs_ref="$(git -C "${project_dir}" rev-parse HEAD)"
build_date="$(date -u +%Y-%m-%dT%H:%M:%SZ)"

build_args=(
  --file "${project_dir}/container/Dockerfile"
  --progress "${DOCKER_BUILD_PROGRESS:-auto}"
  --load
  --tag "${image_tag}"
  --build-arg "VCS_REF=${vcs_ref}"
  --build-arg "BUILD_DATE=${build_date}"
)

if [[ -n "${TARGET_PLATFORM:-}" ]]; then
  build_args+=(--platform "${TARGET_PLATFORM}")
fi

# On a network that inspects HTTPS, provide its PEM trust bundle without
# embedding it in the build context or final image. The Dockerfile consumes it
# only during dependency retrieval. Leave this unset for ordinary builds.
if [[ -n "${DOCKER_CA_CERT:-}" ]]; then
  if [[ ! -r "${DOCKER_CA_CERT}" ]]; then
    printf '%s\n' "DOCKER_CA_CERT is not a readable PEM file: ${DOCKER_CA_CERT}" >&2
    exit 2
  fi
  build_args+=(--secret "id=build_ca,src=${DOCKER_CA_CERT}")
fi

# Buildx makes --load explicit: the resulting single-platform image is usable
# immediately by container/run_docker.sh and container/smoke_test.sh.
docker buildx build "${build_args[@]}" "${project_dir}"
