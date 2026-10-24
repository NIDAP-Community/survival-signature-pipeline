# Container Build and Runtime

This directory packages the maintained root CLI as a Container-as-a-Function
(CAF): a noninteractive process with explicit command-line inputs and an
explicit result directory. The container is a runtime artifact, not a
development shell.

## Canonical entrypoint

The image entrypoint dispatches to the repository-root CLI:

```text
entrypoint.R <workflow> [workflow arguments]
```

Use `entrypoint.R list` to see maintained workflows. Core and driver code
remain in the repository; this directory only defines build and runtime
behavior.

## Distribution model

The canonical artifact is an OCI image published to GitHub Container Registry:

```text
ghcr.io/nidap-community/survival-signature-pipeline@sha256:a25e77f298eb602bf0a5e088a295cce783dd987baf23b7dada4d1cafea9ed35f
```

Deploy callers by digest, not by mutable tag. The full source-revision tag is
available for discovery, but a job record should retain the digest together
with the resulting `RUN_METADATA.md`.

The image targets `linux/amd64`, the same target architecture used by
multi-gene-correlations. Apple Silicon hosts can run it under emulation for
development; that is not a distinct release target.

HPC users can convert the OCI image for Apptainer/Singularity where needed:

```bash
apptainer pull survival-signature-pipeline.sif \
  docker://ghcr.io/nidap-community/survival-signature-pipeline@sha256:a25e77f298eb602bf0a5e088a295cce783dd987baf23b7dada4d1cafea9ed35f
```

No container binaries or `.sif` files belong in Git.

## Apptainer/Singularity runtime helper

`run_apptainer.sh` gives a converted image the same mount contract as
`run_docker.sh`. It chooses `apptainer` when present and otherwise falls back
to `singularity`:

```bash
container/run_apptainer.sh survival-signature-pipeline.sif \
  /path/to/input_directory /path/to/output_directory \
  single-cohort-survival \
  --expr-file /data/expression.tsv \
  --meta-file /data/metadata.csv \
  --output-dir /output
```

The wrapper uses only runtime bind mounts; it does not build images or require
privileged container access.

## Build and validate locally

Build a local `linux/amd64` image and run the canonical smoke test:

```bash
container/build_docker.sh survival-signature-pipeline:local
container/smoke_test.sh survival-signature-pipeline:local
```

On HTTPS-inspecting networks, `DOCKER_CA_CERT` may point to a trusted PEM
bundle. The helper forwards it only as a BuildKit secret for dependency
retrieval; it is not included in an image layer. See
[the container contract](../docs/container_contract.md) for details.

## Runtime contract

Mount caller-provided inputs read-only at `/data` and a new writable result
directory at `/output`. `container/run_docker.sh` applies that layout and
passes workflow arguments through unchanged:

```bash
container/run_docker.sh <image> <data_dir> <output_dir> \
  single-cohort-survival --expr-file /data/expression.tsv \
  --meta-file /data/metadata.csv --output-dir /output
```

The image runs as a non-root user by default and contains no study data,
database credentials, or cloud credentials. The complete runtime and
OncoChat/AWS integration boundary are documented in
[docs/container_contract.md](../docs/container_contract.md) and
[docs/aws_adapter_contract.md](../docs/aws_adapter_contract.md).
