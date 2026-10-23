# Container Contract

The survival-signature-pipeline OCI image packages the maintained root CLI for portable batch execution. It is intentionally independent of AWS infrastructure, database access, credentials, and OncoChat-specific request formats.

## Runtime Interface

The image entrypoint is equivalent to:

```bash
Rscript --vanilla entrypoint.R <workflow> [workflow arguments]
```

The container expects two mount points:

| Container path | Access | Purpose |
| --- | --- | --- |
| `/data` | read-only | Input expression, metadata, and plan files supplied by the caller. |
| `/output` | writable | Complete run-scoped output directory supplied by the caller. |

Pass `/data/...` input paths and `/output` (or a child directory) to `--output-dir`. The image runs as a non-root `survival` user. `container/run_docker.sh` overrides that user with the invoking host UID solely to make bind-mounted output files host-writable.

The image has no embedded clinical data and requires no network access at execution time. Its build stage restores exactly the R dependency versions recorded in `renv.lock`.

## Local Build and Smoke Test

Build natively for the local Docker host:

```bash
container/build_docker.sh survival-signature-pipeline:local
```

For a target platform, make the architecture explicit:

```bash
TARGET_PLATFORM=linux/arm64 container/build_docker.sh survival-signature-pipeline:arm64
TARGET_PLATFORM=linux/amd64 container/build_docker.sh survival-signature-pipeline:amd64
```

`container/build_docker.sh` is intentionally a local-build helper and loads a
single-platform image into the local Docker image store. A release workflow
should build and push a separately validated image by immutable digest; it
should not rely on this local tag as release evidence.

Run the image's canonical smoke test after each target-platform build:

```bash
container/smoke_test.sh survival-signature-pipeline:local
```

The smoke test retains its output directory and verifies that the container creates run metadata, session information, single-cohort result tables/RDS, and a Kaplan-Meier PNG.

## Docker Invocation

```bash
container/run_docker.sh \
  survival-signature-pipeline:local \
  /path/to/input_directory \
  /path/to/output_directory \
  single-cohort-analysis-plan \
  --expr-file /data/expression.tsv \
  --meta-file /data/metadata.csv \
  --plan-file /data/analysis_plan.csv \
  --output-dir /output
```

## Integration Boundary

An AWS integration may stage approved inputs to `/data` and collect `/output` artifacts, but it must remain outside this repository. A caller integrating with a database or an agent runtime is responsible for request validation, data materialization, storage access, IAM, job lifecycle, artifact delivery, and its own infrastructure as code.

For deployment reproducibility, record the immutable OCI image digest together with the existing `RUN_METADATA.md` file. Do not treat a locally successful container run as evidence of a deployed image, AWS task, or release.

## Continuous Validation

`.github/workflows/container.yml` builds and smoke-tests the native Linux image and separately builds the ARM64 image on pull requests and pushes to `main`. It deliberately has no registry credentials and does not publish an image. Image publication, registry retention, signing, vulnerability policy, and AWS deployment remain release/integration responsibilities.
