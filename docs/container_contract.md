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

If a local institutional network transparently inspects HTTPS, retain
certificate verification by providing its PEM trust bundle only for dependency
retrieval. The helper passes this file as a BuildKit secret; it is not copied
into the image or repository:

```bash
DOCKER_CA_CERT=/path/to/trusted-network-root.pem \
  container/build_docker.sh survival-signature-pipeline:local
```

Do not disable TLS verification or add an institution-specific certificate to
the Dockerfile. This setting is unnecessary for standard internet-connected
or AWS builds.

The container follows the same single-architecture convention as
multi-gene-correlation: `linux/amd64`.

```bash
container/build_docker.sh survival-signature-pipeline:local
```

`container/build_docker.sh` is intentionally a local-build helper and loads a
single-platform image into the local Docker image store. A release workflow
should build and push a separately validated image by immutable digest; it
should not rely on this local tag as release evidence.

Run the image's canonical smoke test after each build:

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

For the required adapter-side request validation, input staging, artifact
collection, and acceptance checks, see the
[AWS adapter contract](aws_adapter_contract.md).

## Continuous Validation

`.github/workflows/container.yml` builds and smoke-tests the linux/amd64 image on pull requests and pushes to `main`. It deliberately has no registry credentials and does not publish an image. Image publication, registry retention, signing, vulnerability policy, and AWS deployment remain release/integration responsibilities.
