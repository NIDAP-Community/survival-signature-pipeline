# AWS Adapter Contract

This document defines the handoff between the portable
survival-signature-pipeline image and an AWS-hosted application such as
OncoChat. It is an integration contract, not AWS infrastructure code. The
survival repository remains responsible for scientific analysis behavior and
the command-line interface; the calling application owns request handling,
authorization, data staging, job orchestration, and result presentation.

## Versioned Image and Job Identity

The adapter must launch an immutable OCI image digest, not a mutable tag. For
every submitted job, retain:

- image digest and source commit;
- submitted request after validation;
- input object identifiers and checksums;
- task/job identifier, timestamps, resource limits, and terminal status; and
- the produced RUN_METADATA.md.

The adapter must not infer successful completion from a container start. A
successful job is one that exits successfully and contains the required
run-scoped artifacts.

## Input Staging Contract

The adapter stages only approved files into a job-specific, read-only
/data directory. It allocates a fresh writable /output directory for each job
and passes /output (or an immediate child) as --output-dir.

At a minimum, a single-cohort request supplies:

| Logical input | Required content |
| --- | --- |
| Expression table | Delimited table with a unique gene identifier column and sample columns, for expression/signature requests. |
| Binary biomarker table | Delimited table with a unique feature column and subject columns containing only 0, 1, or missing, for binary-biomarker requests. |
| Metadata table | Delimited table with a unique sample identifier column and paired endpoint event/time columns. |
| Optional analysis plan | CSV plan only when running single-cohort-analysis-plan. |

Before execution, the adapter must verify that each requested path resolves
inside /data; it must reject absolute host paths, .. traversal, duplicate
logical inputs, unreadable objects, and unapproved input locations.

## Initial Application Request Surface

For the first application integration, expose `single-cohort-survival` and,
when mutation calls are available, `single-cohort-binary-biomarker-survival`.
Both have bounded, well-documented input/output contracts. The adapter should
translate validated form values into explicit CLI flags rather than accepting
arbitrary shell arguments.

The application request should include:

| Field | Adapter action |
| --- | --- |
| Expression and metadata input selections | Stage the selected approved files and pass their /data paths. |
| Sample and gene identifier columns | Verify the requested columns exist before launch. |
| Genes or signature members | Parse a bounded list of nonempty identifiers; use a signature only with one or more signature genes. |
| Optional display label | Pass the label after length and character validation. |
| Endpoints | Allow only configured endpoint labels and their configured event/time mappings. The driver defaults to OS,PFI. |
| Covariates and cohort filter | Allow only metadata columns verified for the staged table. |
| Gene resolution | Offer none and l2p; when l2p is selected, surface the generated mapping report. |
| Grouping and analysis settings | Constrain values to the driver’s documented enums and safe numeric ranges. |
| Binary biomarker selections | For the binary workflow, stage an approved 0/1 matrix, verify feature and subject identifiers, and pass only configured labels for 0 and 1. Do not ask the container to query PostgreSQL or derive mutation calls. |

The adapter should set a unique output directory under /output itself. It must
not let callers choose an arbitrary output path.

## Execution and Artifacts

Run the image with /data read-only and /output writable, using the container's
non-root default user. The adapter owns task CPU, memory, timeout, retry,
logging, cancellation, and retention policies.

For a completed single-cohort-survival run, collect at least:

    <job_id>/RUN_METADATA.md
    <job_id>/logs/sessionInfo.txt
    <job_id>/run_manifest.csv
    <job_id>/tables/gene_identifier_resolution.csv
    <job_id>/OS/ and/or PFI/

Each endpoint directory may include Kaplan-Meier figures, model tables,
rendered summaries, and (when enabled) an RDS result object. Preserve the
complete result directory rather than selecting a result based only on a file
name.

The overall forest plot is a downstream single-cohort-reporting-bundle
operation. If the application offers it, retain saved result RDS objects for
the analysis run and submit the reporting bundle as a linked follow-on job
using the first job's result directory as its approved input. Do not make the
reporting bundle silently rerun the survival analysis.

For a completed binary-biomarker run, collect its run-level provenance and
each endpoint directory, including `grouped_biomarker_cox_results.csv`, the
endpoint summary, and any generated KM figures. Do not present an input value
as a mutation unless the staged upstream data contract identifies it as one.

## Required Application Guards

The adapter must:

1. authorize access to the requested datasets before staging;
2. validate schemas and option values before starting a container;
3. associate logs and output objects with the requesting user/job;
4. present failures without exposing host paths, credentials, or unrelated
   tenant data;
5. retain the image digest, provenance metadata, and input checksums with
   results; and
6. require an explicit user acknowledgement that outputs are exploratory
   analyses, not clinical decision support.

## Acceptance Tests for the Integration Owner

Before enabling the application, verify against the target AWS runtime:

1. an approved small fixture completes using a pinned image digest;
2. both default endpoints (OS and PFI) are collected when requested;
3. an l2p request produces a nonempty gene_identifier_resolution.csv and
   resolves a known historic symbol;
4. malformed paths, unsupported options, missing columns, and unauthorized
   datasets are rejected before job launch;
5. a reporting-bundle follow-on job can produce the overall forest plot; and
6. a real approved study fixture reproduces its pre-specified reference
   results within the agreed tolerance.

Passing a local container test or an image build alone does not satisfy these
application acceptance tests.
