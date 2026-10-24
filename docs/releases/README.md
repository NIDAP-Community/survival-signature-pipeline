# Release Documentation

This directory contains the Markdown release notes for formally versioned
survival-signature-pipeline releases. Each GitHub Release should reuse the
matching document as its release body and identify:

- the Git tag and target commit;
- the immutable OCI image digest and any discovery tag;
- the validation fixtures and commands used; and
- user-visible interface, output, or compatibility changes.

The currently published GHCR image is a validated source-revision publication,
not yet a semantic-version GitHub Release. Its immutable deployment reference
is recorded in the top-level README and `container/README.md`.
