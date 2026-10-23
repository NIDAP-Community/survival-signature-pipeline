#!/usr/bin/env bash
# Run the canonical single-cohort fixture inside an already built image and
# assert the minimal portable artifact contract. The temporary output directory
# is retained after completion so failures and generated files are inspectable.

set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(CDPATH= cd -- "${script_dir}/.." && pwd)"
image_tag="${1:-survival-signature-pipeline:local}"
smoke_root="$(mktemp -d "${TMPDIR:-/tmp}/survival-container-smoke.XXXXXX")"
output_dir="${smoke_root}/output"

"${script_dir}/run_docker.sh" \
  "${image_tag}" \
  "${project_dir}/test_data" \
  "${output_dir}" \
  single-cohort-survival \
  --expr-file /data/tcga_expression.tsv \
  --meta-file /data/tcga_metadata_no_subtype.csv \
  --features GENE_D \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir /output

required_files=(
  "${output_dir}/RUN_METADATA.md"
  "${output_dir}/logs/sessionInfo.txt"
  "${output_dir}/run_manifest.csv"
  "${output_dir}/OS/results_univariate.csv"
  "${output_dir}/OS/grouped_gene_cox_results.csv"
  "${output_dir}/OS/OS_results.rds"
  "${output_dir}/OS/KM_GENE_D.png"
)

for required_file in "${required_files[@]}"; do
  if [[ ! -s "${required_file}" ]]; then
    printf 'Container smoke test is missing required artifact: %s\n' "${required_file}" >&2
    exit 1
  fi
done

printf 'Container smoke test passed. Artifacts retained at: %s\n' "${output_dir}"
