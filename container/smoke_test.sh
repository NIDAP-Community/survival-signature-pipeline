#!/usr/bin/env bash
# Run the maintained expression, binary-biomarker, and continuous-biomarker
# fixtures, including immune-score and cell-fraction contexts, inside an
# already built image and assert the minimal portable artifact contract. The
# temporary output directory is retained after completion so failures and
# generated files are inspectable.

set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
project_dir="$(CDPATH= cd -- "${script_dir}/.." && pwd)"
image_tag="${1:-survival-signature-pipeline:local}"
smoke_root="$(mktemp -d "${TMPDIR:-/tmp}/survival-container-smoke.XXXXXX")"
expression_output_dir="${smoke_root}/expression_output"
binary_output_dir="${smoke_root}/binary_output"
continuous_output_dir="${smoke_root}/continuous_output"
estimate_output_dir="${smoke_root}/estimate_output"
fraction_output_dir="${smoke_root}/fraction_output"

"${script_dir}/run_docker.sh" \
  "${image_tag}" \
  "${project_dir}/test_data" \
  "${expression_output_dir}" \
  single-cohort-survival \
  --expr-file /data/tcga_expression.tsv \
  --meta-file /data/tcga_metadata_no_subtype.csv \
  --features GENE_D \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir /output

"${script_dir}/run_docker.sh" \
  "${image_tag}" \
  "${project_dir}/test_data" \
  "${estimate_output_dir}" \
  single-cohort-continuous-biomarker-survival \
  --biomarker-file /data/estimate_score_matrix.tsv \
  --meta-file /data/tcga_metadata_no_subtype.csv \
  --features ImmuneScore \
  --biomarker-context estimate_scores \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir /output

"${script_dir}/run_docker.sh" \
  "${image_tag}" \
  "${project_dir}/test_data" \
  "${fraction_output_dir}" \
  single-cohort-continuous-biomarker-survival \
  --biomarker-file /data/cell_fraction_matrix.tsv \
  --meta-file /data/tcga_metadata_no_subtype.csv \
  --features CD8_T_cells \
  --biomarker-context cell_fractions \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir /output

"${script_dir}/run_docker.sh" \
  "${image_tag}" \
  "${project_dir}/test_data" \
  "${continuous_output_dir}" \
  single-cohort-continuous-biomarker-survival \
  --biomarker-file /data/continuous_biomarker_matrix.tsv \
  --meta-file /data/tcga_metadata_no_subtype.csv \
  --features MYC_log2_copy_number \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir /output

"${script_dir}/run_docker.sh" \
  "${image_tag}" \
  "${project_dir}/test_data" \
  "${binary_output_dir}" \
  single-cohort-binary-biomarker-survival \
  --biomarker-file /data/binary_biomarker_matrix.tsv \
  --meta-file /data/tcga_metadata_no_subtype.csv \
  --features TP53_mutated \
  --endpoints OS \
  --adjust-covariates none \
  --output-dir /output

required_files=(
  "${expression_output_dir}/RUN_METADATA.md"
  "${expression_output_dir}/logs/sessionInfo.txt"
  "${expression_output_dir}/run_manifest.csv"
  "${expression_output_dir}/OS/results_univariate.csv"
  "${expression_output_dir}/OS/grouped_gene_cox_results.csv"
  "${expression_output_dir}/OS/OS_results.rds"
  "${expression_output_dir}/OS/KM_GENE_D.png"
  "${binary_output_dir}/RUN_METADATA.md"
  "${binary_output_dir}/logs/sessionInfo.txt"
  "${binary_output_dir}/binary_biomarker_definition.csv"
  "${binary_output_dir}/run_manifest.csv"
  "${binary_output_dir}/OS/results_univariate.csv"
  "${binary_output_dir}/OS/grouped_biomarker_cox_results.csv"
  "${binary_output_dir}/OS/OS_results.rds"
  "${binary_output_dir}/OS/KM_TP53_mutated.png"
  "${continuous_output_dir}/RUN_METADATA.md"
  "${continuous_output_dir}/logs/sessionInfo.txt"
  "${continuous_output_dir}/continuous_biomarker_definition.csv"
  "${continuous_output_dir}/run_manifest.csv"
  "${continuous_output_dir}/OS/results_univariate.csv"
  "${continuous_output_dir}/OS/grouped_biomarker_cox_results.csv"
  "${continuous_output_dir}/OS/OS_results.rds"
  "${continuous_output_dir}/OS/KM_MYC_log2_copy_number.png"
  "${estimate_output_dir}/RUN_METADATA.md"
  "${estimate_output_dir}/continuous_biomarker_definition.csv"
  "${estimate_output_dir}/OS/KM_ImmuneScore.png"
  "${fraction_output_dir}/RUN_METADATA.md"
  "${fraction_output_dir}/continuous_biomarker_definition.csv"
  "${fraction_output_dir}/OS/KM_CD8_T_cells.png"
)

for required_file in "${required_files[@]}"; do
  if [[ ! -s "${required_file}" ]]; then
    printf 'Container smoke test is missing required artifact: %s\n' "${required_file}" >&2
    exit 1
  fi
done

printf 'Container smoke test passed. Artifacts retained at: %s\n' "${smoke_root}"
