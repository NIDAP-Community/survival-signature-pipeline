#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("src/core/GeneExpressionBasedSurvivalAnalys.r")
source("src/utilities/cli_helpers.R")
source("src/utilities/run_provenance_guardrails.R")

cli_defaults <- list(
  output_dir = default_run_output_dir("single_cohort_survival"),
  expr_file = file.path("test_data", "tcga_expression.tsv"),
  meta_file = file.path("test_data", "tcga_metadata.csv"),
  sample_column = "sampnames",
  gene_column = "Gene",
  subtype_column = "Subtype",
  stage_column = "ajcc_pathologic_tumor_stage",
  age_column = "Diagnosis_Age",
  signature_genes = "SAT1,SRM,ODC1,AMD1",
  os_event_column = "OS",
  os_time_column = "OS_time",
  pfi_event_column = "PFI",
  pfi_time_column = "PFI_time",
  subtype_filter = "ALL",
  primary_filter = "tcga_primary",
  primary_pattern = "_01",
  adjust_covariates = "subtype_stage_age",
  split_method = "tertiles",
  include_groups = "low,high",
  features = "Signature,AMD1,SRM,SAT1,ODC1",
  endpoints = "OS,PFI",
  cut_months = "60",
  run_km_individual = "false",
  save_plots = "false",
  export_csv = "true",
  run_lasso = "false",
  verbose = "true"
)

print_usage <- function() {
  print_standard_cli_help(
    title = "Single-Cohort Survival Driver",
    workflow_name = "single-cohort-survival",
    script_path = file.path("drivers", "run_single_cohort_survival.R"),
    defaults = cli_defaults,
    summary_lines = c(
      "subtype_column is optional unless subtype filtering or subtype-adjusted modeling is requested.",
      "primary_filter=tcga_primary uses primary_pattern to subset sample IDs.",
      "adjust_covariates supports: none, subtype, stage, age, subtype_stage, subtype_age, stage_age, subtype_stage_age."
    )
  )
}

is_truthy <- function(value) {
  if (is.null(value) || is.na(value)) return(FALSE)
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

resolve_endpoint_columns <- function(endpoint_name, options) {
  resolve_standard_endpoint_columns(endpoint_name, options, "single-cohort-survival")
}

determine_covariates <- function(mode, subtype_column, stage_column, age_column) {
  mode <- tolower(mode)
  if (mode == "none") return(character(0))
  covariates <- character(0)
  if (grepl("subtype", mode)) covariates <- c(covariates, subtype_column)
  if (grepl("stage", mode)) covariates <- c(covariates, stage_column)
  if (grepl("age", mode)) covariates <- c(covariates, age_column)
  unique(covariates)
}

write_simple_markdown <- function(lines, file_path) {
  writeLines(lines, con = file_path)
}

build_endpoint_summary <- function(endpoint_name, endpoint_dir, result, feature_list, covariate_mode) {
  lines <- c(
    paste0("# ", endpoint_name, " Single-Cohort Summary"),
    "",
    paste0("- Endpoint: ", endpoint_name),
    paste0("- Covariate mode: ", covariate_mode),
    paste0("- Features: ", paste(feature_list, collapse = ", ")),
    ""
  )

  if (!is.null(result$results_univariate) && nrow(result$results_univariate)) {
    lines <- c(lines, "## Univariate Cox", "", utils::capture.output(print(head(result$results_univariate, 10))), "")
  }
  if (!is.null(result$results_adjusted) && nrow(result$results_adjusted)) {
    lines <- c(lines, "## Adjusted Cox", "", utils::capture.output(print(head(result$results_adjusted, 10))), "")
  }
  if (!is.null(result$results_signature)) {
    lines <- c(lines, "## Signature Cox", "")
    if (!is.null(result$results_signature$unadjusted)) {
      lines <- c(lines, "### Unadjusted", "", utils::capture.output(print(result$results_signature$unadjusted)), "")
    }
    if (!is.null(result$results_signature$adjusted)) {
      lines <- c(lines, "### Adjusted", "", utils::capture.output(print(result$results_signature$adjusted)), "")
    }
  }

  write_simple_markdown(lines, file.path(endpoint_dir, paste0(endpoint_name, "_summary.md")))
}

opt <- parse_standard_cli_args(
  defaults = cli_defaults,
  workflow_name = "single-cohort-survival",
  help_printer = print_usage
)
opt$cut_months <- cli_parse_numeric_flag(opt$cut_months, "cut_months")
opt$adjust_covariates <- cli_parse_enum_flag(
  value = opt$adjust_covariates,
  allowed_values = c("none", "subtype", "stage", "age", "subtype_stage", "subtype_age", "stage_age", "subtype_stage_age"),
  flag_name = "adjust_covariates",
  normalize = "lower"
)
feature_list <- cli_parse_csv_list(opt$features)
endpoint_list <- cli_parse_csv_list(opt$endpoints, uppercase = TRUE)
signature_genes <- cli_parse_csv_list(opt$signature_genes)
include_groups <- cli_parse_csv_list(opt$include_groups, lowercase = TRUE)
has_signature <- any(tolower(feature_list) == "signature")
individual_features <- feature_list[tolower(feature_list) != "signature"]
analysis_genes <- unique(c(individual_features, if (has_signature) signature_genes else character(0)))

assert_output_guardrails(opt$output_dir)

if (!length(analysis_genes)) {
  stop("No analyzable genes resolved from --features and --signature-genes.")
}

cli_assert_existing_input_files(
  input_paths = c(
    expr_file = opt$expr_file,
    meta_file = opt$meta_file
  ),
  workflow_name = "single-cohort-survival"
)

expression_data <- read_expression_table(opt$expr_file, opt$gene_column)
metadata <- read_metadata_table(opt$meta_file)

if (!opt$sample_column %in% names(metadata)) stop("Metadata missing sample column: ", opt$sample_column)

subtype_clean <- tolower(trimws(opt$subtype_filter))
subtype_filter_requested <- !(subtype_clean %in% c("", "all", "overall", "any"))
subtype_adjustment_requested <- grepl("subtype", tolower(opt$adjust_covariates))
has_subtype_column <- opt$subtype_column %in% names(metadata)

if (!has_subtype_column && subtype_filter_requested) {
  stop(
    "Metadata missing subtype column '", opt$subtype_column,
    "'. subtype_column is optional unless subtype filtering or subtype-adjusted modeling is requested. ",
    "Either provide --subtype-column for metadata that includes subtypes, or remove --subtype-filter."
  )
}

if (!has_subtype_column && subtype_adjustment_requested) {
  stop(
    "Metadata missing subtype column '", opt$subtype_column,
    "'. subtype_column is optional unless subtype filtering or subtype-adjusted modeling is requested. ",
    "Either provide --subtype-column for metadata that includes subtypes, or choose an --adjust-covariates mode without subtype."
  )
}

metadata_filtered <- metadata
if (tolower(opt$primary_filter) == "tcga_primary") {
  metadata_filtered <- metadata_filtered %>% filter(grepl(opt$primary_pattern, .data[[opt$sample_column]]))
}

if (has_subtype_column && subtype_filter_requested) {
  selected_subtypes <- cli_parse_csv_list(opt$subtype_filter)
  metadata_filtered <- metadata_filtered %>% filter(.data[[opt$subtype_column]] %in% selected_subtypes)
}

if (!nrow(metadata_filtered)) {
  stop("No metadata rows remain after filtering. Check --primary-filter and --subtype-filter.")
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

covariates <- determine_covariates(opt$adjust_covariates, opt$subtype_column, opt$stage_column, opt$age_column)
covariates <- covariates[covariates %in% names(metadata_filtered)]

steps_to_run <- character(0)
if (length(individual_features)) {
  steps_to_run <- c(steps_to_run, "individual")
  if (length(covariates)) steps_to_run <- c(steps_to_run, "individual+covariates")
}
if (has_signature) {
  steps_to_run <- c(steps_to_run, "signature")
  if (length(covariates)) steps_to_run <- c(steps_to_run, "signature+covariates")
}
if (is_truthy(opt$run_km_individual) && length(individual_features)) {
  steps_to_run <- c(steps_to_run, "km_individual")
}
if (is_truthy(opt$run_lasso) && length(analysis_genes)) {
  steps_to_run <- c(steps_to_run, "lasso")
}

results_manifest <- list()

for (endpoint_name in endpoint_list) {
  endpoint_columns <- resolve_endpoint_columns(endpoint_name, opt)
  endpoint_dir <- file.path(opt$output_dir, endpoint_name)
  dir.create(endpoint_dir, recursive = TRUE, showWarnings = FALSE)

  result <- analyze_survival_expression(
    expression_data = expression_data,
    metadata = metadata_filtered,
    sample_col = opt$sample_column,
    gene_col = opt$gene_column,
    genes = analysis_genes,
    event_col = endpoint_columns$event,
    time_col = endpoint_columns$time,
    clinical_covariates = covariates,
    split_method = opt$split_method,
    include_groups = include_groups,
    steps_to_run = steps_to_run,
    run_lasso = is_truthy(opt$run_lasso),
    run_ph_checks = TRUE,
    run_time_dependent = TRUE,
    cut_months = opt$cut_months,
    save_plots = is_truthy(opt$save_plots),
    export_csv = is_truthy(opt$export_csv),
    output_dir = endpoint_dir,
    plot_params = list(plot_title = paste0(endpoint_name, " Single-Cohort Survival")),
    verbose = is_truthy(opt$verbose)
  )

  saveRDS(result, file = file.path(endpoint_dir, paste0(endpoint_name, "_results.rds")))
  build_endpoint_summary(endpoint_name, endpoint_dir, result, feature_list, opt$adjust_covariates)

  results_manifest[[length(results_manifest) + 1]] <- data.frame(
    endpoint = endpoint_name,
    endpoint_dir = endpoint_dir,
    metadata_rows = nrow(metadata_filtered),
    stringsAsFactors = FALSE
  )
}

manifest_df <- bind_rows(results_manifest)
readr::write_csv(manifest_df, file.path(opt$output_dir, "run_manifest.csv"))

write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Run single-cohort survival analysis through the reusable core analysis engine.",
  driver = "drivers/run_single_cohort_survival.R",
  inputs = list(
    expr_file = opt$expr_file,
    meta_file = opt$meta_file
  ),
  parameters = list(
    sample_column = opt$sample_column,
    gene_column = opt$gene_column,
    subtype_column = opt$subtype_column,
    stage_column = opt$stage_column,
    age_column = opt$age_column,
    signature_genes = paste(signature_genes, collapse = ","),
    features = opt$features,
    endpoints = opt$endpoints,
    analysis_genes = paste(analysis_genes, collapse = ","),
    subtype_filter = opt$subtype_filter,
    primary_filter = opt$primary_filter,
    primary_pattern = opt$primary_pattern,
    adjust_covariates = opt$adjust_covariates,
    split_method = opt$split_method,
    include_groups = paste(include_groups, collapse = ","),
    steps_to_run = paste(steps_to_run, collapse = ","),
    cut_months = opt$cut_months,
    run_km_individual = opt$run_km_individual,
    save_plots = opt$save_plots,
    export_csv = opt$export_csv,
    run_lasso = opt$run_lasso
  ),
  scope = "exploratory"
)

message("Single-cohort survival analysis written to: ", opt$output_dir)