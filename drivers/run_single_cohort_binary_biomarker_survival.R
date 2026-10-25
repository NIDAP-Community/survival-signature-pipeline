#!/usr/bin/env Rscript

# Maintained workflow for a prepared, patient-level binary biomarker matrix.
# The calling application is responsible for database extraction, authorization,
# and harmonizing one value per analysis subject before invoking this driver.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("src/core/GeneExpressionBasedSurvivalAnalys.r")
source("src/utilities/cli_helpers.R")
source("src/utilities/reporting_helpers.R")
source("src/utilities/run_provenance_guardrails.R")

cli_defaults <- list(
  output_dir = default_run_output_dir("single_cohort_binary_biomarker_survival"),
  biomarker_file = file.path("test_data", "binary_biomarker_matrix.tsv"),
  meta_file = file.path("test_data", "tcga_metadata_no_subtype.csv"),
  sample_column = "sampnames",
  feature_column = "Feature",
  features = "TP53_mutated",
  positive_label = "Mutated",
  negative_label = "Wild type",
  os_event_column = "OS",
  os_time_column = "OS_time",
  pfi_event_column = "PFI",
  pfi_time_column = "PFI_time",
  endpoints = "OS,PFI",
  adjust_covariates = "none",
  group_column = "",
  group_filter = "ALL",
  run_ph_checks = "true",
  run_time_dependent = "true",
  save_results_rds = "true",
  save_plots = "true",
  export_csv = "true",
  verbose = "true"
)

print_usage <- function() {
  print_standard_cli_help(
    title = "Single-Cohort Binary Biomarker Survival Driver",
    workflow_name = "single-cohort-binary-biomarker-survival",
    script_path = file.path("drivers", "run_single_cohort_binary_biomarker_survival.R"),
    defaults = cli_defaults,
    summary_lines = c(
      "Use this workflow for a prepared feature-by-subject 0/1 matrix, such as mutation status.",
      "0 is the reference group and 1 is the comparison group; missing values mean the feature was not evaluable for that subject.",
      "The driver does not query a database and does not infer mutation calls from expression data.",
      "adjust_covariates accepts explicit metadata column names or none.",
      "group_column/group_filter optionally restrict the prepared cohort before modeling."
    )
  )
}

is_truthy <- function(value) {
  !is.null(value) && !is.na(value) && tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

is_filter_requested <- function(value) {
  !(tolower(trimws(as.character(value))) %in% c("", "all", "overall", "any"))
}

resolve_endpoint_columns <- function(endpoint_name, options) {
  resolve_standard_endpoint_columns(endpoint_name, options, "single-cohort-binary-biomarker-survival")
}

validate_binary_matrix <- function(biomarker_data, feature_column, selected_features) {
  if (!feature_column %in% names(biomarker_data)) {
    stop("Biomarker table missing feature column: ", feature_column)
  }
  if (anyDuplicated(biomarker_data[[feature_column]])) {
    stop("Biomarker table contains duplicate feature identifiers in column: ", feature_column)
  }
  if (!all(selected_features %in% biomarker_data[[feature_column]])) {
    stop(
      "Requested features are absent from the biomarker table: ",
      paste(setdiff(selected_features, biomarker_data[[feature_column]]), collapse = ", ")
    )
  }

  sample_columns <- setdiff(names(biomarker_data), feature_column)
  if (!length(sample_columns)) stop("Biomarker table must include one or more subject columns.")

  selected <- biomarker_data[biomarker_data[[feature_column]] %in% selected_features, sample_columns, drop = FALSE]
  raw_values <- as.matrix(selected)
  numeric_values <- suppressWarnings(as.numeric(raw_values))
  invalid <- !is.na(raw_values) & (!is.finite(numeric_values) | !(numeric_values %in% c(0, 1)))
  if (any(invalid)) {
    stop("Binary biomarker values must be exactly 0, 1, or missing; fractional values are not accepted.")
  }
  invisible(sample_columns)
}

write_binary_endpoint_summary <- function(endpoint_name, endpoint_dir, result, feature_list,
                                          negative_label, positive_label, covariates) {
  grouped <- result$results_grouped_gene
  lines <- c(
    paste0("# ", endpoint_name, " Binary Biomarker Survival Summary"),
    "",
    paste0("- Endpoint: ", endpoint_name),
    paste0("- Binary comparison: ", positive_label, " (1) vs ", negative_label, " (0)"),
    paste0("- Features: ", paste(feature_list, collapse = ", ")),
    paste0("- Requested covariates: ", if (length(covariates)) paste(covariates, collapse = ", ") else "none"),
    "",
    "The Kaplan-Meier figure is an unadjusted display. Grouped Cox rows estimate the same binary comparison; adjusted rows include the retained complete-case covariates.",
    ""
  )
  if (!is.null(grouped) && nrow(grouped)) {
    display <- grouped %>%
      transmute(
        Feature = feature,
        Model = model_type,
        Comparison = paste0(positive_label, " vs ", negative_label),
        `Hazard ratio` = report_fmt_num(HR, 3),
        `95% CI` = report_fmt_ci(CI_lower, CI_upper, 3),
        `p value` = report_fmt_p(p_value),
        N = N,
        Events = events,
        `Retained covariates` = covariates_retained
      )
    lines <- c(lines, "## Grouped Cox Results", "", markdown_table_lines(display), "")
  }
  writeLines(lines, file.path(endpoint_dir, paste0(endpoint_name, "_summary.md")))
}

write_binary_grouped_results <- function(endpoint_name, endpoint_dir, result,
                                         negative_label, positive_label) {
  grouped <- result$results_grouped_gene
  if (is.null(grouped) || !nrow(grouped)) return(invisible(NULL))
  output <- grouped %>%
    mutate(
      endpoint = endpoint_name,
      comparison = paste0(positive_label, " vs ", negative_label),
      reference_group = negative_label,
      comparison_group = positive_label
    ) %>%
    select(endpoint, feature, comparison, reference_group, comparison_group,
           HR, CI_lower, CI_upper, p_value, N, events, model_type,
           covariates_requested, covariates_retained)
  readr::write_csv(output, file.path(endpoint_dir, "grouped_biomarker_cox_results.csv"))
}

opt <- parse_standard_cli_args(
  defaults = cli_defaults,
  workflow_name = "single-cohort-binary-biomarker-survival",
  help_printer = print_usage
)

feature_list <- cli_parse_csv_list(opt$features)
endpoint_list <- cli_parse_csv_list(opt$endpoints, uppercase = TRUE)
requested_covariates <- setdiff(cli_parse_csv_list(opt$adjust_covariates, allow_empty = TRUE), "none")

if (!length(feature_list)) stop("--features must name at least one binary biomarker.")
if (!nzchar(trimws(opt$positive_label)) || !nzchar(trimws(opt$negative_label))) {
  stop("--positive-label and --negative-label must be nonempty.")
}
if (identical(trimws(opt$positive_label), trimws(opt$negative_label))) {
  stop("--positive-label and --negative-label must be different.")
}

assert_output_guardrails(opt$output_dir)
cli_assert_existing_input_files(
  input_paths = c(biomarker_file = opt$biomarker_file, meta_file = opt$meta_file),
  workflow_name = "single-cohort-binary-biomarker-survival"
)

biomarker_data <- readr::read_tsv(opt$biomarker_file, show_col_types = FALSE)
if (!opt$feature_column %in% names(biomarker_data)) names(biomarker_data)[1] <- opt$feature_column
metadata <- read_metadata_table(opt$meta_file)
sample_columns <- validate_binary_matrix(biomarker_data, opt$feature_column, feature_list)

if (!opt$sample_column %in% names(metadata)) stop("Metadata missing sample column: ", opt$sample_column)
if (anyDuplicated(metadata[[opt$sample_column]])) stop("Metadata contains duplicate values in sample column: ", opt$sample_column)

if (is_filter_requested(opt$group_filter)) {
  if (!nzchar(trimws(opt$group_column))) stop("--group-filter requires --group-column.")
  if (!opt$group_column %in% names(metadata)) stop("Metadata missing group column: ", opt$group_column)
  group_values <- cli_parse_csv_list(opt$group_filter)
  metadata <- metadata %>% filter(.data[[opt$group_column]] %in% group_values)
  if (!nrow(metadata)) stop("No metadata rows remain after the requested group filter.")
}

covariates <- requested_covariates[requested_covariates %in% names(metadata)]
missing_covariates <- setdiff(requested_covariates, covariates)
if (length(missing_covariates)) {
  warning("Requested covariate columns not found in metadata and skipped: ", paste(missing_covariates, collapse = ", "))
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(
  data.frame(
    feature = feature_list,
    value_0_label = trimws(opt$negative_label),
    value_1_label = trimws(opt$positive_label),
    stringsAsFactors = FALSE
  ),
  file.path(opt$output_dir, "binary_biomarker_definition.csv")
)

manifest <- list()
for (endpoint_name in endpoint_list) {
  endpoint_columns <- resolve_endpoint_columns(endpoint_name, opt)
  endpoint_dir <- file.path(opt$output_dir, endpoint_name)
  dir.create(endpoint_dir, recursive = TRUE, showWarnings = FALSE)

  result <- analyze_survival_expression(
    expression_data = biomarker_data,
    metadata = metadata,
    sample_col = opt$sample_column,
    gene_col = opt$feature_column,
    genes = feature_list,
    event_col = endpoint_columns$event,
    time_col = endpoint_columns$time,
    clinical_covariates = covariates,
    split_method = "binary",
    include_groups = c("low", "high"),
    steps_to_run = c("individual", if (length(covariates)) "individual+covariates", "km_individual"),
    run_lasso = FALSE,
    run_ph_checks = is_truthy(opt$run_ph_checks),
    run_time_dependent = is_truthy(opt$run_time_dependent),
    save_plots = is_truthy(opt$save_plots),
    export_csv = is_truthy(opt$export_csv),
    output_dir = endpoint_dir,
    plot_params = list(
      plot_title = paste0(endpoint_name, " Binary Biomarker Survival"),
      binary_negative_label = trimws(opt$negative_label),
      binary_positive_label = trimws(opt$positive_label)
    ),
    verbose = is_truthy(opt$verbose)
  )

  if (is_truthy(opt$save_results_rds)) saveRDS(result, file.path(endpoint_dir, paste0(endpoint_name, "_results.rds")))
  write_binary_grouped_results(endpoint_name, endpoint_dir, result, trimws(opt$negative_label), trimws(opt$positive_label))
  write_binary_endpoint_summary(endpoint_name, endpoint_dir, result, feature_list, trimws(opt$negative_label), trimws(opt$positive_label), covariates)
  manifest[[length(manifest) + 1]] <- data.frame(endpoint = endpoint_name, endpoint_dir = endpoint_dir, metadata_rows = nrow(metadata))
}

readr::write_csv(dplyr::bind_rows(manifest), file.path(opt$output_dir, "run_manifest.csv"))
write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Run binary molecular-biomarker survival analysis through the reusable core engine.",
  driver = "drivers/run_single_cohort_binary_biomarker_survival.R",
  inputs = list(biomarker_file = opt$biomarker_file, meta_file = opt$meta_file),
  parameters = list(
    sample_column = opt$sample_column,
    feature_column = opt$feature_column,
    features = paste(feature_list, collapse = ","),
    binary_reference = paste0("0 = ", trimws(opt$negative_label)),
    binary_comparison = paste0("1 = ", trimws(opt$positive_label)),
    endpoints = paste(endpoint_list, collapse = ","),
    adjust_covariates = if (length(covariates)) paste(covariates, collapse = ",") else "none",
    group_column = opt$group_column,
    group_filter = opt$group_filter,
    run_ph_checks = opt$run_ph_checks,
    run_time_dependent = opt$run_time_dependent
  ),
  notes = "Input values are validated as 0, 1, or missing; missing values are excluded per biomarker.",
  scope = "exploratory"
)

message("Binary biomarker survival analysis written to: ", opt$output_dir)
