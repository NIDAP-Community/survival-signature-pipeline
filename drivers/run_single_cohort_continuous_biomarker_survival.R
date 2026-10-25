#!/usr/bin/env Rscript

# Maintained workflow for a prepared, patient-level continuous molecular-feature
# matrix. Copy-number log2 ratios and immune-score matrices are intended use
# cases, but the driver deliberately does not calculate source measurements.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("src/core/GeneExpressionBasedSurvivalAnalys.r")
source("src/utilities/cli_helpers.R")
source("src/utilities/reporting_helpers.R")
source("src/utilities/run_provenance_guardrails.R")

cli_defaults <- list(
  output_dir = default_run_output_dir("single_cohort_continuous_biomarker_survival"),
  biomarker_file = file.path("test_data", "continuous_biomarker_matrix.tsv"),
  meta_file = file.path("test_data", "tcga_metadata_no_subtype.csv"),
  sample_column = "sampnames",
  feature_column = "Feature",
  features = "MYC_log2_copy_number",
  biomarker_context = "generic",
  fraction_sum_tolerance = "0.02",
  split_method = "tertiles",
  include_groups = "low,mid,high",
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
    title = "Single-Cohort Continuous Biomarker Survival Driver",
    workflow_name = "single-cohort-continuous-biomarker-survival",
    script_path = file.path("drivers", "run_single_cohort_continuous_biomarker_survival.R"),
    defaults = cli_defaults,
    summary_lines = c(
      "Use this workflow for a prepared feature-by-subject numeric matrix, including log2 copy-number values or immune scores.",
      "Continuous Cox results are reported per within-cohort standard deviation; low/mid/high KM groups are display strata only.",
      "biomarker_context=estimate_scores records a prepared ESTIMATE-score run; biomarker_context=cell_fractions requires a complete 0-1 composition that sums to one per subject.",
      "The driver does not query a database, normalize copy number, calculate immune scores, or infer biological cutoffs.",
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
  resolve_standard_endpoint_columns(endpoint_name, options, "single-cohort-continuous-biomarker-survival")
}

validate_continuous_matrix <- function(biomarker_data, feature_column, selected_features,
                                       biomarker_context, fraction_sum_tolerance) {
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

  raw_values <- as.matrix(biomarker_data[, sample_columns, drop = FALSE])
  numeric_values <- suppressWarnings(as.numeric(raw_values))
  invalid <- !is.na(raw_values) & (!is.finite(numeric_values))
  if (any(invalid)) {
    stop("Continuous biomarker values must be finite numeric values or missing.")
  }

  for (column in sample_columns) {
    biomarker_data[[column]] <- suppressWarnings(as.numeric(biomarker_data[[column]]))
  }

  if (identical(biomarker_context, "cell_fractions")) {
    fraction_matrix <- as.matrix(biomarker_data[, sample_columns, drop = FALSE])
    outside_unit_interval <- !is.na(fraction_matrix) & (fraction_matrix < 0 | fraction_matrix > 1)
    if (any(outside_unit_interval)) {
      stop("cell_fractions context requires every non-missing value to be between 0 and 1.")
    }
    complete_subjects <- colSums(is.na(fraction_matrix)) == 0
    if (!all(complete_subjects)) {
      incomplete_subjects <- colnames(fraction_matrix)[!complete_subjects]
      stop(
        "cell_fractions context requires every staged component for every subject. Incomplete subjects: ",
        paste(head(incomplete_subjects, 10), collapse = ", ")
      )
    }
    composition_totals <- colSums(fraction_matrix)
    bad_subjects <- names(composition_totals)[abs(composition_totals - 1) > fraction_sum_tolerance]
    if (length(bad_subjects)) {
      stop(
        "cell_fractions context requires complete subject-wise fractions to sum to 1 within --fraction-sum-tolerance. Failed subjects: ",
        paste(head(bad_subjects, 10), collapse = ", ")
      )
    }
  }

  selected_numeric <- biomarker_data[biomarker_data[[feature_column]] %in% selected_features, sample_columns, drop = FALSE]
  feature_sd <- apply(selected_numeric, 1, stats::sd, na.rm = TRUE)
  if (any(!is.finite(feature_sd) | feature_sd == 0)) {
    bad_features <- biomarker_data[[feature_column]][biomarker_data[[feature_column]] %in% selected_features][!is.finite(feature_sd) | feature_sd == 0]
    stop("Requested features must have nonzero finite variation: ", paste(bad_features, collapse = ", "))
  }
  biomarker_data
}

write_continuous_endpoint_summary <- function(endpoint_name, endpoint_dir, result, feature_list,
                                              split_method, covariates, biomarker_context) {
  grouped <- result$results_grouped_gene
  lines <- c(
    paste0("# ", endpoint_name, " Continuous Biomarker Survival Summary"),
    "",
    paste0("- Endpoint: ", endpoint_name),
    paste0("- Biomarker context: ", gsub("_", " ", biomarker_context)),
    paste0("- Kaplan-Meier grouping: ", split_method),
    "- Continuous Cox unit: one within-cohort standard deviation of the submitted feature.",
    paste0("- Features: ", paste(feature_list, collapse = ", ")),
    paste0("- Requested covariates: ", if (length(covariates)) paste(covariates, collapse = ", ") else "none"),
    "",
    "Kaplan-Meier curves and log-rank testing are unadjusted displays. Continuous Cox results estimate the primary standardized feature effect; grouped Cox rows compare the displayed strata.",
    if (identical(biomarker_context, "cell_fractions")) "Cell fractions are compositional: these are one-feature exploratory associations, not a multivariable compositional survival model." else NULL,
    ""
  )
  if (!is.null(grouped) && nrow(grouped)) {
    display <- grouped %>%
      transmute(
        Feature = feature,
        Model = model_type,
        Comparison = comparison,
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

write_continuous_grouped_results <- function(endpoint_name, endpoint_dir, result) {
  grouped <- result$results_grouped_gene
  if (is.null(grouped) || !nrow(grouped)) return(invisible(NULL))
  output <- grouped %>%
    mutate(endpoint = endpoint_name) %>%
    select(endpoint, feature, comparison, reference_group, comparison_group,
           HR, CI_lower, CI_upper, p_value, N, events, model_type,
           covariates_requested, covariates_retained)
  readr::write_csv(output, file.path(endpoint_dir, "grouped_biomarker_cox_results.csv"))
}

opt <- parse_standard_cli_args(
  defaults = cli_defaults,
  workflow_name = "single-cohort-continuous-biomarker-survival",
  help_printer = print_usage
)

feature_list <- cli_parse_csv_list(opt$features)
endpoint_list <- cli_parse_csv_list(opt$endpoints, uppercase = TRUE)
requested_covariates <- cli_parse_csv_list(opt$adjust_covariates, allow_empty = TRUE)
requested_covariates <- requested_covariates[tolower(requested_covariates) != "none"]
include_groups <- cli_parse_csv_list(opt$include_groups, lowercase = TRUE)
opt$biomarker_context <- cli_parse_enum_flag(
  opt$biomarker_context,
  allowed_values = c("generic", "estimate_scores", "cell_fractions"),
  flag_name = "biomarker_context",
  normalize = "lower"
)
opt$fraction_sum_tolerance <- cli_parse_numeric_flag(opt$fraction_sum_tolerance, "fraction_sum_tolerance")
opt$split_method <- cli_parse_enum_flag(
  opt$split_method,
  allowed_values = c("tertiles", "quartiles", "median"),
  flag_name = "split_method",
  normalize = "lower"
)

if (!length(feature_list)) stop("--features must name at least one continuous biomarker.")
if (anyDuplicated(feature_list)) stop("--features must not repeat a biomarker identifier.")
if (!all(include_groups %in% c("low", "mid", "high"))) stop("--include-groups must contain only low, mid, and/or high.")
if (opt$split_method == "median" && "mid" %in% include_groups) {
  stop("split_method=median supports only low and high groups.")
}
if (opt$fraction_sum_tolerance < 0 || opt$fraction_sum_tolerance >= 1) {
  stop("--fraction-sum-tolerance must be at least 0 and less than 1.")
}

assert_output_guardrails(opt$output_dir)
cli_assert_existing_input_files(
  input_paths = c(biomarker_file = opt$biomarker_file, meta_file = opt$meta_file),
  workflow_name = "single-cohort-continuous-biomarker-survival"
)

biomarker_data <- readr::read_tsv(opt$biomarker_file, show_col_types = FALSE)
if (!opt$feature_column %in% names(biomarker_data)) names(biomarker_data)[1] <- opt$feature_column
metadata <- read_metadata_table(opt$meta_file)
biomarker_data <- validate_continuous_matrix(
  biomarker_data,
  opt$feature_column,
  feature_list,
  opt$biomarker_context,
  opt$fraction_sum_tolerance
)

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
    biomarker_context = opt$biomarker_context,
    continuous_cox_unit = "one within-cohort standard deviation",
    km_split_method = opt$split_method,
    stringsAsFactors = FALSE
  ),
  file.path(opt$output_dir, "continuous_biomarker_definition.csv")
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
    split_method = opt$split_method,
    include_groups = include_groups,
    steps_to_run = c("individual", if (length(covariates)) "individual+covariates", "km_individual"),
    run_lasso = FALSE,
    run_ph_checks = is_truthy(opt$run_ph_checks),
    run_time_dependent = is_truthy(opt$run_time_dependent),
    save_plots = is_truthy(opt$save_plots),
    export_csv = is_truthy(opt$export_csv),
    output_dir = endpoint_dir,
    plot_params = list(
      plot_title = paste0(
        endpoint_name,
        " ",
        if (opt$biomarker_context == "estimate_scores") "ESTIMATE Score" else if (opt$biomarker_context == "cell_fractions") "Cell Fraction" else "Continuous Biomarker",
        " Survival"
      )
    ),
    verbose = is_truthy(opt$verbose)
  )

  if (is_truthy(opt$save_results_rds)) saveRDS(result, file.path(endpoint_dir, paste0(endpoint_name, "_results.rds")))
  write_continuous_grouped_results(endpoint_name, endpoint_dir, result)
  write_continuous_endpoint_summary(endpoint_name, endpoint_dir, result, feature_list, opt$split_method, covariates, opt$biomarker_context)
  manifest[[length(manifest) + 1]] <- data.frame(endpoint = endpoint_name, endpoint_dir = endpoint_dir, metadata_rows = nrow(metadata))
}

readr::write_csv(dplyr::bind_rows(manifest), file.path(opt$output_dir, "run_manifest.csv"))
write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Run continuous molecular-biomarker survival analysis through the reusable core engine.",
  driver = "drivers/run_single_cohort_continuous_biomarker_survival.R",
  inputs = list(biomarker_file = opt$biomarker_file, meta_file = opt$meta_file),
  parameters = list(
    sample_column = opt$sample_column,
    feature_column = opt$feature_column,
    features = paste(feature_list, collapse = ","),
    biomarker_context = opt$biomarker_context,
    fraction_sum_tolerance = opt$fraction_sum_tolerance,
    continuous_cox_unit = "one within-cohort standard deviation",
    km_split_method = opt$split_method,
    include_groups = paste(include_groups, collapse = ","),
    endpoints = paste(endpoint_list, collapse = ","),
    adjust_covariates = if (length(covariates)) paste(covariates, collapse = ",") else "none",
    group_column = opt$group_column,
    group_filter = opt$group_filter,
    run_ph_checks = opt$run_ph_checks,
    run_time_dependent = opt$run_time_dependent
  ),
  notes = if (opt$biomarker_context == "cell_fractions") {
    "Input is a complete 0-1 composition within the supplied tolerance; results are exploratory one-feature associations, not a multivariable compositional model."
  } else if (opt$biomarker_context == "estimate_scores") {
    "Input is a prepared ESTIMATE-score matrix; score calculation and provenance remain upstream."
  } else {
    "Input values are finite numeric measurements or missing; no amplification/deletion threshold is inferred."
  },
  scope = "exploratory"
)

message("Continuous biomarker survival analysis written to: ", opt$output_dir)
