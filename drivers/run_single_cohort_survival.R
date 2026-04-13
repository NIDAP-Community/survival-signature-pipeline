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
  group_column = "",
  group_filter = "ALL",
  group_type = "subgroup",
  subtype_column = "Subtype",
  stage_column = "ajcc_pathologic_tumor_stage",
  age_column = "Diagnosis_Age",
  race_column = "race_collapsed",
  signature_genes = "GENE_A,GENE_B,GENE_C,GENE_D",
  os_event_column = "OS",
  os_time_column = "OS_time",
  pfi_event_column = "PFI",
  pfi_time_column = "PFI_time",
  subtype_filter = "ALL",
  adjust_covariates = "none",
  split_method = "tertiles",
  include_groups = "low,mid,high",
  features = "Signature,GENE_A,GENE_B,GENE_C,GENE_D",
  endpoints = "OS,PFI",
  cut_months = "60",
  run_km_individual = "true",
  save_plots = "true",
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
      "group_column/group_filter provide the maintained generic cohorting interface for subgroup analyses.",
      "subtype_column/subtype_filter are retained as compatibility flags for subtype-shaped datasets.",
      "Inputs are expected to be post-filtered upstream; the maintained driver does not apply an internal primary-only sample filter.",
      "adjust_covariates should generally be supplied as explicit metadata column names for the current dataset.",
      "Optional aliases are supported as a convenience only: subtype, stage, age, race.",
      "Requested covariates missing from a dataset are skipped automatically; active group filtering requires the chosen group column."
    ),
    aliases = c(
      subgroup_column = "group_column",
      subgroup_filter = "group_filter"
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

parse_covariate_request <- function(value, subtype_column, stage_column, age_column, race_column) {
  alias_map <- c(
    subtype = subtype_column,
    stage = stage_column,
    age = age_column,
    race = race_column
  )
  legacy_pattern <- "^(subtype|stage|age|race)(_(subtype|stage|age|race))*$"
  tokens <- cli_parse_csv_list(value, allow_empty = TRUE)
  if (!length(tokens)) {
    return(list(requested = character(0), columns = character(0), used_aliases = character(0)))
  }

  requested_columns <- character(0)
  used_aliases <- character(0)

  for (token in tokens) {
    token_lower <- tolower(token)
    if (token_lower == "none") {
      next
    }

    if (token_lower %in% names(alias_map)) {
      requested_columns <- c(requested_columns, alias_map[[token_lower]])
      used_aliases <- c(used_aliases, token_lower)
      next
    }

    if (grepl(legacy_pattern, token_lower)) {
      parts <- unlist(strsplit(token_lower, "_", fixed = TRUE))
      requested_columns <- c(requested_columns, unname(alias_map[parts]))
      used_aliases <- c(used_aliases, parts)
      next
    }

    requested_columns <- c(requested_columns, token)
  }

  list(
    requested = tokens,
    columns = unique(requested_columns[nzchar(requested_columns)]),
    used_aliases = unique(used_aliases)
  )
}

write_simple_markdown <- function(lines, file_path) {
  writeLines(lines, con = file_path)
}

filter_requested <- function(value) {
  !(tolower(trimws(as.character(value))) %in% c("", "all", "overall", "any"))
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
  help_printer = print_usage,
  aliases = c(
    subgroup_column = "group_column",
    subgroup_filter = "group_filter"
  )
)
opt$cut_months <- cli_parse_numeric_flag(opt$cut_months, "cut_months")
feature_list <- cli_parse_csv_list(opt$features)
endpoint_list <- cli_parse_csv_list(opt$endpoints, uppercase = TRUE)
signature_genes <- cli_parse_csv_list(opt$signature_genes)
include_groups <- cli_parse_csv_list(opt$include_groups, lowercase = TRUE)
covariate_request <- parse_covariate_request(
  value = opt$adjust_covariates,
  subtype_column = opt$subtype_column,
  stage_column = opt$stage_column,
  age_column = opt$age_column,
  race_column = opt$race_column
)
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

group_filter_requested <- filter_requested(opt$group_filter)
subtype_filter_requested <- filter_requested(opt$subtype_filter)
group_column_requested <- nzchar(trimws(opt$group_column))
subtype_adjustment_requested <- any(covariate_request$used_aliases == "subtype") || opt$subtype_column %in% covariate_request$columns
has_subtype_column <- opt$subtype_column %in% names(metadata)

if (group_filter_requested && !group_column_requested) {
  stop("--group-filter requires --group-column.")
}

if (group_filter_requested && subtype_filter_requested) {
  stop("Use either --group-filter/--group-column or --subtype-filter, not both in the same run.")
}

effective_group_column <- NULL
effective_group_values <- character(0)
effective_group_type <- if (nzchar(trimws(opt$group_type))) trimws(opt$group_type) else "subgroup"

if (group_filter_requested) {
  effective_group_column <- trimws(opt$group_column)
  effective_group_values <- cli_parse_csv_list(opt$group_filter)
} else if (subtype_filter_requested) {
  effective_group_column <- opt$subtype_column
  effective_group_values <- cli_parse_csv_list(opt$subtype_filter)
  effective_group_type <- "subtype"
}

if (!has_subtype_column && subtype_filter_requested) {
  stop(
    "Metadata missing subtype column '", opt$subtype_column,
    "'. subtype_column is optional unless subtype filtering or subtype-adjusted modeling is requested. ",
    "Either provide --subtype-column for metadata that includes subtypes, or remove --subtype-filter."
  )
}

if (!is.null(effective_group_column) && !effective_group_column %in% names(metadata)) {
  stop("Metadata missing group column: ", effective_group_column)
}

metadata_filtered <- metadata

if (!is.null(effective_group_column) && length(effective_group_values)) {
  metadata_filtered <- metadata_filtered %>% filter(.data[[effective_group_column]] %in% effective_group_values)
}

if (!nrow(metadata_filtered)) {
  stop("No metadata rows remain after filtering. Check the active group filter and upstream dataset preparation.")
}

dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

covariates_requested <- covariate_request$columns
covariates <- covariates_requested[covariates_requested %in% names(metadata_filtered)]
missing_covariates <- setdiff(covariates_requested, covariates)

if (length(missing_covariates)) {
  warning(
    "Requested covariate columns not found in metadata and skipped: ",
    paste(missing_covariates, collapse = ", ")
  )
}

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
    group_column = if (!is.null(effective_group_column)) effective_group_column else opt$group_column,
    group_filter = if (length(effective_group_values)) paste(effective_group_values, collapse = ",") else opt$group_filter,
    group_type = effective_group_type,
    subtype_column = opt$subtype_column,
    stage_column = opt$stage_column,
    age_column = opt$age_column,
    race_column = opt$race_column,
    signature_genes = paste(signature_genes, collapse = ","),
    features = opt$features,
    endpoints = opt$endpoints,
    analysis_genes = paste(analysis_genes, collapse = ","),
    subtype_filter = opt$subtype_filter,
    subgroup_column = if (!is.null(effective_group_column)) effective_group_column else opt$group_column,
    subgroup_filter = if (length(effective_group_values)) paste(effective_group_values, collapse = ",") else opt$group_filter,
    adjust_covariates = opt$adjust_covariates,
    covariates_requested = paste(covariates_requested, collapse = ","),
    covariates_resolved = paste(covariates, collapse = ","),
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