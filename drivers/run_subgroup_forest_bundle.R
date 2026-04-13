#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("src/utilities/cli_helpers.R")
source("src/utilities/run_provenance_guardrails.R")
source("src/utilities/reporting_helpers.R")

cli_defaults <- list(
  subgroup_summary_table = "",
  output_dir = default_run_output_dir("subgroup_forest_bundle"),
  subgroup_column = "subgroup",
  subgroup_type_column = "subgroup_type",
  feature_column = "feature",
  endpoint_column = "endpoint",
  comparison_column = "comparison",
  hr_column = "HR",
  ci_low_column = "CI_low",
  ci_high_column = "CI_high",
  p_value_column = "p_value",
  p_value_bh_column = "p_value_bh",
  subgroup_type_filter = "",
  features = "",
  endpoints = "",
  comparisons = "",
  subgroup_order = "",
  subgroup_labels = "",
  q_threshold = "0.10",
  p_threshold = "0.05",
  plot_title = "Subgroup Effects",
  plot_subtitle = "Hazard ratios across subgroup-defined analyses"
)

print_usage <- function() {
  print_standard_cli_help(
    title = "Subgroup Forest Bundle Driver",
    workflow_name = "subgroup-forest-bundle",
    script_path = file.path("drivers", "run_subgroup_forest_bundle.R"),
    defaults = cli_defaults,
    summary_lines = c(
      "This workflow is downstream-only and consumes the shared subgroup-summary contract.",
      "Use subgroup_type_filter to select one cohorting scheme such as subtype, response_status, or gleason_risk.",
      "subgroup_order and subgroup_labels are optional CSV lists for custom axis ordering and display labels.",
      "The figure is generic: rows are subgroups, points are hazard ratios, and panels are split by endpoint, comparison, and feature."
    ),
    required_flags = c("subgroup_summary_table")
  )
}

resolve_subgroup_label_map <- function(subgroup_order, subgroup_labels) {
  if (!length(subgroup_labels)) {
    return(NULL)
  }
  if (!length(subgroup_order)) {
    stop("--subgroup-labels requires --subgroup-order so labels can be mapped deterministically.")
  }
  if (length(subgroup_order) != length(subgroup_labels)) {
    stop("--subgroup-order and --subgroup-labels must contain the same number of values.")
  }
  stats::setNames(subgroup_labels, subgroup_order)
}

opt <- parse_standard_cli_args(
  defaults = cli_defaults,
  workflow_name = "subgroup-forest-bundle",
  help_printer = print_usage,
  required_flags = c("subgroup_summary_table")
)

opt$q_threshold <- cli_parse_numeric_flag(opt$q_threshold, "q_threshold")
opt$p_threshold <- cli_parse_numeric_flag(opt$p_threshold, "p_threshold")

assert_output_guardrails(opt$output_dir)
ensure_reporting_output_dirs(opt$output_dir)
cli_assert_existing_input_files(
  input_paths = c(subgroup_summary_table = opt$subgroup_summary_table),
  workflow_name = "subgroup-forest-bundle"
)

subgroup_input_df <- readr::read_csv(opt$subgroup_summary_table, show_col_types = FALSE)
normalized_subgroup_df <- normalize_subgroup_summary_table(
  subgroup_df = subgroup_input_df,
  subgroup_column = opt$subgroup_column,
  subgroup_type_column = opt$subgroup_type_column,
  feature_column = opt$feature_column,
  endpoint_column = opt$endpoint_column,
  comparison_column = opt$comparison_column,
  hr_column = opt$hr_column,
  ci_low_column = opt$ci_low_column,
  ci_high_column = opt$ci_high_column,
  p_value_column = opt$p_value_column,
  p_value_bh_column = opt$p_value_bh_column,
  allow_legacy_pooled_columns = TRUE
)
validate_subgroup_summary_table(
  subgroup_df = normalized_subgroup_df,
  context_label = "subgroup_summary_table",
  require_dataset_window = FALSE
)

if (nzchar(opt$subgroup_type_filter)) {
  subgroup_type_values <- cli_parse_csv_list(opt$subgroup_type_filter)
  normalized_subgroup_df <- normalized_subgroup_df %>% filter(.data$subgroup_type %in% subgroup_type_values)
}
if (nzchar(opt$features)) {
  feature_values <- cli_parse_csv_list(opt$features)
  normalized_subgroup_df <- normalized_subgroup_df %>% filter(.data$feature %in% feature_values)
}
if (nzchar(opt$endpoints)) {
  endpoint_values <- cli_parse_csv_list(opt$endpoints, uppercase = TRUE)
  normalized_subgroup_df <- normalized_subgroup_df %>% filter(toupper(.data$endpoint) %in% endpoint_values)
}
if (nzchar(opt$comparisons)) {
  comparison_values <- cli_parse_csv_list(opt$comparisons)
  normalized_subgroup_df <- normalized_subgroup_df %>% filter(.data$comparison %in% comparison_values)
}

if (!nrow(normalized_subgroup_df)) {
  stop("No subgroup rows remain after applying the requested filters.")
}

subgroup_order <- cli_parse_csv_list(opt$subgroup_order, allow_empty = TRUE)
if (!length(subgroup_order)) {
  subgroup_order <- unique(as.character(normalized_subgroup_df$subgroup))
}
subgroup_labels <- cli_parse_csv_list(opt$subgroup_labels, allow_empty = TRUE)
subgroup_label_map <- resolve_subgroup_label_map(subgroup_order, subgroup_labels)
subgroup_transform <- if (is.null(subgroup_label_map)) {
  identity
} else {
  function(values) {
    mapped <- subgroup_label_map[as.character(values)]
    mapped[is.na(mapped)] <- as.character(values)[is.na(mapped)]
    unname(mapped)
  }
}

plot_df <- build_subgroup_forest(
  subgroup_df = normalized_subgroup_df,
  output_stub = file.path(opt$output_dir, "figures", "figure1_subgroup_forest"),
  subgroup_column = "subgroup",
  feature_column = "feature",
  endpoint_column = "endpoint",
  comparison_column = "comparison",
  hr_column = "HR",
  ci_low_column = "CI_low",
  ci_high_column = "CI_high",
  p_value_column = "p_value",
  p_value_bh_column = "p_value_bh",
  subgroup_transform = subgroup_transform,
  subgroup_levels = subgroup_transform(subgroup_order),
  feature_levels = unique(as.character(normalized_subgroup_df$feature)),
  endpoint_levels = unique(as.character(normalized_subgroup_df$endpoint)),
  comparison_levels = unique(as.character(normalized_subgroup_df$comparison)),
  q_threshold = opt$q_threshold,
  p_threshold = opt$p_threshold,
  title = opt$plot_title,
  subtitle = opt$plot_subtitle
)

table_df <- plot_df %>%
  transmute(
    subgroup = as.character(.data$subgroup_label),
    subgroup_type = as.character(normalized_subgroup_df$subgroup_type[match(as.character(.data$subgroup_label), subgroup_transform(as.character(normalized_subgroup_df$subgroup)))]),
    feature = as.character(.data$feature_label),
    endpoint = as.character(.data$endpoint_label),
    comparison = as.character(.data$comparison_label),
    HR = .data$HR,
    CI_low = .data$CI_low,
    CI_high = .data$CI_high,
    p_value = .data$p_value,
    p_value_bh = .data$p_value_bh,
    significance = .data$significance_label
  )

table_md <- table_df %>%
  mutate(
    HR = report_fmt_num(HR, 3),
    `95% CI` = report_fmt_ci(CI_low, CI_high),
    p = report_fmt_p(p_value),
    q = report_fmt_p(p_value_bh)
  ) %>%
  select(subgroup, subgroup_type, feature, endpoint, comparison, HR, `95% CI`, p, q, significance)

write_csv_and_md_table(
  csv_df = table_df,
  csv_path = file.path(opt$output_dir, "tables", "consumed_subgroup_forest_rows.csv"),
  md_path = file.path(opt$output_dir, "tables", "consumed_subgroup_forest_rows.md"),
  md_df = table_md
)

write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Build generic subgroup forest plots from a shared subgroup-summary table.",
  driver = "drivers/run_subgroup_forest_bundle.R",
  inputs = list(subgroup_summary_table = opt$subgroup_summary_table),
  parameters = as.list(opt),
  scope = "exploratory"
)

message("Subgroup forest bundle written to: ", opt$output_dir)