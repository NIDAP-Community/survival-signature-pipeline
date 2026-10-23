#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source("src/utilities/cli_helpers.R")
source("src/utilities/run_provenance_guardrails.R")
source("src/utilities/reporting_helpers.R")

cli_defaults <- list(
  input_run_dir = "",
  output_dir = default_run_output_dir("single_cohort_reporting_bundle"),
  signature_label = "Signature",
  subgroup_summary_table = "",
  subgroup_column = "subgroup",
  subgroup_type_column = "subgroup_type",
  feature_column = "feature",
  endpoint_column = "endpoint",
  comparison_column = "comparison",
  hr_column = "HR",
  ci_low_column = "CI_low",
  ci_high_column = "CI_high",
  p_value_column = "p_value",
  p_value_bh_column = "p_value_bh"
)

print_usage <- function() {
  print_standard_cli_help(
    title = "Single-Cohort Reporting Bundle Driver",
    workflow_name = "single-cohort-reporting-bundle",
    script_path = file.path("drivers", "run_single_cohort_reporting_bundle.R"),
    defaults = cli_defaults,
    summary_lines = c(
      "This workflow is downstream-only and expects a prior single-cohort run directory with run_manifest.csv and endpoint RDS outputs.",
      "It always builds an overall forest-style summary from available endpoint result objects.",
      "signature_label controls the displayed signature name in the overall forest plot and summary table.",
      "It only builds a heatmap when subgroup_summary_table is provided with explicit subgroup rows.",
      "The shared subgroup-summary contract uses subgroup, subgroup_type, feature, endpoint, comparison, HR, CI_low, CI_high, p_value, and p_value_bh."
    ),
    required_flags = c("input_run_dir"),
    aliases = c(subgroup_summary_csv = "subgroup_summary_table")
  )
}

extract_signature_summary <- function(signature_summary, endpoint_name, model_label, feature_label = "Signature") {
  if (is.null(signature_summary) || !inherits(signature_summary, "summary.coxph")) return(NULL)
  row_name <- if ("clusterC" %in% rownames(signature_summary$coefficients)) "clusterC" else rownames(signature_summary$coefficients)[1]
  if (is.na(row_name) || !nzchar(row_name)) return(NULL)

  data.frame(
    endpoint = endpoint_name,
    feature = feature_label,
    model = model_label,
    HR = signature_summary$conf.int[row_name, "exp(coef)"],
    CI_low = signature_summary$conf.int[row_name, "lower .95"],
    CI_high = signature_summary$conf.int[row_name, "upper .95"],
    p_value = signature_summary$coefficients[row_name, "Pr(>|z|)"],
    q_value = NA_real_,
    comparison = "High vs Low",
    stringsAsFactors = FALSE
  )
}

# Format a coefficient name from coxph without assuming any disease- or
# dataset-specific covariate names.  Factor coefficients are encoded by
# survival::coxph as <metadata_column><level>; the retained column names saved
# with the signature plot let us restore a readable "level vs reference" label.
format_adjusted_model_term <- function(term_name, signature_label, split_method, covariates) {
  display_name <- function(value) {
    value <- gsub("_", " ", value, fixed = TRUE)
    value <- gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", value, perl = TRUE)
    trimws(value)
  }

  if (identical(term_name, "clusterC")) {
    comparison <- if (identical(split_method, "detectable")) {
      "Detectable vs Undetectable"
    } else {
      "High vs Low"
    }
    return(paste0(signature_label, ": ", comparison))
  }
  if (identical(term_name, "clusterB")) {
    return(paste0(signature_label, ": Mid vs Low"))
  }

  covariates <- sort(unique(as.character(covariates)), decreasing = TRUE)
  matched_covariate <- covariates[startsWith(term_name, covariates)]
  if (length(matched_covariate)) {
    covariate <- matched_covariate[[1]]
    level <- substring(term_name, nchar(covariate) + 1L)
    if (nzchar(level)) {
      return(paste0(display_name(covariate), ": ", display_name(level), " vs reference"))
    }
    return(paste0(display_name(covariate), " (per unit)"))
  }

  # A model supplied by an older run may not retain covariate metadata. Keep
  # its coefficient visible rather than silently dropping it from the report.
  display_name(term_name)
}

extract_adjusted_model_terms <- function(signature_summary, endpoint_name, signature_label,
                                         split_method, covariates) {
  if (is.null(signature_summary) || !inherits(signature_summary, "summary.coxph")) return(NULL)

  coefficient_names <- rownames(signature_summary$coefficients)
  if (is.null(coefficient_names) || !length(coefficient_names)) return(NULL)

  data.frame(
    endpoint = endpoint_name,
    term = vapply(
      coefficient_names,
      format_adjusted_model_term,
      character(1),
      signature_label = signature_label,
      split_method = split_method,
      covariates = covariates
    ),
    raw_term = coefficient_names,
    effect_type = ifelse(grepl("^cluster", coefficient_names), "Signature", "Covariate"),
    HR = unname(signature_summary$conf.int[, "exp(coef)"]),
    CI_low = unname(signature_summary$conf.int[, "lower .95"]),
    CI_high = unname(signature_summary$conf.int[, "upper .95"]),
    p_value = unname(signature_summary$coefficients[, "Pr(>|z|)"]),
    N = signature_summary$n,
    events = signature_summary$nevent,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

extract_gene_table <- function(result_table, endpoint_name, model_label) {
  if (is.null(result_table) || !is.data.frame(result_table) || !nrow(result_table)) return(NULL)

  gene_column <- if ("Gene" %in% names(result_table)) "Gene" else names(result_table)[1]
  p_column <- if ("P_Value" %in% names(result_table)) "P_Value" else if ("p_value" %in% names(result_table)) "p_value" else NULL
  q_column <- if ("FDR" %in% names(result_table)) "FDR" else if ("p_value_bh" %in% names(result_table)) "p_value_bh" else NULL
  ci_low_column <- if ("CI_Lower" %in% names(result_table)) "CI_Lower" else if ("CI_low" %in% names(result_table)) "CI_low" else NULL
  ci_high_column <- if ("CI_Upper" %in% names(result_table)) "CI_Upper" else if ("CI_high" %in% names(result_table)) "CI_high" else NULL

  if (!("HR" %in% names(result_table)) || is.null(p_column) || is.null(ci_low_column) || is.null(ci_high_column)) return(NULL)

  result_table %>%
    transmute(
      endpoint = endpoint_name,
      feature = .data[[gene_column]],
      model = model_label,
      HR = .data[["HR"]],
      CI_low = .data[[ci_low_column]],
      CI_high = .data[[ci_high_column]],
      p_value = .data[[p_column]],
      q_value = if (!is.null(q_column)) .data[[q_column]] else NA_real_,
      comparison = "Continuous"
    )
}

# Use one data-driven presentation convention for every forest plot.  The
# effect text means a PNG remains interpretable after being copied out of its
# accompanying CSV/Markdown table, while the fixed familiar log ticks make OS,
# PFI, and other endpoints directly comparable.
format_forest_effect_label <- function(hr, ci_low, ci_high, p_value) {
  paste0(
    report_fmt_num(hr, 2), " [", report_fmt_num(ci_low, 2), ", ", report_fmt_num(ci_high, 2), "]",
    "; p=", report_fmt_p(p_value)
  )
}

forest_log_breaks <- function(lower, upper) {
  candidates <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.25, 0.5, 1, 2, 4, 5, 10, 20, 50, 100)
  selected <- candidates[candidates >= lower & candidates <= upper]
  unique(c(selected, if (lower <= 1 && upper >= 1) 1 else numeric(0)))
}

prepare_forest_display <- function(plot_df, label_column = "effect_label") {
  lower <- min(plot_df$CI_low, na.rm = TRUE)
  upper <- max(plot_df$CI_high, na.rm = TRUE)
  # Reserve a narrow right-hand column for effect estimates without placing the
  # text over confidence intervals.  This is logarithmic because HR is shown
  # on a log scale.
  text_x <- upper * 1.22
  x_upper <- text_x * 2
  x_lower <- lower / 1.15

  plot_df[[label_column]] <- format_forest_effect_label(
    plot_df$HR, plot_df$CI_low, plot_df$CI_high, plot_df$p_value
  )
  attr(plot_df, "forest_x_lower") <- x_lower
  attr(plot_df, "forest_x_upper") <- x_upper
  attr(plot_df, "forest_text_x") <- text_x
  attr(plot_df, "forest_breaks") <- forest_log_breaks(x_lower, upper)
  plot_df
}

build_overall_forest <- function(forest_df, output_stub) {
  plot_df <- forest_df %>%
    filter(is.finite(.data$CI_low), is.finite(.data$CI_high), .data$CI_low > 0, .data$CI_high > 0) %>%
    mutate(
      feature_label = factor(paste(.data$feature, .data$model, sep = " | "), levels = rev(unique(paste(.data$feature, .data$model, sep = " | ")))),
      model_type = case_when(
        grepl("unadjusted", .data$model, ignore.case = TRUE) ~ "Unadjusted",
        grepl("adjusted", .data$model, ignore.case = TRUE) ~ "Adjusted",
        TRUE ~ "Other"
      )
    )
  plot_df <- prepare_forest_display(plot_df)
  x_lower <- attr(plot_df, "forest_x_lower")
  x_upper <- attr(plot_df, "forest_x_upper")
  text_x <- attr(plot_df, "forest_text_x")
  x_breaks <- attr(plot_df, "forest_breaks")

  forest_plot <- ggplot(plot_df, aes(x = .data$HR, y = .data$feature_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(xmin = .data$CI_low, xmax = .data$CI_high), orientation = "y", width = 0.2, color = "grey35", linewidth = 0.5) +
    geom_point(aes(shape = .data$model_type, color = .data$model_type), size = 2.9, stroke = 0.7) +
    geom_text(aes(x = text_x, label = .data$effect_label), hjust = 0, size = 3.1, color = "grey20") +
    facet_grid(endpoint ~ ., scales = "free_y", space = "free_y") +
    scale_x_log10(breaks = x_breaks) +
    coord_cartesian(xlim = c(x_lower, x_upper), clip = "off") +
    scale_shape_manual(values = c("Unadjusted" = 1, "Adjusted" = 16, "Other" = 15), name = NULL) +
    scale_color_manual(values = c("Unadjusted" = "#5C5C5C", "Adjusted" = "#0072B2", "Other" = "#5C5C5C"), name = NULL) +
    labs(
      title = "Single-Cohort Overall Effects",
      subtitle = "Available gene and signature hazard ratios by endpoint; right labels show HR [95% CI] and p value",
      x = "Hazard ratio (log scale)",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(t = 5.5, r = 85, b = 5.5, l = 5.5),
      legend.position = "bottom"
    )

  plot_height <- max(5.5, 2.3 + 1.05 * max(table(plot_df$endpoint)))
  ggsave(paste0(output_stub, ".png"), forest_plot, width = 11.5, height = plot_height, dpi = 300)
  ggsave(paste0(output_stub, ".pdf"), forest_plot, width = 11.5, height = plot_height)
}

build_full_adjusted_model_forest <- function(model_terms_df, output_stub) {
  plot_df <- model_terms_df %>%
    filter(is.finite(.data$CI_low), is.finite(.data$CI_high), .data$CI_low > 0, .data$CI_high > 0) %>%
    mutate(
      endpoint = factor(.data$endpoint, levels = unique(.data$endpoint)),
      effect_type = factor(.data$effect_type, levels = c("Signature", "Covariate")),
      term_label = factor(.data$term, levels = rev(unique(.data$term)))
    )

  if (!nrow(plot_df)) {
    warning("No finite adjusted-model coefficient rows were available for the full forest plot.")
    return(invisible(NULL))
  }
  plot_df <- prepare_forest_display(plot_df)
  x_lower <- attr(plot_df, "forest_x_lower")
  x_upper <- attr(plot_df, "forest_x_upper")
  text_x <- attr(plot_df, "forest_text_x")
  x_breaks <- attr(plot_df, "forest_breaks")

  forest_plot <- ggplot(plot_df, aes(x = .data$HR, y = .data$term_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_errorbar(
      aes(xmin = .data$CI_low, xmax = .data$CI_high),
      orientation = "y",
      width = 0.2,
      color = "grey35",
      linewidth = 0.5
    ) +
    geom_point(aes(shape = .data$effect_type, color = .data$effect_type), size = 2.8, stroke = 0.7) +
    geom_text(aes(x = text_x, label = .data$effect_label), hjust = 0, size = 3.1, color = "grey20") +
    facet_grid(endpoint ~ ., scales = "free_y", space = "free_y") +
    scale_x_log10(breaks = x_breaks) +
    coord_cartesian(xlim = c(x_lower, x_upper), clip = "off") +
    scale_shape_manual(values = c("Signature" = 18, "Covariate" = 16), name = NULL) +
    scale_color_manual(values = c("Signature" = "#D55E00", "Covariate" = "#0072B2"), name = NULL) +
    labs(
      title = "Full Adjusted Cox Model",
      subtitle = "Signature comparison and every retained covariate; right labels show HR [95% CI] and p value",
      x = "Hazard ratio (log scale)",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold"),
      plot.margin = margin(t = 5.5, r = 85, b = 5.5, l = 5.5),
      legend.position = "bottom"
    )

  plot_height <- max(7.5, 2.5 + 0.32 * max(table(plot_df$endpoint)))
  ggsave(paste0(output_stub, ".png"), forest_plot, width = 11.5, height = plot_height, dpi = 300)
  ggsave(paste0(output_stub, ".pdf"), forest_plot, width = 11.5, height = plot_height)

  invisible(plot_df)
}

build_subgroup_heatmap <- function(subgroup_df, output_stub, subgroup_column, feature_column, endpoint_column, comparison_column, hr_column, p_value_column, p_value_bh_column) {
  plot_df <- prepare_subgroup_heatmap_data(
    subgroup_df = subgroup_df,
    subgroup_column = subgroup_column,
    feature_column = feature_column,
    endpoint_column = endpoint_column,
    hr_column = hr_column,
    p_value_column = p_value_column,
    p_value_bh_column = p_value_bh_column,
    comparison_column = comparison_column
  )

  heatmap_plot <- ggplot(plot_df, aes(x = .data$subgroup_label, y = .data$feature_label, fill = .data$log2_hr)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_point(
      data = plot_df %>% filter(.data$significance_label != "NS"),
      aes(shape = .data$significance_label),
      size = 2.4,
      stroke = 0.5,
      fill = "#ffd166",
      color = "black"
    ) +
    facet_grid(endpoint_label ~ comparison_label) +
    scale_shape_manual(values = c("p<0.05" = 21, "q<0.10" = 23), name = NULL) +
    subgroup_heatmap_fill_scale() +
    labs(
      title = "Single-Cohort Subgroup Effects",
      subtitle = "Optional subgroup heatmap generated from explicit subgroup summary rows",
      x = NULL,
      y = NULL
    ) +
    subgroup_heatmap_theme(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(paste0(output_stub, ".png"), heatmap_plot, width = 10, height = 6, dpi = 300)
  ggsave(paste0(output_stub, ".pdf"), heatmap_plot, width = 10, height = 6)
}

opt <- parse_standard_cli_args(
  defaults = cli_defaults,
  workflow_name = "single-cohort-reporting-bundle",
  help_printer = print_usage,
  required_flags = c("input_run_dir"),
  aliases = c(subgroup_summary_csv = "subgroup_summary_table")
)

assert_output_guardrails(opt$output_dir)

manifest_path <- file.path(opt$input_run_dir, "run_manifest.csv")
if (!file.exists(manifest_path)) {
  stop("input_run_dir is missing run_manifest.csv: ", manifest_path)
}

manifest_df <- readr::read_csv(manifest_path, show_col_types = FALSE)
if (!all(c("endpoint", "endpoint_dir") %in% names(manifest_df))) {
  stop("run_manifest.csv must contain endpoint and endpoint_dir columns.")
}

ensure_reporting_output_dirs(opt$output_dir)

forest_rows <- list()
adjusted_model_term_rows <- list()

for (row_index in seq_len(nrow(manifest_df))) {
  endpoint_name <- as.character(manifest_df$endpoint[[row_index]])
  endpoint_dir <- as.character(manifest_df$endpoint_dir[[row_index]])
  result_path <- file.path(endpoint_dir, paste0(endpoint_name, "_results.rds"))
  if (!file.exists(result_path)) {
    stop("Missing endpoint result object: ", result_path)
  }

  result_obj <- readRDS(result_path)

  forest_rows[[length(forest_rows) + 1]] <- extract_gene_table(result_obj$results_univariate, endpoint_name, "gene_unadjusted")
  forest_rows[[length(forest_rows) + 1]] <- extract_gene_table(result_obj$results_adjusted, endpoint_name, "gene_adjusted")
  forest_rows[[length(forest_rows) + 1]] <- extract_signature_summary(
    result_obj$results_signature$unadjusted,
    endpoint_name,
    "Unadjusted",
    feature_label = opt$signature_label
  )
  forest_rows[[length(forest_rows) + 1]] <- extract_signature_summary(
    result_obj$results_signature$adjusted,
    endpoint_name,
    "Adjusted",
    feature_label = opt$signature_label
  )

  signature_plot <- result_obj$plots$signature
  adjusted_model_term_rows[[length(adjusted_model_term_rows) + 1]] <- extract_adjusted_model_terms(
    signature_summary = result_obj$results_signature$adjusted,
    endpoint_name = endpoint_name,
    signature_label = opt$signature_label,
    split_method = signature_plot$split_method %||% "quartiles",
    covariates = signature_plot$covariates_retained_grouped %||% character(0)
  )
}

forest_df <- bind_rows(Filter(Negate(is.null), forest_rows))
if (!nrow(forest_df)) {
  stop("No overall effect rows could be extracted from the single-cohort run outputs.")
}

forest_md <- forest_df %>%
  transmute(
    endpoint = endpoint,
    feature = feature,
    model = model,
    comparison = comparison,
    HR = report_fmt_num(HR, 3),
    `95% CI` = report_fmt_ci(CI_low, CI_high),
    p = report_fmt_p(p_value),
    q = report_fmt_p(q_value)
  )

write_csv_and_md_table(
  csv_df = forest_df,
  csv_path = file.path(opt$output_dir, "tables", "single_cohort_overall_effects.csv"),
  md_path = file.path(opt$output_dir, "tables", "single_cohort_overall_effects.md"),
  md_df = forest_md
)

build_overall_forest(forest_df, file.path(opt$output_dir, "figures", "figure1_single_cohort_overall_forest"))

adjusted_model_terms_df <- bind_rows(Filter(Negate(is.null), adjusted_model_term_rows))
if (nrow(adjusted_model_terms_df)) {
  adjusted_model_terms_md <- adjusted_model_terms_df %>%
    transmute(
      endpoint = endpoint,
      term = term,
      type = effect_type,
      `raw coefficient` = raw_term,
      HR = report_fmt_num(HR, 3),
      `95% CI` = report_fmt_ci(CI_low, CI_high),
      p = report_fmt_p(p_value),
      N = N,
      events = events
    )
  write_csv_and_md_table(
    csv_df = adjusted_model_terms_df,
    csv_path = file.path(opt$output_dir, "tables", "single_cohort_full_adjusted_model_terms.csv"),
    md_path = file.path(opt$output_dir, "tables", "single_cohort_full_adjusted_model_terms.md"),
    md_df = adjusted_model_terms_md
  )
  build_full_adjusted_model_forest(
    adjusted_model_terms_df,
    file.path(opt$output_dir, "figures", "figure2_single_cohort_full_adjusted_model_forest")
  )
}

notes <- character(0)

if (nzchar(opt$subgroup_summary_table)) {
  if (!file.exists(opt$subgroup_summary_table)) {
    stop("subgroup_summary_table does not exist: ", opt$subgroup_summary_table)
  }
  subgroup_df <- readr::read_csv(opt$subgroup_summary_table, show_col_types = FALSE)
  normalized_subgroup_df <- normalize_subgroup_summary_table(
    subgroup_df = subgroup_df,
    subgroup_column = opt$subgroup_column,
    subgroup_type_column = opt$subgroup_type_column,
    feature_column = opt$feature_column,
    endpoint_column = opt$endpoint_column,
    comparison_column = opt$comparison_column,
    hr_column = opt$hr_column,
    ci_low_column = opt$ci_low_column,
    ci_high_column = opt$ci_high_column,
    p_value_column = opt$p_value_column,
    p_value_bh_column = opt$p_value_bh_column
  )
  validate_subgroup_summary_table(
    subgroup_df = normalized_subgroup_df,
    context_label = "subgroup_summary_table"
  )
  readr::write_csv(normalized_subgroup_df, file.path(opt$output_dir, "tables", "consumed_subgroup_summary.csv"))
  build_subgroup_heatmap(
    normalized_subgroup_df,
    file.path(opt$output_dir, "figures", "figure2_single_cohort_subgroup_heatmap"),
    "subgroup",
    "feature",
    "endpoint",
    "comparison",
    "HR",
    "p_value",
    "p_value_bh"
  )
} else {
  notes <- c(notes, "No subgroup_summary_csv was provided, so the subgroup heatmap was intentionally skipped.")
  notes <- sub("subgroup_summary_csv", "subgroup_summary_table", notes, fixed = TRUE)
  write_reporting_text_block(
    lines = c(
      "# Heatmap Skipped",
      "",
      "No subgroup summary table was provided.",
      "The single-cohort reporting bundle therefore generated its available overall forest figures only."
    ),
    file_path = file.path(opt$output_dir, "figures", "figure3_single_cohort_subgroup_heatmap_skipped.md")
  )
}

write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Assemble single-cohort reporting artifacts from a prior single-cohort analysis run.",
  driver = "drivers/run_single_cohort_reporting_bundle.R",
  inputs = list(
    input_run_dir = opt$input_run_dir,
    subgroup_summary_table = if (nzchar(opt$subgroup_summary_table)) opt$subgroup_summary_table else ""
  ),
  parameters = as.list(opt),
  notes = notes,
  scope = "exploratory"
)

message("Single-cohort reporting bundle written to: ", opt$output_dir)
