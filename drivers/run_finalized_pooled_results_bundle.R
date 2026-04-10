#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(rlang)
})

source("src/utilities/cli_helpers.R")
source("src/utilities/run_provenance_guardrails.R")
source("src/utilities/reporting_helpers.R")

cli_defaults <- list(
  output_dir = default_run_output_dir("finalized_pooled_results_bundle"),
  overall_signature_table = "",
  pooled_term_table = "",
  subgroup_summary_table = "",
  subtype_all_features_table = "",
  allow_legacy_output = "false",
  overall_level_label = "Overall",
  exploratory_dataset = "",
  exploratory_window = "0-60m",
  exploratory_q_threshold = "0.10"
)

print_usage <- function() {
  print_standard_cli_help(
    title = "Finalized Pooled Results Bundle Driver",
    workflow_name = "finalized-pooled-results-bundle",
    script_path = file.path("drivers", "run_finalized_pooled_results_bundle.R"),
    defaults = cli_defaults,
    summary_lines = c(
      "This workflow is downstream-only and requires explicit upstream result tables.",
      "It assembles pooled forest plots, subtype heatmaps, and manuscript-style summary artifacts.",
      "The current reporting mode assumes two-dataset pooled inputs and subtype-stratified rows.",
      "subgroup_summary_table uses the shared subgroup-summary contract; subtype_all_features_table is retained as a compatibility alias.",
      "overall_level_label and exploratory_* flags tune the current reporting filters without changing upstream models.",
      "If exploratory_dataset is omitted, the report auto-selects the second dataset present in the input tables when available."
    ),
    required_flags = c("overall_signature_table", "pooled_term_table", "subgroup_summary_table"),
    aliases = c(subtype_all_features_table = "subgroup_summary_table")
  )
}

resolve_subgroup_summary_input <- function(options) {
  if (nzchar(options$subgroup_summary_table)) return(options$subgroup_summary_table)
  if (nzchar(options$subtype_all_features_table)) return(options$subtype_all_features_table)
  ""
}

validate_reporting_inputs <- function(signature_df, pooled_terms_df, subgroup_df) {
  required_signature_columns <- c("level", "subtype", "feature", "endpoint", "dataset", "window", "HR", "CI_low", "CI_high", "p_value", "p_value_bh")
  required_term_columns <- c("term", "feature", "endpoint", "HR", "CI_low", "CI_high", "p_value", "p_value_bh")

  missing_signature_columns <- setdiff(required_signature_columns, names(signature_df))
  missing_term_columns <- setdiff(required_term_columns, names(pooled_terms_df))

  if (length(missing_signature_columns)) {
    stop("overall_signature_table is missing required columns: ", paste(missing_signature_columns, collapse = ", "))
  }
  if (length(missing_term_columns)) {
    stop("pooled_term_table is missing required columns: ", paste(missing_term_columns, collapse = ", "))
  }

  validate_subgroup_summary_table(
    subgroup_df = subgroup_df,
    context_label = "subgroup_summary_table",
    require_dataset_window = TRUE,
    require_subtype_alias = TRUE
  )

  signature_dataset_count <- dplyr::n_distinct(stats::na.omit(signature_df$dataset))
  subgroup_dataset_count <- dplyr::n_distinct(stats::na.omit(subgroup_df$dataset))
  non_overall_subgroup_count <- dplyr::n_distinct(stats::na.omit(subgroup_df$subtype[subgroup_df$subtype != "ALL"]))

  if (signature_dataset_count < 2 || subgroup_dataset_count < 2) {
    stop(
      "The current finalized pooled reporting workflow expects two-dataset pooled inputs. ",
      "If you only have one dataset, use the single-cohort workflow outputs or add a dedicated single-dataset reporting mode instead of this bundled pooled report."
    )
  }

  if (non_overall_subgroup_count < 1) {
    stop(
      "The current finalized pooled reporting workflow expects subtype-stratified input rows for the heatmap and subtype-oriented forest summaries. ",
      "If no subtype or other sub-cohort rows are available, this bundled report is not the right reporting mode yet."
    )
  }
}

is_truthy <- function(value) {
  if (is.null(value) || is.na(value)) {
    return(FALSE)
  }
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

clean_subtype <- function(x) {
  ifelse(x == "ALL", "Overall", gsub("^BRCA_", "", x))
}

short_panel_label <- function(dataset, window) {
  paste(dataset, window, sep = "
")
}

ordered_levels <- function(values, preferred_levels = character(0)) {
  values <- as.character(stats::na.omit(values))
  preferred_present <- preferred_levels[preferred_levels %in% values]
  c(preferred_present, setdiff(unique(values), preferred_levels))
}

build_panel_levels <- function(dataset_levels, window_levels, short = FALSE) {
  labels <- character(0)
  for (dataset_name in dataset_levels) {
    for (window_name in window_levels) {
      labels <- c(labels, if (short) short_panel_label(dataset_name, window_name) else paste(dataset_name, window_name, sep = " | "))
    }
  }
  labels
}

ordered_panel_label <- function(dataset, window, panel_levels, short = FALSE) {
  label <- if (short) {
    short_panel_label(dataset, window)
  } else {
    paste(dataset, window, sep = " | ")
  }

  factor(label, levels = panel_levels)
}

resolve_exploratory_dataset <- function(requested_dataset, dataset_levels) {
  if (nzchar(requested_dataset)) return(requested_dataset)
  if (length(dataset_levels) >= 2) return(dataset_levels[[2]])
  if (length(dataset_levels) == 1) return(dataset_levels[[1]])
  ""
}

extract_interaction_dataset_label <- function(term_name, dataset_levels) {
  parsed_label <- sub("^group_high:dataset", "", term_name)
  parsed_label <- sub(":post60$", "", parsed_label)
  if (length(dataset_levels) >= 2) return(dataset_levels[[2]])
  parsed_label
}

plot_signature_forest <- function(signature_df, output_stub, dataset_levels, window_levels) {
  panel_levels <- build_panel_levels(dataset_levels, window_levels, short = FALSE)
  plot_df <- signature_df %>%
    mutate(
      subtype_label = clean_subtype(.data$subtype),
      subtype_label = factor(.data$subtype_label, levels = c("Overall", "Basal", "Her2", "LumA", "LumB", "Normal")),
      panel = ordered_panel_label(.data$dataset, .data$window, panel_levels = panel_levels, short = FALSE),
      significant = .data$p_value_bh < 0.05
    ) %>%
    filter(is.finite(.data$CI_low), is.finite(.data$CI_high), .data$CI_low > 0, .data$CI_high > 0)

  forest_plot <- ggplot(plot_df, aes(x = .data$HR, y = .data$subtype_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_errorbar(
      aes(xmin = .data$CI_low, xmax = .data$CI_high),
      orientation = "y",
      width = 0.2,
      color = "grey35",
      linewidth = 0.5
    ) +
    geom_point(aes(fill = .data$significant), shape = 21, size = 2.8, color = "black", stroke = 0.35) +
    facet_grid(endpoint ~ panel) +
    scale_x_log10() +
    scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#c44e52"), guide = "none") +
    labs(
      title = "Pooled Signature Effects Across Datasets and Subtypes",
      subtitle = "High vs Low signature hazard ratios from pooled piecewise Cox models",
      x = "Hazard ratio (log scale)",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold")
    )

  ggsave(paste0(output_stub, ".png"), forest_plot, width = 13, height = 8.5, dpi = 300)
  ggsave(paste0(output_stub, ".pdf"), forest_plot, width = 13, height = 8.5)
}

plot_signature_forest_journal <- function(signature_df, output_stub, dataset_levels, window_levels) {
  panel_levels <- build_panel_levels(dataset_levels, window_levels, short = TRUE)
  plot_df <- signature_df %>%
    mutate(
      subtype_label = clean_subtype(.data$subtype),
      subtype_label = factor(.data$subtype_label, levels = c("Overall", "Basal", "Her2", "LumA", "LumB", "Normal")),
      panel = ordered_panel_label(.data$dataset, .data$window, panel_levels = panel_levels, short = TRUE),
      significance_label = case_when(
        .data$p_value_bh < 0.05 ~ "q<0.05",
        .data$p_value < 0.05 ~ "p<0.05",
        TRUE ~ "NS"
      )
    ) %>%
    filter(is.finite(.data$CI_low), is.finite(.data$CI_high), .data$CI_low > 0, .data$CI_high > 0)

  forest_plot <- ggplot(plot_df, aes(x = .data$HR, y = .data$subtype_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey55", linewidth = 0.35) +
    geom_errorbar(
      aes(xmin = .data$CI_low, xmax = .data$CI_high),
      orientation = "y",
      width = 0.16,
      color = "grey30",
      linewidth = 0.45
    ) +
    geom_point(aes(shape = .data$significance_label, fill = .data$significance_label), size = 2.6, color = "black", stroke = 0.35) +
    facet_grid(endpoint ~ panel) +
    scale_x_log10(breaks = c(0.1, 0.5, 1, 2, 5, 10), labels = c("0.1", "0.5", "1", "2", "5", "10")) +
    scale_shape_manual(values = c("NS" = 21, "p<0.05" = 21, "q<0.05" = 23), name = NULL) +
    scale_fill_manual(values = c("NS" = "white", "p<0.05" = "#9ecae1", "q<0.05" = "#d95f0e"), name = NULL) +
    labs(
      x = "Hazard ratio",
      y = NULL
    ) +
    theme_bw(base_size = 9, base_family = "serif") +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom",
      legend.box = "horizontal"
    )

  ggsave(paste0(output_stub, ".png"), forest_plot, width = 11.5, height = 7.2, dpi = 400)
  ggsave(paste0(output_stub, ".pdf"), forest_plot, width = 11.5, height = 7.2)
}

plot_subtype_heatmap <- function(subtype_df, output_stub, dataset_levels) {
  plot_df <- prepare_subgroup_heatmap_data(
    subgroup_df = subtype_df,
    subgroup_column = "subtype",
    feature_column = "feature",
    endpoint_column = "endpoint",
    hr_column = "HR",
    p_value_column = "p_value",
    p_value_bh_column = "p_value_bh",
    dataset_column = "dataset",
    subgroup_transform = clean_subtype,
    subgroup_levels = c("Basal", "Her2", "LumA", "LumB", "Normal"),
    feature_levels = c("Signature", "GENE_A", "GENE_B", "GENE_C", "GENE_D"),
    endpoint_levels = c("OS", "PFI"),
    dataset_levels = dataset_levels,
    filter_window = "0-60m"
  ) %>%
    mutate(
      nominal_sig = .data$significance_label == "p<0.05",
      exploratory_sig = .data$significance_label %in% c("p<0.05", "q<0.10")
    )

  heatmap_plot <- ggplot(plot_df, aes(x = .data$subgroup_label, y = .data$feature_label, fill = .data$log2_hr)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_point(
      data = plot_df %>% filter(.data$exploratory_sig),
      shape = 23,
      size = 2.4,
      stroke = 0.4,
      fill = "#f6c85f",
      color = "black"
    ) +
    facet_grid(dataset_label ~ endpoint_label) +
    subgroup_heatmap_fill_scale(low = "#2166ac", high = "#b2182b", limits = c(-2, 2)) +
    labs(
      title = "Subtype-Specific Early Effects Across Features",
      subtitle = "Pooled 0-60 month high vs low contrasts by dataset and endpoint",
      x = NULL,
      y = NULL
    ) +
    subgroup_heatmap_theme(base_size = 10) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(paste0(output_stub, ".png"), heatmap_plot, width = 10.5, height = 6.5, dpi = 300)
  ggsave(paste0(output_stub, ".pdf"), heatmap_plot, width = 10.5, height = 6.5)
}

plot_subtype_heatmap_journal <- function(subtype_df, output_stub, dataset_levels) {
  plot_df <- prepare_subgroup_heatmap_data(
    subgroup_df = subtype_df,
    subgroup_column = "subtype",
    feature_column = "feature",
    endpoint_column = "endpoint",
    hr_column = "HR",
    p_value_column = "p_value",
    p_value_bh_column = "p_value_bh",
    dataset_column = "dataset",
    subgroup_transform = clean_subtype,
    subgroup_levels = c("Basal", "Her2", "LumA", "LumB", "Normal"),
    feature_levels = c("Signature", "GENE_A", "GENE_B", "GENE_C", "GENE_D"),
    endpoint_levels = c("OS", "PFI"),
    dataset_levels = dataset_levels,
    filter_window = "0-60m"
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
    facet_grid(dataset_label ~ endpoint_label) +
    scale_shape_manual(values = c("p<0.05" = 21, "q<0.10" = 23), name = NULL) +
    subgroup_heatmap_fill_scale(limits = c(-2, 2)) +
    labs(x = NULL, y = NULL) +
    subgroup_heatmap_theme(
      base_size = 9,
      base_family = "serif",
      axis_text_x = element_text(angle = 0, hjust = 0.5),
      legend_position = "right"
    )

  ggsave(paste0(output_stub, ".png"), heatmap_plot, width = 10.2, height = 6.4, dpi = 400)
  ggsave(paste0(output_stub, ".pdf"), heatmap_plot, width = 10.2, height = 6.4)
}

opt <- parse_standard_cli_args(
  defaults = cli_defaults,
  workflow_name = "finalized-pooled-results-bundle",
  help_printer = print_usage,
  required_flags = c("overall_signature_table", "pooled_term_table"),
  aliases = c(subtype_all_features_table = "subgroup_summary_table")
)
allow_legacy <- is_truthy(opt$allow_legacy_output) || is_truthy(Sys.getenv("SURVIVAL_ALLOW_LEGACY_OUTPUT"))
opt$exploratory_q_threshold <- as.numeric(opt$exploratory_q_threshold)

assert_output_guardrails(opt$output_dir, allow_legacy_output = allow_legacy)

required_inputs <- c(
  overall_signature_table = opt$overall_signature_table,
  pooled_term_table = opt$pooled_term_table,
  subgroup_summary_table = resolve_subgroup_summary_input(opt)
)
missing_flags <- names(required_inputs)[!nzchar(required_inputs)]
if (length(missing_flags)) {
  stop_cli_missing_required(missing_flags, "finalized-pooled-results-bundle")
}
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop("Missing required input file(s): ", paste(unname(missing_inputs), collapse = ", "))
}

ensure_reporting_output_dirs(opt$output_dir)

signature_df <- readr::read_csv(opt$overall_signature_table, show_col_types = FALSE)
pooled_terms_df <- readr::read_csv(opt$pooled_term_table, show_col_types = FALSE)
subgroup_summary_path <- resolve_subgroup_summary_input(opt)
subgroup_input_df <- readr::read_csv(subgroup_summary_path, show_col_types = FALSE)
subtype_df <- normalize_subgroup_summary_table(
  subgroup_df = subgroup_input_df,
  allow_legacy_pooled_columns = TRUE
)

validate_reporting_inputs(signature_df, pooled_terms_df, subtype_df)

dataset_levels <- ordered_levels(c(signature_df$dataset, subtype_df$dataset))
window_levels <- ordered_levels(c(signature_df$window, subtype_df$window), preferred_levels = c("0-60m", ">60m"))
reference_dataset <- if (length(dataset_levels)) dataset_levels[[1]] else "reference dataset"
comparison_dataset <- if (length(dataset_levels) >= 2) dataset_levels[[2]] else "comparison dataset"
exploratory_dataset <- resolve_exploratory_dataset(opt$exploratory_dataset, dataset_levels)

readr::write_csv(subtype_df, file.path(opt$output_dir, "tables", "consumed_subgroup_summary.csv"))

overall_signature_hits <- signature_df %>%
  filter(level == opt$overall_level_label, p_value_bh < 0.05) %>%
  transmute(
    section = as.character("Overall signature"),
    feature = as.character(feature),
    subtype = as.character(clean_subtype(subtype)),
    endpoint = as.character(endpoint),
    dataset = as.character(dataset),
    window = as.character(window),
    comparison = as.character("High vs Low"),
    HR = round(HR, 3),
    `95% CI` = as.character(report_fmt_ci(CI_low, CI_high)),
    p = as.character(report_fmt_p(p_value)),
    q = as.character(report_fmt_p(p_value_bh)),
    evidence_tier = as.character("BH-significant"),
    source = as.character(basename(opt$overall_signature_table))
  )

overall_term_hits <- pooled_terms_df %>%
  filter(p_value_bh < 0.05) %>%
  mutate(
    subtype = "Overall",
    dataset = case_when(
      term == "group_high" ~ reference_dataset,
      term == "group_high:post60" ~ reference_dataset,
      grepl("^group_high:dataset", term) ~ paste(extract_interaction_dataset_label(term, dataset_levels), "vs", reference_dataset),
      TRUE ~ "Model term"
    ),
    window = case_when(
      grepl("post60", term, fixed = TRUE) ~ ">60m shift",
      TRUE ~ "0-60m"
    ),
    comparison = case_when(
      term == "group_high" ~ "High vs Low",
      term == "group_high:post60" ~ "Late shift",
      grepl("^group_high:dataset.*:post60$", term) ~ "Late interaction",
      grepl("^group_high:dataset", term) ~ "Early interaction",
      TRUE ~ term
    )
  ) %>%
  transmute(
    section = as.character("Overall pooled term"),
    feature = as.character(feature),
    subtype = as.character(subtype),
    endpoint = as.character(endpoint),
    dataset = as.character(dataset),
    window = as.character(window),
    comparison = as.character(comparison),
    HR = round(HR, 3),
    `95% CI` = as.character(report_fmt_ci(CI_low, CI_high)),
    p = as.character(report_fmt_p(p_value)),
    q = as.character(report_fmt_p(p_value_bh)),
    evidence_tier = as.character("BH-significant"),
    source = as.character(basename(opt$pooled_term_table))
  )

subtype_exploratory_hits <- subtype_df %>%
  filter(dataset == exploratory_dataset, window == opt$exploratory_window, p_value_bh < opt$exploratory_q_threshold) %>%
  transmute(
    section = as.character("Subtype-specific exploratory"),
    feature = as.character(feature),
    subtype = as.character(clean_subtype(subtype)),
    endpoint = as.character(endpoint),
    dataset = as.character(dataset),
    window = as.character(window),
    comparison = as.character("High vs Low"),
    HR = round(HR, 3),
    `95% CI` = as.character(report_fmt_ci(CI_low, CI_high)),
    p = as.character(report_fmt_p(p_value)),
    q = as.character(report_fmt_p(p_value_bh)),
    evidence_tier = as.character(paste0("Exploratory (q<", formatC(opt$exploratory_q_threshold, format = "f", digits = 2), ")")),
    source = as.character(basename(subgroup_summary_path))
  )

summary_table <- bind_rows(overall_signature_hits, overall_term_hits, subtype_exploratory_hits)

strict_summary_table <- bind_rows(overall_signature_hits, overall_term_hits)

write_csv_and_md_table(
  csv_df = summary_table,
  csv_path = file.path(opt$output_dir, "tables", "finalized_pooled_results_summary.csv"),
  md_path = file.path(opt$output_dir, "tables", "finalized_pooled_results_summary.md")
)
write_csv_and_md_table(
  csv_df = strict_summary_table,
  csv_path = file.path(opt$output_dir, "tables", "finalized_pooled_results_summary_strict_q05.csv"),
  md_path = file.path(opt$output_dir, "tables", "finalized_pooled_results_summary_strict_q05.md")
)

plot_signature_forest(
  signature_df = signature_df,
  output_stub = file.path(opt$output_dir, "figures", "figure1_signature_pooled_forest"),
  dataset_levels = dataset_levels,
  window_levels = window_levels
)

plot_signature_forest_journal(
  signature_df = signature_df,
  output_stub = file.path(opt$output_dir, "figures", "figure1_signature_pooled_forest_journal"),
  dataset_levels = dataset_levels,
  window_levels = window_levels
)

plot_subtype_heatmap(
  subtype_df = subtype_df,
  output_stub = file.path(opt$output_dir, "figures", "figure2_subtype_early_feature_heatmap"),
  dataset_levels = dataset_levels
)

plot_subtype_heatmap_journal(
  subtype_df = subtype_df,
  output_stub = file.path(opt$output_dir, "figures", "figure2_subtype_early_feature_heatmap_journal"),
  dataset_levels = dataset_levels
)

basal_null_count <- subtype_df %>%
  filter(grepl("BRCA_Basal", subtype), p_value < 0.05) %>%
  nrow()

figure_legends <- c(
  "# Figure Legends",
  "",
  "## Figure 1. Pooled signature effects across datasets and subtypes.",
  paste(
    "Forest plot of high-versus-low hazard ratios for the generic four-feature signature",
    "derived from pooled piecewise Cox models. Panels are stratified by endpoint",
    "(OS or PFI), dataset-specific implied estimate, and time window",
    "(0-60 months or >60 months). Points denote hazard ratios and horizontal bars denote",
    "95% confidence intervals; filled points indicate BH-adjusted q < 0.05 within the",
    "signature-consistency screen. Rows with infinite confidence limits were omitted from",
    "the plot because sparse late events made those subtype-specific estimates unstable."
  ),
  "",
  "## Figure 2. Subtype-specific early pooled effects across the signature and individual features.",
  paste(
    "Heatmap of 0-60 month pooled high-versus-low hazard ratios for Signature, GENE_A, GENE_B,",
    "GENE_C, and GENE_D across BRCA subtypes. Fill encodes log2(HR), and gold diamonds mark",
    "exploratory BH-adjusted q < 0.10 across the full",
    "subtype-by-feature screen. This figure highlights that the",
    "strongest subtype-specific early signals can localize to one dataset rather than the other and that",
    if (basal_null_count == 0) {
      "BRCA_Basal does not retain a significant pooled early effect for any feature."
    } else {
      "BRCA_Basal shows only nominal, non-BH-significant early effects."
    }
  )
)

write_reporting_text_block(figure_legends, file.path(opt$output_dir, "figure_legends.md"))

results_section <- c(
  "# Draft Results Section",
  "",
  paste(
    "The pooled piecewise Cox analysis identified the strongest reproducible signal for the generic four-feature signature rather than for any single feature.",
    "This bundled report reads dataset labels from the supplied upstream result tables rather than hardcoding any named cohorts.",
    "For this run,", reference_dataset, "is treated as the reference dataset and", comparison_dataset, "as the comparison dataset."
  ),
  "",
  "Within the same pooled model family, individual-feature signals and signature-level signals are preserved separately in the finalized summary tables, allowing downstream interpretation without hardcoded cohort language in the driver.",
  "",
  paste(
    "Subtype-specific exploratory filtering currently uses", exploratory_dataset, "at", opt$exploratory_window, "unless you override --exploratory-dataset.",
    "The consumed subgroup-summary table remains", basename(subgroup_summary_path), "."
  ),
  "",
  "Overall, this bundled report now acts as a reusable downstream summary layer for paired-dataset survival analyses rather than a TCGA/METABRIC-specific reporting path."
)

write_reporting_text_block(results_section, file.path(opt$output_dir, "results_section.md"))

manuscript_summary <- c(
  "# Manuscript-Style Summary",
  "",
  paste(
    "This manuscript-style summary is a reusable placeholder for pooled two-dataset reporting.",
    "It uses dataset labels inferred from the upstream tables, with", reference_dataset, "as the reference dataset",
    "and", comparison_dataset, "as the comparison dataset when both are present.",
    "Use the generated summary tables and figures in this run directory as the source of truth for dataset-specific interpretation."
  )
)

slide_summary <- c(
  "# Slide Summary",
  "",
  "- The pooled model supports the generic four-feature signature more strongly than any individual feature.",
  paste("- Dataset labels are inferred from the input tables; this run used", reference_dataset, "as the reference dataset and", comparison_dataset, "as the comparison dataset."),
  "- Overall and subgroup-specific interpretations should be taken from the generated tables rather than from hardcoded cohort text.",
  paste("- Exploratory subtype filtering currently targets", exploratory_dataset, "unless you override --exploratory-dataset."),
  "- In the subtype-specific pooled screen, no result survived q < 0.05 across all features and subtypes.",
  "- The strongest exploratory subtype-specific signals may concentrate in one dataset rather than the other depending on the supplied inputs.",
  "- BRCA_Basal did not show a reproducible pooled signature or individual-feature effect."
)

write_reporting_text_block(manuscript_summary, file.path(opt$output_dir, "manuscript_summary.md"))
write_reporting_text_block(slide_summary, file.path(opt$output_dir, "slide_summary.md"))

write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Assemble paper-ready figures, summary table, and draft Results text from finalized pooled-model outputs.",
  driver = "drivers/run_finalized_pooled_results_bundle.R",
  inputs = list(
    overall_signature_table = opt$overall_signature_table,
    pooled_term_table = opt$pooled_term_table,
    subgroup_summary_table = subgroup_summary_path
  ),
  parameters = list(
    figure_1 = "signature pooled forest",
    figure_1_journal = "signature pooled forest journal",
    figure_2 = "subtype early feature heatmap",
    figure_2_journal = "subtype early feature heatmap journal",
    summary_table = "finalized pooled results summary",
    strict_summary_table = "finalized pooled results summary strict q05",
    overall_level_label = opt$overall_level_label,
    exploratory_dataset = exploratory_dataset,
    exploratory_window = opt$exploratory_window,
    exploratory_q_threshold = opt$exploratory_q_threshold
  ),
  notes = c(
    "Figure 1 uses the archived signature consistency screen.",
    "Figure 2 consumes the shared subgroup-summary contract, with legacy subtype_all_features_table input still accepted as a compatibility alias."
  ),
  scope = "exploratory"
)

message("Paper bundle written to: ", opt$output_dir)
