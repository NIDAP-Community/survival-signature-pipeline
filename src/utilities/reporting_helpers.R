suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
})

ensure_reporting_output_dirs <- function(output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)
}

write_markdown_table <- function(df, file_path) {
  header <- paste0("| ", paste(names(df), collapse = " | "), " |")
  separator <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(row) paste0("| ", paste(row, collapse = " | "), " |"))
  writeLines(c(header, separator, rows), con = file_path)
}

write_csv_and_md_table <- function(csv_df, csv_path, md_path, md_df = csv_df) {
  readr::write_csv(csv_df, csv_path)
  write_markdown_table(md_df, md_path)
}

write_reporting_text_block <- function(lines, file_path) {
  writeLines(lines, con = file_path)
}

prepare_subgroup_heatmap_data <- function(
    subgroup_df,
    subgroup_column,
    feature_column,
    endpoint_column,
    hr_column,
    p_value_column,
    p_value_bh_column,
    comparison_column = NULL,
    dataset_column = NULL,
    subgroup_transform = identity,
    subgroup_levels = NULL,
    feature_levels = NULL,
    endpoint_levels = NULL,
    dataset_levels = NULL,
    filter_window = NULL,
    q_threshold = 0.10,
    p_threshold = 0.05,
    q_label = "q<0.10",
    p_label = "p<0.05") {
  plot_df <- subgroup_df

  if (!is.null(filter_window) && "window" %in% names(plot_df)) {
    plot_df <- plot_df %>% dplyr::filter(.data$window == filter_window)
  }

  plot_df <- plot_df %>%
    dplyr::mutate(
      subgroup_label = subgroup_transform(as.character(.data[[subgroup_column]])),
      feature_label = as.character(.data[[feature_column]]),
      endpoint_label = as.character(.data[[endpoint_column]]),
      log2_hr = log2(.data[[hr_column]]),
      significance_label = dplyr::case_when(
        !is.na(.data[[p_value_bh_column]]) & .data[[p_value_bh_column]] < q_threshold ~ q_label,
        !is.na(.data[[p_value_column]]) & .data[[p_value_column]] < p_threshold ~ p_label,
        TRUE ~ "NS"
      )
    )

  if (!is.null(comparison_column) && comparison_column %in% names(plot_df)) {
    plot_df <- plot_df %>% dplyr::mutate(comparison_label = as.character(.data[[comparison_column]]))
  }

  if (!is.null(dataset_column) && dataset_column %in% names(plot_df)) {
    plot_df <- plot_df %>% dplyr::mutate(dataset_label = as.character(.data[[dataset_column]]))
  }

  if (!is.null(subgroup_levels)) {
    plot_df <- plot_df %>% dplyr::mutate(subgroup_label = factor(.data$subgroup_label, levels = subgroup_levels))
  }
  if (!is.null(feature_levels)) {
    plot_df <- plot_df %>% dplyr::mutate(feature_label = factor(.data$feature_label, levels = feature_levels))
  }
  if (!is.null(endpoint_levels)) {
    plot_df <- plot_df %>% dplyr::mutate(endpoint_label = factor(.data$endpoint_label, levels = endpoint_levels))
  }
  if (!is.null(dataset_levels) && "dataset_label" %in% names(plot_df)) {
    plot_df <- plot_df %>% dplyr::mutate(dataset_label = factor(.data$dataset_label, levels = dataset_levels))
  }

  plot_df
}

subgroup_heatmap_fill_scale <- function(
    low = "#2b6cb0",
    mid = "#f7f7f7",
    high = "#b22222",
    midpoint = 0,
    limits = NULL,
    name = "log2(HR)") {
  ggplot2::scale_fill_gradient2(
    low = low,
    mid = mid,
    high = high,
    midpoint = midpoint,
    limits = limits,
    oob = scales::squish,
    name = name
  )
}

subgroup_heatmap_theme <- function(
    base_size = 10,
    base_family = NULL,
    axis_text_x = NULL,
    legend_position = NULL) {
  theme <- ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold")
    )

  if (!is.null(axis_text_x)) {
    theme <- theme + ggplot2::theme(axis.text.x = axis_text_x)
  }
  if (!is.null(legend_position)) {
    theme <- theme + ggplot2::theme(legend.position = legend_position)
  }

  theme
}