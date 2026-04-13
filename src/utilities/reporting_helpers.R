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

report_fmt_num <- function(x, digits = 3) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "f", digits = digits))
}

report_fmt_p <- function(x) {
  ifelse(
    is.na(x),
    NA_character_,
    ifelse(x < 0.001, formatC(x, format = "e", digits = 2), formatC(x, format = "f", digits = 3))
  )
}

report_fmt_ci <- function(low, high, digits = 3) {
  paste0("[", report_fmt_num(low, digits), ", ", report_fmt_num(high, digits), "]")
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

prepare_subgroup_forest_data <- function(
    subgroup_df,
    subgroup_column,
    feature_column,
    endpoint_column,
    comparison_column,
    hr_column,
    ci_low_column,
    ci_high_column,
    p_value_column,
    p_value_bh_column,
    subgroup_transform = identity,
    subgroup_levels = NULL,
    feature_levels = NULL,
    endpoint_levels = NULL,
    comparison_levels = NULL,
    q_threshold = 0.10,
    p_threshold = 0.05,
    q_label = "q<0.10",
    p_label = "p<0.05") {
  plot_df <- subgroup_df %>%
    dplyr::mutate(
      subgroup_label = subgroup_transform(as.character(.data[[subgroup_column]])),
      feature_label = as.character(.data[[feature_column]]),
      endpoint_label = as.character(.data[[endpoint_column]]),
      comparison_label = if (!is.null(comparison_column) && comparison_column %in% names(subgroup_df)) {
        as.character(.data[[comparison_column]])
      } else {
        "Overall"
      },
      HR = .data[[hr_column]],
      CI_low = .data[[ci_low_column]],
      CI_high = .data[[ci_high_column]],
      p_value = .data[[p_value_column]],
      p_value_bh = .data[[p_value_bh_column]],
      significance_label = dplyr::case_when(
        !is.na(.data$p_value_bh) & .data$p_value_bh < q_threshold ~ q_label,
        !is.na(.data$p_value) & .data$p_value < p_threshold ~ p_label,
        TRUE ~ "NS"
      )
    ) %>%
    dplyr::filter(
      is.finite(.data$HR),
      is.finite(.data$CI_low),
      is.finite(.data$CI_high),
      .data$HR > 0,
      .data$CI_low > 0,
      .data$CI_high > 0
    )

  if (!is.null(subgroup_levels)) {
    plot_df <- plot_df %>% dplyr::mutate(subgroup_label = factor(.data$subgroup_label, levels = subgroup_levels))
  }
  if (!is.null(feature_levels)) {
    plot_df <- plot_df %>% dplyr::mutate(feature_label = factor(.data$feature_label, levels = feature_levels))
  }
  if (!is.null(endpoint_levels)) {
    plot_df <- plot_df %>% dplyr::mutate(endpoint_label = factor(.data$endpoint_label, levels = endpoint_levels))
  }
  if (!is.null(comparison_levels)) {
    plot_df <- plot_df %>% dplyr::mutate(comparison_label = factor(.data$comparison_label, levels = comparison_levels))
  }

  plot_df
}

build_subgroup_forest <- function(
    subgroup_df,
    output_stub,
    subgroup_column,
    feature_column,
    endpoint_column,
    comparison_column,
    hr_column,
    ci_low_column,
    ci_high_column,
    p_value_column,
    p_value_bh_column,
    subgroup_transform = identity,
    subgroup_levels = NULL,
    feature_levels = NULL,
    endpoint_levels = NULL,
    comparison_levels = NULL,
    q_threshold = 0.10,
    p_threshold = 0.05,
    title = "Subgroup Effects",
    subtitle = "Hazard ratios across subgroup-defined analyses",
    width = 12,
    height = 8) {
  plot_df <- prepare_subgroup_forest_data(
    subgroup_df = subgroup_df,
    subgroup_column = subgroup_column,
    feature_column = feature_column,
    endpoint_column = endpoint_column,
    comparison_column = comparison_column,
    hr_column = hr_column,
    ci_low_column = ci_low_column,
    ci_high_column = ci_high_column,
    p_value_column = p_value_column,
    p_value_bh_column = p_value_bh_column,
    subgroup_transform = subgroup_transform,
    subgroup_levels = subgroup_levels,
    feature_levels = feature_levels,
    endpoint_levels = endpoint_levels,
    comparison_levels = comparison_levels,
    q_threshold = q_threshold,
    p_threshold = p_threshold
  )

  forest_plot <- ggplot(plot_df, aes(x = .data$HR, y = .data$subgroup_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_errorbar(
      aes(xmin = .data$CI_low, xmax = .data$CI_high),
      orientation = "y",
      width = 0.2,
      color = "grey35",
      linewidth = 0.5
    ) +
    geom_point(
      aes(fill = .data$significance_label),
      shape = 21,
      size = 2.8,
      color = "black",
      stroke = 0.35
    ) +
    facet_grid(endpoint_label ~ comparison_label + feature_label, scales = "free_y", space = "free_y") +
    scale_x_log10() +
    scale_fill_manual(
      values = c("NS" = "white", "p<0.05" = "#9ecae1", "q<0.10" = "#d95f0e"),
      name = NULL
    ) +
    labs(
      title = title,
      subtitle = subtitle,
      x = "Hazard ratio (log scale)",
      y = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_rect(fill = "grey95"),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  ggsave(paste0(output_stub, ".png"), forest_plot, width = width, height = height, dpi = 300)
  ggsave(paste0(output_stub, ".pdf"), forest_plot, width = width, height = height)

  invisible(plot_df)
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