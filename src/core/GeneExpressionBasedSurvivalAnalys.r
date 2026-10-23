#' Gene Expression-Based Survival Analysis with Time-Dependent Effects
#'
#' Performs comprehensive survival analysis on gene expression data including
#' univariate and multivariate Cox regression, LASSO selection, Kaplan-Meier
#' plotting, and time-dependent proportional hazards analysis.
#'
#' @param expression_data Data frame or matrix with genes in rows and samples in columns.
#'   First column should contain gene identifiers.
#' @param metadata Data frame with sample metadata including survival outcomes and covariates.
#' @param sample_col Character. Column name in metadata containing sample identifiers.
#' @param gene_col Character. Column name in expression_data containing gene names.
#' @param genes Character vector. Gene names to analyze.
#' @param event_col Character. Column name for event indicator (1=event, 0=censored or "dead"/"alive").
#' @param time_col Character. Column name for survival/follow-up time.
#' @param clinical_covariates Character vector. Column names of clinical covariates to adjust for.
#' @param split_method Character. Method to split samples: "median", "quartiles", "tertiles", "percentiles", or "detectable".
#' @param percentile_cutoffs Numeric vector of length 2. Percentile cutoffs if split_method="percentiles".
#' @param detection_threshold Numeric. Expression threshold used when split_method="detectable".
#' @param include_groups Character vector. Which groups to include: "low", "mid", "high".
#' @param steps_to_run Character vector. Analysis steps to perform. Options:
#'   "individual+covariates", "individual", "signature+covariates", "signature",
#'   "lasso", "km_individual".
#' @param fdr_threshold Numeric. FDR threshold for significance (default 0.05).
#' @param run_lasso Logical. Whether to perform LASSO-COX selection.
#' @param run_ph_checks Logical. Whether to test proportional hazards assumption.
#' @param run_time_dependent Logical. Whether to fit time-dependent models if PH violated.
#' @param ph_alpha Numeric. Significance level for PH test (default 0.05).
#' @param cut_months Numeric. Time cutpoint (in months) for early/late analysis.
#' @param save_plots Logical. Whether to save individual KM plots as PNG files.
#' @param output_dir Character. Directory to save plots and results.
#' @param plot_params List. Plotting parameters (colors, labels, sizes, etc.).
#' @param export_csv Logical. Whether to export results as CSV files.
#' @param verbose Logical. Whether to print detailed progress messages.
#'
#' @return A list containing:
#'   \item{results_univariate}{Data frame with univariate Cox results per gene}
#'   \item{results_adjusted}{Data frame with covariate-adjusted Cox results per gene}
#'   \item{results_grouped_gene}{Data frame with grouped per-gene Cox results used by KM-style comparisons}
#'   \item{results_signature}{Cox results for gene signature}
#'   \item{lasso_selected}{Genes selected by LASSO}
#'   \item{ph_diagnostics}{Proportional hazards test results}
#'   \item{sample_groups}{Data frame mapping samples to expression groups}
#'   \item{plots}{List of ggsurvplot objects}
#'
#' @examples
#' \dontrun{
#' # Basic analysis
#' results <- analyze_survival_expression(
#'   expression_data = expr_matrix,
#'   metadata = clinical_data,
#'   sample_col = "patient_id",
#'   gene_col = "gene_name",
#'   genes = c("TP53", "BRCA1", "BRCA2"),
#'   event_col = "vital_status",
#'   time_col = "survival_months",
#'   clinical_covariates = c("age", "stage", "subtype")
#' )
#'
#' # Access results
#' print(results$results_adjusted)
#' print(results$plots$signature)
#' }
#'
#' @export
#' @importFrom dplyr filter mutate select rename case_when
#' @importFrom survival Surv coxph survfit survdiff cox.zph
#' @importFrom survminer ggsurvplot
#' @importFrom ggplot2 ggplot theme element_text scale_x_continuous ggsave
#' @importFrom glmnet cv.glmnet coef.glmnet
#' @importFrom stats median quantile var p.adjust complete.cases model.matrix
analyze_survival_expression <- function(
  expression_data,
  metadata,
  sample_col = "sample_id",
  gene_col = "gene",
  genes,
  event_col = "event",
  time_col = "time",
  clinical_covariates = NULL,
  split_method = c("quartiles", "tertiles", "median", "percentiles", "detectable"),
  percentile_cutoffs = c(25, 75),
  detection_threshold = 0,
  signature_aggregation = c("sum", "mean"),
  signature_continuous_cohort = c("grouped", "all"),
  include_groups = c("low", "high"),
  steps_to_run = c("individual+covariates", "individual", "signature+covariates", 
                   "signature", "lasso", "km_individual"),
  fdr_threshold = 0.05,
  p_adjust_method = "BH",
  run_lasso = FALSE,
  lasso_alpha = 1,
  lasso_nfolds = 10,
  lasso_lambda = "1se",
  run_ph_checks = TRUE,
  run_time_dependent = TRUE,
  ph_alpha = 0.05,
  cut_months = 50,
  save_plots = FALSE,
  output_dir = NULL,
  plot_params = list(),
  export_csv = FALSE,
  combine_gene_plots = TRUE,
  max_cols_combined = 3,
  max_plots_per_page = 9,
  gene_early_late = TRUE,
  verbose = TRUE
) {
  
  # ============================================================================
  # PARAMETER VALIDATION
  # ============================================================================
  
  split_method <- match.arg(split_method)
  signature_aggregation <- match.arg(signature_aggregation)
  signature_continuous_cohort <- match.arg(signature_continuous_cohort)
  
  if (!is.data.frame(expression_data) && !is.matrix(expression_data)) {
    stop("expression_data must be a data frame or matrix")
  }
  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }
  if (!sample_col %in% colnames(metadata)) {
    stop("sample_col '", sample_col, "' not found in metadata")
  }
  if (!gene_col %in% colnames(expression_data)) {
    stop("gene_col '", gene_col, "' not found in expression_data")
  }
  if (!event_col %in% colnames(metadata)) {
    stop("event_col '", event_col, "' not found in metadata")
  }
  if (!time_col %in% colnames(metadata)) {
    stop("time_col '", time_col, "' not found in metadata")
  }
  if (length(genes) == 0) {
    stop("genes vector cannot be empty")
  }
  if (save_plots && is.null(output_dir)) {
    stop("output_dir must be specified when save_plots=TRUE")
  }
  
  # Create output directory if needed
  if (save_plots && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ============================================================================
  # SETUP PLOTTING PARAMETERS
  # ============================================================================
  
  default_plot_params <- list(
        # Colors & labels (kept compatible with your version)
        color_low = "#377eb8", # Vivid Blue
        color_mid = "#4daf4a", # Green
        color_high = "#e41a1c", # Deep Red
  label_low = "Low",
  label_mid = "Mid",
  label_high = "High",


        # Titles & text sizes
        plot_title = "Survival Analysis",
        x_label = "Time (months)",
        y_label = "Survival Probability",
        title_size = 15,
        axis_title_size = 15,
        axis_text_size = 12,
  legend_title = "Group",
  legend_title_size = 11,
  legend_text_size = 10,
  legend_row_spacing = 2.5,  # additional control for spacing between legend entries


        # Legend box (mimic original)
  legend_position = c(0.03, 0.22),
        legend_just = c(0, 0),
        legend_box_linewidth = 0.3,


        # Risk table, median line, censor plot toggles
        show_confidence_interval = FALSE,
        show_risk_table = FALSE,
        risk_table_fontsize = 5,
        show_censor_plot = FALSE, # (aka ncensor.plot)
        show_median_survival = TRUE,
        include_sample_size = TRUE,


        # Line & censor glyphs (match original)
        line_size = 0.9,
        censor_shape = 124,
        censor_size = 2.4,


        # Scale & axis behavior
        surv_scale_percent = TRUE, # original used percent scale
        break_time_by = 12, # pre-breaks before custom axis
        months_to_years_threshold = 60, # switch at >= 60 months
        target_labels_years = 7,
        target_labels_months = 6,


  # Layout/margins + optional x compression (reduced side margins to widen plot area)
  margin_top = 0.6, margin_right = 1.2, margin_bottom = 0.8, margin_left = 1.2, # cm
        compress_x = 1, # 1 = full width; >1 squeezes x-range


        # Footer (HR/p lines) spacing and text size
        footer_offset = 1,
    hr_text_size = 4,
  stats_panel_mode = "caption", # 'caption', 'side', or 'bottom_table'
  stats_panel_mode_signature = "bottom_table", # for signature plots place table beneath main plot
    stats_panel_title = "Statistics",
    stats_panel_font_family = NULL
    )
  
  plot_params <- modifyList(default_plot_params, plot_params)
  
  # ============================================================================
  # HELPER FUNCTIONS
  # ============================================================================
  # Null coalescing helper
  `%||%` <- function(a, b) { if (!is.null(a) && length(a) > 0 && !all(is.na(a))) a else b }

  # Nice step helper for axis tick calculation (moved to top-level scope)
  .nice_step <- function(range_units, target_labels = 7) {
    raw <- range_units / target_labels
    base <- 10 ^ floor(log10(raw))
    mult <- raw / base
    step <- if (mult <= 1) 1 else if (mult <= 2) 2 else if (mult <= 5) 5 else 10
    base * step
  }

  wrap_title <- function(txt, width = 42) {
    if (is.null(txt) || length(txt) == 0 || !nzchar(txt)) return(txt)
    paste(strwrap(txt, width = width), collapse = "\n")
  }
  
  # Parse time fields (handles character, factor, numeric with commas)
  parse_time_numeric <- function(x) {
    if (is.numeric(x)) return(x)
    x0 <- as.character(x)
    x0[x0 %in% c("", "NA", "NaN", "null", "_Not_Available_", "_Not_Applicable_",
                 "[Not Available]", "[Not Applicable]")] <- NA
    x0 <- trimws(gsub(",", "", x0))
    x0 <- gsub("[^0-9eE+\\-\\.]", "", x0)
    out <- suppressWarnings(as.numeric(x0))
    newly_na <- sum(is.na(out) & !is.na(x))
    if (newly_na > 0 && verbose) {
      message(sprintf("parse_time_numeric: %d values could not be parsed and became NA.", newly_na))
    }
    out
  }
  
  # Clean covariate columns
  clean_covariates <- function(data, covariates) {
    special_missing <- c("", "NA", "NaN", "null", "_Not_Applicable_", "_Not_Available_",
                        "[Not Available]", "[Not Applicable]")
    
    # Ensure columns exist
    for (cv in covariates) {
      if (!(cv %in% colnames(data))) data[[cv]] <- NA_character_
    }

    # Helper for nice ticks:
    .nice_step <- function(range_units, target_labels = 7) {
        raw <- range_units / target_labels
        base <- 10 ^ floor(log10(raw))
        mult <- raw / base
        step <- if (mult <= 1) 1 else if (mult <= 2) 2 else if (mult <= 5) 5 else 10
        base * step
    }
    
    # Clean
    data <- data %>%
      dplyr::mutate(across(everything(), ~ if (is.factor(.)) as.character(.) else .)) %>%
      dplyr::mutate(dplyr::across(all_of(covariates), ~ {
        x <- as.character(.)
        x <- trimws(x)
        x[x %in% special_missing] <- NA_character_
        x
      }))
    
    # Convert numeric columns back to numeric
    for (cv in covariates) {
      x <- data[[cv]]
      if (!all(is.na(x))) {
        x_char <- as.character(x)
        nonmiss <- x_char[!is.na(x_char)]
        if (length(nonmiss) > 0 && all(grepl("^[-+]?[0-9]+(\\.[0-9]+)?$", nonmiss))) {
          data[[cv]] <- as.numeric(x_char)
        } else {
          data[[cv]] <- as.character(x_char)
        }
      }
    }
    
    data
  }

  # Align expression to metadata samples
  get_expr_aligned <- function(gene, gene_mat, merged_base, available_samples) {
    if (!(gene %in% rownames(gene_mat))) {
      return(rep(NA_real_, nrow(merged_base)))
    }
    idx <- match(merged_base$sample_id, available_samples)
    out <- rep(NA_real_, nrow(merged_base))
    ok  <- !is.na(idx)
    if (any(ok)) {
      out[ok] <- as.numeric(gene_mat[gene, available_samples[idx[ok]], drop = TRUE])
    }
    out
  }
  
  cluster_display_settings <- function(split_method, local_plot_params) {
    if (split_method == "detectable") {
      return(list(
        color_map = c(A = local_plot_params$color_low, C = local_plot_params$color_high),
        label_map = c(A = "Undetectable", C = "Detectable")
      ))
    }

    if (split_method == "median") {
      return(list(
        color_map = c(A = local_plot_params$color_low, B = local_plot_params$color_high),
        label_map = c(A = local_plot_params$label_low, B = local_plot_params$label_high)
      ))
    }

    list(
      color_map = c(
        A = local_plot_params$color_low,
        B = local_plot_params$color_mid,
        C = local_plot_params$color_high
      ),
      label_map = c(
        A = local_plot_params$label_low,
        B = local_plot_params$label_mid,
        C = local_plot_params$label_high
      )
    )
  }

  # Make expression groups
  make_clusters <- function(x_num, split_method, percentile_cutoffs, include_groups, detection_threshold) {
    if (split_method == "median") {
      med <- stats::median(x_num, na.rm = TRUE)
      cl <- ifelse(x_num < med, "A", ifelse(x_num >= med, "B", NA))
      levels <- c("A", "B")
    } else if (split_method == "quartiles") {
      # Use cut() to avoid dropping samples at boundaries
      keep <- c()
      if ("low"  %in% include_groups) keep <- c(keep, "A")
      if ("mid"  %in% include_groups) keep <- c(keep, "B")
      if ("high" %in% include_groups) keep <- c(keep, "C")
      q25 <- stats::quantile(x_num, 0.25, na.rm = TRUE)
      q75 <- stats::quantile(x_num, 0.75, na.rm = TRUE)
      if (!is.finite(q25) || !is.finite(q75) || q25 >= q75) {
        return(factor(rep(NA_character_, length(x_num)), levels = keep))
      }
      breaks <- c(-Inf, q25, q75, Inf)
      cl <- cut(x_num, breaks = breaks, labels = c("A", "B", "C"), 
                include.lowest = TRUE, right = TRUE)
      cl <- as.character(cl)
      cl[!(cl %in% keep)] <- NA
      levels <- keep
    } else if (split_method == "tertiles") {
      # Use cut() to avoid dropping samples at boundaries
      keep <- c()
      if ("low"  %in% include_groups) keep <- c(keep, "A")
      if ("mid"  %in% include_groups) keep <- c(keep, "B")
      if ("high" %in% include_groups) keep <- c(keep, "C")
      t1 <- stats::quantile(x_num, 1/3, na.rm = TRUE)
      t2 <- stats::quantile(x_num, 2/3, na.rm = TRUE)
      if (!is.finite(t1) || !is.finite(t2) || t1 >= t2) {
        return(factor(rep(NA_character_, length(x_num)), levels = keep))
      }
      breaks <- c(-Inf, t1, t2, Inf)
      cl <- cut(x_num, breaks = breaks, labels = c("A", "B", "C"), 
                include.lowest = TRUE, right = TRUE)
      cl <- as.character(cl)
      cl[!(cl %in% keep)] <- NA
      levels <- keep
    } else if (split_method == "percentiles") {
      stopifnot(length(percentile_cutoffs) == 2)
      lo <- stats::quantile(x_num, percentile_cutoffs[1]/100, na.rm = TRUE)
      hi <- stats::quantile(x_num, percentile_cutoffs[2]/100, na.rm = TRUE)
      cl <- ifelse(x_num <= lo, "A",
                  ifelse(x_num >= hi, "C",
                        ifelse(x_num > lo & x_num < hi, "B", NA)))
      keep <- c()
      if ("low"  %in% include_groups) keep <- c(keep, "A")
      if ("mid"  %in% include_groups) keep <- c(keep, "B")
      if ("high" %in% include_groups) keep <- c(keep, "C")
      cl[!(cl %in% keep)] <- NA
      levels <- keep
    } else if (split_method == "detectable") {
      cl <- ifelse(is.na(x_num), NA, ifelse(x_num <= detection_threshold, "A", "C"))
      keep <- c()
      if ("low" %in% include_groups) keep <- c(keep, "A")
      if ("high" %in% include_groups) keep <- c(keep, "C")
      cl[!(cl %in% keep)] <- NA
      levels <- keep
    }
    factor(cl, levels = levels)
  }
  
  # Format time axis (months vs years)
  format_time_axis <- function(p, surv_df, plot_params) {
    max_t <- suppressWarnings(max(surv_df$time, na.rm = TRUE))
    if (!is.finite(max_t) || max_t <= 0) return(p)

    thr <- plot_params$months_to_years_threshold
    if (is.null(thr)) thr <- 60

    if (max_t >= thr) {
        # YEARS mode
        max_years_exact <- max_t / 12
        step_y <- .nice_step(max_years_exact, target_labels = plot_params$target_labels_years %||% 7)
        brks_y <- seq(0, floor(max_years_exact / step_y) * step_y, by = step_y)
        brks_m <- brks_y * 12


        # Replace *existing* x scale only once to suppress warning spam.
        p$plot$scales$scales <- Filter(function(s) {
          !inherits(s, "ScaleContinuousPosition") || !identical(s$aesthetics, "x")
        }, p$plot$scales$scales)
        p$plot <- p$plot + ggplot2::scale_x_continuous(
          limits = c(0, max_t), breaks = brks_m, labels = paste0(brks_y, "y"),
          expand = ggplot2::expansion(mult = c(0, 0.02)),
          guide = if (utils::packageVersion("ggplot2") >= "3.3.0") ggplot2::guide_axis(check.overlap = TRUE) else NULL
        ) + ggplot2::xlab("Time (years)")
    } else {
        # MONTHS mode
        brks_m <- pretty(c(0, max_t), n = plot_params$target_labels_months %||% 6)
        p$plot$scales$scales <- Filter(function(s) {
          !inherits(s, "ScaleContinuousPosition") || !identical(s$aesthetics, "x")
        }, p$plot$scales$scales)
        p$plot <- p$plot + ggplot2::scale_x_continuous(
          limits = c(0, max_t), breaks = brks_m, labels = as.integer(brks_m),
          expand = ggplot2::expansion(mult = c(0, 0.02)),
          guide = if (utils::packageVersion("ggplot2") >= "3.3.0") ggplot2::guide_axis(check.overlap = TRUE) else NULL
        ) + ggplot2::xlab("Time (months)")
    }
    p
}
  
  # Helper: suppress only benign duplicate scale/coord messages while preserving other messages
  muffle_scale_messages <- function(expr) {
    withCallingHandlers(expr, message = function(m) {
      msg <- conditionMessage(m)
      if (grepl("^Scale for x is already present", msg) ||
          grepl("^Adding another scale for x", msg) ||
          grepl("^Coordinate system already present", msg)) {
        invokeRestart("muffleMessage")
      }
    })
  }

  # Draw KM plot (now accepts local_plot_params for overrides)
  draw_km_plot <- function(fit, surv_df, cox_unadj = NULL, cox_adj = NULL,
                           cont_unadj = NULL, cont_adj = NULL,
                           main_title = "", local_plot_params = NULL,
                           attempted_adj = FALSE, adj_reason = NULL,
                           retained_covariates = character(0)) {
    if (is.null(local_plot_params)) local_plot_params <- plot_params
    stats_mode <- local_plot_params$stats_panel_mode %||% "caption"
    if (getOption("km_debug", FALSE)) message("[draw_km_plot] stats_mode initial=", stats_mode)
    strata_levels <- character(0)
    if (!is.null(fit$strata)) {
      strata_levels <- sub("^.*=", "", names(fit$strata))
    }
    if (!length(strata_levels)) {
      strata_levels <- unique(as.character(stats::na.omit(surv_df$cluster)))
    }
    recognized_levels <- c("A", "B", "C")
    if (length(strata_levels) && all(strata_levels %in% recognized_levels)) {
      display_settings <- cluster_display_settings(split_method, local_plot_params)
      color_map <- display_settings$color_map
      label_map <- display_settings$label_map
      colors <- unname(color_map[strata_levels])
      labels <- unname(label_map[strata_levels])
    } else {
      fallback_palette <- c(
        local_plot_params$color_low,
        local_plot_params$color_mid,
        local_plot_params$color_high
      )
      colors <- fallback_palette[seq_along(strata_levels)]
      labels <- strata_levels
    }

    if (isTRUE(local_plot_params$include_sample_size)) {
      group_sizes <- table(factor(as.character(surv_df$cluster), levels = strata_levels))[strata_levels]
      labels <- paste0(labels, " (n = ", as.integer(group_sizes), ")")
    }
    p <- muffle_scale_messages(survminer::ggsurvplot(
      fit, data = surv_df,
      size = local_plot_params$line_size %||% 0.9,
      censor.shape = local_plot_params$censor_shape %||% 124,
      censor.size = local_plot_params$censor_size %||% 2.4,
      conf.int = local_plot_params$show_confidence_interval,
      risk.table = local_plot_params$show_risk_table,
      risk.table.col = "strata",
      risk.table.fontsize = local_plot_params$risk_table_fontsize %||% 5,
      title = main_title,
      xlab = local_plot_params$x_label,
      ylab = local_plot_params$y_label,
      break.time.by = NULL,
      surv.median.line = if (isTRUE(local_plot_params$show_median_survival)) "v" else "none",
      ggtheme = ggplot2::theme_classic(),
      palette = colors,
      ncensor.plot = isTRUE(local_plot_params$show_censor_plot),
      legend.title = local_plot_params$legend_title,
      legend.labs = labels
    ))
    # Add percent scale only if not already present to avoid duplicate scale warnings
  if (isTRUE(local_plot_params$surv_scale_percent)) {
      existing_y <- any(vapply(p$plot$scales$scales, function(s) {
        "y" %in% tryCatch(s$aesthetics, error = function(...) character())
      }, logical(1)))
      if (!existing_y) {
        p$plot <- p$plot + ggplot2::scale_y_continuous(labels = scales::percent_format())
      }
    }
  if (isTRUE(local_plot_params$show_risk_table) && !is.null(p$table)) {
      p$table <- p$table + ggplot2::theme_classic()
    }
    # Capture user-defined margins first so we can extend bottom margin dynamically for caption annotations later.
    base_margin <- list(
      top = plot_params$margin_top,
      right = plot_params$margin_right,
      bottom = plot_params$margin_bottom,
      left = plot_params$margin_left
    )

    p$plot <- p$plot +
      ggplot2::theme(
        plot.title = ggplot2::element_text(size = local_plot_params$title_size, hjust = 0.5, face = "bold"),
        axis.title = ggplot2::element_text(size = local_plot_params$axis_title_size),
        axis.text = ggplot2::element_text(size = local_plot_params$axis_text_size),
        legend.position.inside = c(0.02, 0.02),
        legend.justification = c(0, 0),
        legend.key.height = grid::unit(6, "pt"),
        legend.margin = ggplot2::margin(3, 6, 3, 6),
        legend.box.margin = ggplot2::margin(3, 4, 3, 4),
        legend.key = ggplot2::element_blank(),
        legend.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.3),
        legend.title = ggplot2::element_text(size = local_plot_params$legend_title_size * 0.75),
        legend.text = ggplot2::element_text(size = local_plot_params$legend_text_size * 0.75, vjust = 0.5, margin = ggplot2::margin(r = 3)),
        plot.margin = ggplot2::margin(t = base_margin$top, r = base_margin$right, b = base_margin$bottom + 1, l = base_margin$left, unit = "cm")
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(ncol = 1, title.position = "top", label.position = "right",
                                                    byrow = TRUE, override.aes = list(linetype = 1, size = 0.6, shape = NA)))
    # Avoid re-adding a coord if one already present (ggsurvplot returns one). Only add if compression requested.
    comp <- plot_params$compress_x %||% 1
    if (is.numeric(comp) && comp > 1) {
      # Replace existing coord_cartesian by editing limits directly
      range_x <- c(min(surv_df$time, na.rm = TRUE), max(surv_df$time, na.rm = TRUE) / comp)
      p$plot$coordinates$limits$x <- range_x
    }
    # Format time axis (suppress duplicate scale warnings from ggplot2)
    # Extract existing breaks to decide if replacement needed
    existing_breaks <- NULL
    if (length(p$plot$scales$scales)) {
      sx <- Filter(function(s) inherits(s, 'ScaleContinuousPosition') && any('x' %in% s$aesthetics), p$plot$scales$scales)
      if (length(sx)) {
        existing_breaks <- tryCatch(sx[[1]]$get_breaks(), error = function(e) NULL)
      }
    }
  p2 <- muffle_scale_messages(suppressMessages(suppressWarnings(format_time_axis(p, surv_df, plot_params))))
    # Compare breaks (if both numeric and identical length + values skip reassignment to avoid message)
    new_scale <- FALSE
    if (!is.null(existing_breaks)) {
      new_scale <- TRUE
    }
    p <- p2

    primary_contrast_label <- if (split_method == "detectable") "Detectable vs Undetectable" else "High vs Low"
    secondary_contrast_label <- if (split_method %in% c("detectable", "median")) NULL else "Mid vs Low"

    # Helper to build comparison lines for the active grouping mode.
    build_lines <- function(cox_summary, prefix) {
      if (is.null(cox_summary)) return(character(0))
      rn <- rownames(cox_summary$coefficients)
      out <- character(0)
      add_line <- function(term, label) {
        idx <- match(term, rn)
        if (!is.na(idx)) {
          hr <- cox_summary$conf.int[idx, "exp(coef)"]
          lo <- cox_summary$conf.int[idx, "lower .95"]
          hi <- cox_summary$conf.int[idx, "upper .95"]
          pval <- cox_summary$coefficients[idx, "Pr(>|z|)"]
          out <<- c(out, sprintf("%s %s HR = %.2f (%.2f-%.2f), p = %s", prefix, label, hr, lo, hi, format.pval(pval, digits = 3)))
        }
      }
      if (split_method == "median" && "clusterB" %in% rn) add_line("clusterB", primary_contrast_label)
      if (split_method != "median" && "clusterC" %in% rn) add_line("clusterC", primary_contrast_label)
      if (!is.null(secondary_contrast_label) && "clusterB" %in% rn) add_line("clusterB", secondary_contrast_label)
      if (!length(out)) { # fallback any cluster term
        idx_any <- grep("^cluster", rn)
        if (length(idx_any)) {
          idx <- idx_any[1]
          hr <- cox_summary$conf.int[idx, "exp(coef)"]
          lo <- cox_summary$conf.int[idx, "lower .95"]
          hi <- cox_summary$conf.int[idx, "upper .95"]
          pval <- cox_summary$coefficients[idx, "Pr(>|z|)"]
          out <- c(out, sprintf("%s HR = %.2f (%.2f-%.2f), p = %s", prefix, hr, lo, hi, format.pval(pval, digits = 3)))
        }
      }
      out
    }
  lines_unadj <- build_lines(cox_unadj, "Cox (unadj)")
  lines_adj   <- build_lines(cox_adj,   "Cox (adj)")

    # Global log-rank p-value (retain numeric for table rather than re-parsing formatted string)
    global_line <- NULL
    global_p_numeric <- NA_real_
    tmp_glob <- tryCatch({
      sd <- survival::survdiff(survival::Surv(time, event) ~ cluster, data = surv_df)
      df <- length(sd$n) - 1
      if (df > 0) {
        p_glob <- stats::pchisq(sd$chisq, df = df, lower.tail = FALSE)
        global_p_numeric <<- p_glob
        global_line <<- sprintf("Log-rank (overall) p = %s", format.pval(p_glob, digits = 3))
      }
      TRUE
    }, error = function(e) NULL)

    # Build structured table rows (Test, HR, CI, p).
    build_table_rows <- function(prefix, cox_summary) {
      # If adjusted model was attempted but not available, emit placeholder rows so user sees absence explicitly
      if (is.null(cox_summary)) {
        if (identical(prefix, "adj") && isTRUE(attempted_adj)) {
          # Provide placeholders for both comparisons
          pl <- data.frame(
            Test = if (split_method == "detectable") "adj Detectable vs Undetectable" else if (split_method == "median") "adj H vs L" else c("adj H vs L", "adj M vs L"),
            HR = NA_real_,
            CI = "",
            p = NA_real_
          )
          if (!is.null(adj_reason)) attr(pl, "adj_reason") <- adj_reason
          return(pl)
        } else {
          return(data.frame())
        }
      }
      rn <- rownames(cox_summary$coefficients)
      rows <- list()
      add_row <- function(term, comp_label) {
        idx <- match(term, rn)
        if (!is.na(idx)) {
          hr  <- cox_summary$conf.int[idx, "exp(coef)"]
          lo  <- cox_summary$conf.int[idx, "lower .95"]
          hi  <- cox_summary$conf.int[idx, "upper .95"]
          pvl <- cox_summary$coefficients[idx, "Pr(>|z|)"]
          # Use shorter test names
          test_name <- paste0(gsub("Cox \\(", "", gsub("\\)", "", prefix)), " ", comp_label)
          rows[[length(rows)+1]] <<- data.frame(Test = test_name, HR = hr, CI = sprintf("%.2f-%.2f", lo, hi), p = pvl)
        }
      }
      if (split_method == "median" && "clusterB" %in% rn) add_row("clusterB", "H vs L")
      if (split_method != "median" && "clusterC" %in% rn) add_row("clusterC", if (split_method == "detectable") "Detectable vs Undetectable" else "H vs L")
      if (!split_method %in% c("detectable", "median") && "clusterB" %in% rn) add_row("clusterB", "M vs L")
      if (!length(rows)) {
        idx_any <- grep("^cluster", rn)
        if (length(idx_any)) {
          i <- idx_any[1]
          test_name <- gsub("Cox \\(", "", gsub("\\)", "", prefix))
          rows[[1]] <- data.frame(Test = test_name, HR = cox_summary$conf.int[i, "exp(coef)"], CI = sprintf("%.2f-%.2f", cox_summary$conf.int[i, "lower .95"], cox_summary$conf.int[i, "upper .95"]), p = cox_summary$coefficients[i, "Pr(>|z|)"])
        }
      }
      do.call(rbind, rows)
    }

    tbl_unadj <- build_table_rows("unadj", cox_unadj)
  tbl_adj   <- build_table_rows("adj", cox_adj)
    tbl_lr <- if (!is.null(global_line)) {
      data.frame(Test = "Log-rank (global)", HR = NA_real_, CI = "-", p = global_p_numeric)
    } else data.frame()
    # Continuous model rows (if supplied). They correspond to results_univariate / results_adjusted (continuous expression) for the gene.
    build_cont_row <- function(sm, label) {
      if (is.null(sm)) return(NULL)
      rn <- rownames(sm$coefficients)
      if (!length(rn)) return(NULL)
      idx <- 1
      hr  <- sm$conf.int[idx, "exp(coef)"]
      lo  <- sm$conf.int[idx, "lower .95"]
      hi  <- sm$conf.int[idx, "upper .95"]
      pvl <- sm$coefficients[idx, "Pr(>|z|)"]
      data.frame(Test = label, HR = hr, CI = sprintf("%.2f-%.2f", lo, hi), p = pvl)
    }
    cont_rows <- do.call(rbind, list(
      build_cont_row(cont_unadj, "cont (unadj)"),
      build_cont_row(cont_adj,   "cont (adj)")
    ))
    stats_tbl <- rbind(tbl_lr, cont_rows, tbl_unadj, tbl_adj)
    # Initialize sample size column
    if (nrow(stats_tbl)) {
      stats_tbl$N <- NA_integer_
      # Use number of rows in surv_df for global and continuous unadj if not provided by summary
      stats_tbl$N[stats_tbl$Test == "Log-rank (global)"] <- nrow(surv_df)
      if (!is.null(cont_unadj)) stats_tbl$N[stats_tbl$Test == "cont (unadj)"] <- cont_unadj$n
      if (!is.null(cont_adj))   stats_tbl$N[stats_tbl$Test == "cont (adj)"]   <- cont_adj$n
      # Pairwise grouped models: use n (total) from the corresponding Cox summaries
      if (!is.null(cox_unadj)) stats_tbl$N[grepl("^unadj ", stats_tbl$Test)] <- cox_unadj$n
      if (!is.null(cox_adj))   stats_tbl$N[grepl("^adj ",   stats_tbl$Test)] <- cox_adj$n
    }
    if (nrow(stats_tbl) && identical(stats_mode, "side")) {
      if (getOption("km_debug", FALSE)) message("[draw_km_plot] building side panel with ", nrow(stats_tbl), " rows")
      if (requireNamespace("patchwork", quietly = TRUE)) {
    # Format columns for display (side panel)
    stats_tbl$HR_disp <- ifelse(is.na(stats_tbl$HR), "-", sprintf("%.2f", stats_tbl$HR))
    stats_tbl$N_disp  <- ifelse(is.na(stats_tbl$N), "-", as.character(stats_tbl$N))
    stats_tbl$P_disp  <- ifelse(is.na(stats_tbl$p), "-", ifelse(stats_tbl$p < 1e-3, formatC(stats_tbl$p, format = "e", digits = 2), format.pval(stats_tbl$p, digits = 3)))
        stats_tbl$Row <- seq_len(nrow(stats_tbl))
        panel <- ggplot2::ggplot(stats_tbl, ggplot2::aes(y = rev(Row))) +
          # Column headers
          ggplot2::annotate("text", x = 0, y = max(stats_tbl$Row)+1.5, label = "Test", hjust = 0, fontface="bold", size = 3.8) +
          ggplot2::annotate("text", x = 1.3, y = max(stats_tbl$Row)+1.5, label = "HR", hjust = 0.5, fontface="bold", size = 3.8) +
          ggplot2::annotate("text", x = 2.2, y = max(stats_tbl$Row)+1.5, label = "95% CI", hjust = 0.5, fontface="bold", size = 3.8) +
          ggplot2::annotate("text", x = 3.2, y = max(stats_tbl$Row)+1.5, label = "N", hjust = 0.5, fontface="bold", size = 3.8) +
          ggplot2::annotate("text", x = 4.0, y = max(stats_tbl$Row)+1.5, label = "p", hjust = 0.5, fontface="bold", size = 3.8) +
          # Data rows
          ggplot2::geom_text(ggplot2::aes(x = 0, label = Test), hjust = 0, size = 3.5, fontface="plain") +
          ggplot2::geom_text(ggplot2::aes(x = 1.3, label = HR_disp), hjust = 0.5, size = 3.5) +
          ggplot2::geom_text(ggplot2::aes(x = 2.2, label = CI), hjust = 0.5, size = 3.5) +
          ggplot2::geom_text(ggplot2::aes(x = 3.2, label = N_disp), hjust = 0.5, size = 3.5) +
          ggplot2::geom_text(ggplot2::aes(x = 4.0, label = P_disp), hjust = 0.5, size = 3.5) +
          # Title
          ggplot2::annotate("text", x = 0, y = max(stats_tbl$Row)+2.5, label = local_plot_params$stats_panel_title %||% "Statistics", hjust = 0, fontface="bold", size = 4.2) +
          ggplot2::xlim(c(0, 4.4)) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult=c(0.02,0.12))) +
          ggplot2::theme_void() +
          ggplot2::theme(plot.margin = ggplot2::margin(t = base_margin$top, r = 2, b = base_margin$bottom, l = 4))
        
        # Ensure legend is positioned correctly before combining
        p$plot <- p$plot + ggplot2::theme(
          legend.position.inside = c(0.05, 0.05),
          legend.justification = c(0, 0),
          legend.background = ggplot2::element_rect(fill = "white", color = "black", linewidth = 0.5)
        )
        
        combined <- patchwork::wrap_plots(p$plot, panel, widths = c(4, 1.8))
        p$plot <- combined + ggplot2::labs(caption = NULL)
      } else {
        warning("patchwork not installed: falling back to caption mode for statistics panel.")
        stats_mode <- "caption"
      }
    }
    # New mode: bottom_table (place a compact table beneath the survival curve spanning full width)
    if (nrow(stats_tbl) && identical(stats_mode, "bottom_table")) {
      if (getOption("km_debug", FALSE)) message("[draw_km_plot] building bottom table with ", nrow(stats_tbl), " rows")
    stats_tbl$HR_disp <- ifelse(is.na(stats_tbl$HR), "-", sprintf("%.2f", stats_tbl$HR))
    stats_tbl$N_disp  <- ifelse(is.na(stats_tbl$N), "-", as.character(stats_tbl$N))
    stats_tbl$P_disp  <- ifelse(is.na(stats_tbl$p), "-", ifelse(stats_tbl$p < 1e-3, formatC(stats_tbl$p, format = "e", digits = 2), format.pval(stats_tbl$p, digits = 3)))
      stats_tbl$Row <- seq_len(nrow(stats_tbl))
      # Allow override of column x-positions via plot_params$stats_col_positions
      if (!is.null(local_plot_params$stats_col_positions)) {
        col_x <- local_plot_params$stats_col_positions
      } else {
        # Increased spacing between N and p columns to reduce overlap in multi-panel layouts
        col_x <- c(Test = 0.01, HR = 0.33, CI = 0.56, N = 0.74, p = 0.98)
      }
      n_rows <- nrow(stats_tbl)
      header_y <- n_rows + 0.9
      title_y  <- n_rows + 1.8
      # Post-process display fields (Option 2): cap extreme magnitudes instead of scientific notation.
      cap_large  <- function(x) {
        ifelse(is.na(x), NA_character_, ifelse(x >= 1e4, ">1e4", sprintf("%.2f", x)))
      }
      cap_small  <- function(x) {
        ifelse(is.na(x), NA_character_, ifelse(x > 0 & x < 1e-6, "<1e-6", sprintf("%.2g", x)))
      }
      # HR display: replace existing numeric with capped string (keep '-' for NA)
      stats_tbl$HR_disp <- ifelse(is.na(stats_tbl$HR), "-", cap_large(stats_tbl$HR))
      # Confidence interval: parse and cap; handle infinite or zero-expansion gracefully
      if ("CI" %in% names(stats_tbl)) {
        ci_raw <- stats_tbl$CI
        stats_tbl$CI <- vapply(seq_along(ci_raw), function(i) {
          ci_entry <- ci_raw[i]
            if (is.na(ci_entry) || !nzchar(ci_entry)) return("")
            parts <- strsplit(ci_entry, "-", fixed = TRUE)[[1]]
            if (length(parts) != 2) return(ci_entry)
            lo <- suppressWarnings(as.numeric(parts[1])); hi <- suppressWarnings(as.numeric(parts[2]))
            lo_lab <- if (is.infinite(lo) && lo < 0) "-Inf" else if (is.infinite(lo) && lo > 0) "+Inf" else if (!is.na(lo) && lo == 0) "0" else cap_large(lo)
            hi_lab <- if (is.infinite(hi) && hi < 0) "-Inf" else if (is.infinite(hi) && hi > 0) "+Inf" else if (!is.na(hi) && hi == 0) "0" else cap_large(hi)
            if (!is.na(lo) && !is.na(hi) && lo == 0 && is.infinite(hi)) return("0-Inf")
            paste0(lo_lab, "-", hi_lab)
        }, character(1))
      }
      # p-value display using caps
      stats_tbl$P_disp <- ifelse(is.na(stats_tbl$p), "-", cap_small(stats_tbl$p))
      # Adjust font size if any tokens exceed 8 chars
      font_size_data <- if (max(nchar(c(stats_tbl$HR_disp, stats_tbl$CI, stats_tbl$P_disp)), na.rm = TRUE) > 10) 3.0 else 3.3
      table_plot <- ggplot2::ggplot(stats_tbl, ggplot2::aes(y = rev(Row))) +
        ggplot2::annotate("text", x = col_x["Test"], y = title_y, label = local_plot_params$stats_panel_title %||% "Statistics", hjust = 0, fontface = "bold", size = 4) +
        ggplot2::annotate("text", x = col_x["Test"], y = header_y, label = "Test", hjust = 0, fontface = "bold", size = 3.5) +
        ggplot2::annotate("text", x = col_x["HR"],   y = header_y, label = "HR",   hjust = 0.5, fontface = "bold", size = 3.5) +
  ggplot2::annotate("text", x = col_x["CI"],   y = header_y, label = "95% CI", hjust = 0.5, fontface = "bold", size = 3.5) +
  ggplot2::annotate("text", x = col_x["N"],    y = header_y, label = "N", hjust = 0.5, fontface = "bold", size = 3.5) +
  ggplot2::annotate("text", x = col_x["p"],    y = header_y, label = "p",    hjust = 1,   fontface = "bold", size = 3.5) +
    ggplot2::geom_text(ggplot2::aes(x = col_x["Test"], label = Test), hjust = 0, size = font_size_data) +
    ggplot2::geom_text(ggplot2::aes(x = col_x["HR"],   label = HR_disp), hjust = 0.5, size = font_size_data) +
  ggplot2::geom_text(ggplot2::aes(x = col_x["CI"],   label = CI), hjust = 0.5, size = font_size_data) +
  ggplot2::geom_text(ggplot2::aes(x = col_x["N"],    label = N_disp), hjust = 0.5, size = font_size_data) +
  ggplot2::geom_text(ggplot2::aes(x = col_x["p"],    label = P_disp), hjust = 1, size = font_size_data) +
        ggplot2::xlim(0,1.05) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.18, 0.25))) +
        ggplot2::coord_cartesian(clip = "off") +
        ggplot2::theme_void() +
        ggplot2::theme(
          plot.background = ggplot2::element_rect(fill = "white", color = NA),
          panel.background = ggplot2::element_rect(fill = "white", color = NA),
          plot.margin = ggplot2::margin(t = 0, r = base_margin$right * 4, b = 30, l = base_margin$left * 4)
        )
      p$plot <- p$plot + ggplot2::theme(legend.position = c(0.02, 0.02), legend.justification = c(0,0), plot.background = ggplot2::element_rect(fill = "white", color = NA))
      combined_success <- FALSE
      table_fraction <- min(0.30 + 0.035 * n_rows, 0.50)
      if (requireNamespace("cowplot", quietly = TRUE)) {
        p$plot <- cowplot::plot_grid(p$plot, table_plot, ncol = 1, rel_heights = c(1, table_fraction), align = "v")
        combined_success <- TRUE
      } else if (requireNamespace("gridExtra", quietly = TRUE)) {
        p$plot <- gridExtra::arrangeGrob(p$plot, table_plot, ncol = 1, heights = c(10, max(5, ceiling(10*table_fraction))))
        combined_success <- TRUE
      }
      if (!combined_success) {
        warning("Neither cowplot nor gridExtra available; falling back to caption mode.")
        stats_mode <- "caption"
      }
    }
    if (identical(stats_mode, "caption") && nrow(stats_tbl)) {
      caption_lines <- apply(stats_tbl, 1, function(r) paste(r["Test"], "HR=", ifelse(is.na(r["HR"]), "-", sprintf("%.2f", as.numeric(r["HR"]))), ifelse(nchar(r["CI"])>0, paste0("(", r["CI"], ")"), ""), "p=", ifelse(is.na(r["p"]), "-", format.pval(as.numeric(r["p"]), digits=3))))
      caption_text <- paste(caption_lines, collapse = "\n")
      extra_cm <- 0.25 + 0.18 * length(caption_lines)
      p$plot <- p$plot + ggplot2::theme(
        plot.margin = ggplot2::margin(t = base_margin$top, r = base_margin$right, b = base_margin$bottom + extra_cm + 0.8, l = base_margin$left, unit = "cm"),
        plot.caption = ggplot2::element_text(size = plot_params$hr_text_size * 1.05, hjust = 0, lineheight = 1.05, margin = ggplot2::margin(t = 10))
      ) + ggplot2::labs(caption = caption_text)
    }
    # Add Cox model summaries to return object for later extraction
    p$cox_unadj <- cox_unadj
    p$cox_adj <- cox_adj
    p$cont_unadj <- cont_unadj
    p$cont_adj <- cont_adj
    p$covariates_requested <- clinical_covariates
    p$covariates_retained_grouped <- retained_covariates
    p$covariates_retained_continuous <- if (exists("keep_covs_here2", inherits = FALSE)) keep_covs_here2 else character(0)
    p$split_method <- split_method
    p
  }
  
  # ============================================================================
  # DATA PREPARATION
  # ============================================================================
  
  if (verbose) cat("\n=== Preparing Data ===\n")
  
  # Filter samples with complete survival data
  complete_samples <- metadata %>%
    dplyr::filter(
      !is.na(.data[[event_col]]) & !is.na(.data[[time_col]]) &
      !.data[[event_col]] %in% c("_Not_Applicable_", "_Not_Available_") &
      !.data[[time_col]] %in% c("_Not_Applicable_", "_Not_Available_")
    ) %>%
    dplyr::pull(.data[[sample_col]])
  
  # Intersect with expression data samples
  expr_samples <- setdiff(colnames(expression_data), gene_col)
  available_samples <- intersect(expr_samples, complete_samples)
  
  if (verbose) {
    cat(sprintf("Total samples in metadata: %d\n", nrow(metadata)))
    cat(sprintf("Samples with complete survival data: %d\n", length(complete_samples)))
    cat(sprintf("Samples with expression data: %d\n", length(expr_samples)))
    cat(sprintf("Final samples for analysis: %d\n", length(available_samples)))
  }
  
  if (length(available_samples) == 0) {
    stop("No samples with both expression and survival metadata!")
  }
  
  # Select and z-score expression data
  selected_data <- expression_data %>%
    dplyr::select(dplyr::all_of(c(gene_col, available_samples))) %>%
    dplyr::filter(.data[[gene_col]] %in% genes)
  
  raw_gene_mat <- as.matrix(selected_data[, -1])
  rownames(raw_gene_mat) <- selected_data[[gene_col]]
  gene_mat <- t(scale(t(raw_gene_mat)))  # Z-score by gene for continuous modeling
  grouped_gene_mat <- if (split_method == "detectable") raw_gene_mat else gene_mat
  
  # Both forms use gene-wise z-scores. Sum preserves the established default;
  # mean is useful when reproducing analyses that report a per-gene average.
  signature_expression <- if (signature_aggregation == "mean") {
    colMeans(gene_mat, na.rm = TRUE)
  } else {
    colSums(gene_mat, na.rm = TRUE)
  }
  
  # Prepare survival data frame
  surv_df <- metadata %>%
    dplyr::filter(.data[[sample_col]] %in% available_samples)
  
  # Clean covariates if provided
  if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
    surv_df <- clean_covariates(surv_df, clinical_covariates)
  }
  
  # Parse time and event
  surv_df <- surv_df %>%
    dplyr::mutate(
      .time_raw = parse_time_numeric(.data[[time_col]]),
      event = dplyr::case_when(
        tolower(as.character(.data[[event_col]])) == "dead" ~ 1,
        tolower(as.character(.data[[event_col]])) == "alive" ~ 0,
        .data[[event_col]] == 1 ~ 1,
        .data[[event_col]] == 0 ~ 0,
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::filter(!is.na(.time_raw) & !is.na(event))
  
  # Convert days to months if needed
  raw_med <- median(surv_df$.time_raw, na.rm = TRUE)
  raw_max <- max(surv_df$.time_raw, na.rm = TRUE)
  is_days <- is.finite(raw_max) && (raw_max > 1200 || raw_med > 365)
  
  if (is_days) {
    if (verbose) message(sprintf("Detected DAYS (median=%.1f, max=%.1f) → converting to MONTHS.", raw_med, raw_max))
    surv_df <- surv_df %>% dplyr::mutate(time = .time_raw / 30.44)
  } else {
    if (verbose) message(sprintf("Detected MONTHS (median=%.1f, max=%.1f).", raw_med, raw_max))
    surv_df <- surv_df %>% dplyr::mutate(time = .time_raw)
  }
  
  surv_df <- surv_df %>%
    dplyr::select(-.time_raw) %>%
    dplyr::filter(!is.na(time), time > 0, !is.na(event))
  
  # Create base merged data frame
  merged_base <- surv_df %>%
    dplyr::rename(sample_id = !!sample_col) %>%
    dplyr::select(sample_id, time, event, 
                 dplyr::any_of(clinical_covariates))
  
  merged_with_sig <- merged_base %>%
    dplyr::mutate(signature_score = signature_expression[sample_id]) %>%
    dplyr::filter(!is.na(signature_score))  # Remove samples with NA signature scores
  
  if (verbose && nrow(merged_base) != nrow(merged_with_sig)) {
    cat(sprintf("  Note: %d samples removed due to NA signature scores\n", 
                nrow(merged_base) - nrow(merged_with_sig)))
  }
  
  # ============================================================================
  # INITIALIZE RESULTS
  # ============================================================================
  
  results <- list(
    results_univariate = NULL,
    results_adjusted = NULL,
    results_grouped_gene = NULL,
    results_signature = NULL,
    lasso_selected = NULL,
    ph_diagnostics = NULL,
    sample_groups = NULL,
    plots = list()
  )

  build_grouped_gene_export_rows <- function(cox_summary, feature_name, model_type,
                                             split_method, requested_covariates,
                                             retained_covariates) {
    if (is.null(cox_summary)) {
      return(NULL)
    }

    coefficient_rows <- rownames(cox_summary$coefficients)
    cluster_rows <- grep("^cluster", coefficient_rows)
    if (!length(cluster_rows)) {
      return(NULL)
    }

    describe_term <- function(term_name, split_method) {
      if (split_method == "detectable") {
        return(list(
          comparison = "Detectable vs Undetectable",
          reference_group = "Undetectable",
          comparison_group = "Detectable"
        ))
      }

      if (term_name == "clusterC") {
        return(list(
          comparison = "High vs Low",
          reference_group = "Low",
          comparison_group = "High"
        ))
      }

      if (term_name == "clusterB") {
        if (split_method == "median") {
          return(list(
            comparison = "High vs Low",
            reference_group = "Low",
            comparison_group = "High"
          ))
        }

        return(list(
          comparison = "Mid vs Low",
          reference_group = "Low",
          comparison_group = "Mid"
        ))
      }

      list(
        comparison = term_name,
        reference_group = NA_character_,
        comparison_group = NA_character_
      )
    }

    rows <- lapply(cluster_rows, function(idx) {
      term_name <- coefficient_rows[[idx]]
      labels <- describe_term(term_name, split_method)

      data.frame(
        feature = feature_name,
        comparison = labels$comparison,
        reference_group = labels$reference_group,
        comparison_group = labels$comparison_group,
        HR = unname(cox_summary$conf.int[idx, "exp(coef)"]),
        CI_lower = unname(cox_summary$conf.int[idx, "lower .95"]),
        CI_upper = unname(cox_summary$conf.int[idx, "upper .95"]),
        p_value = unname(cox_summary$coefficients[idx, "Pr(>|z|)"]),
        N = cox_summary$n,
        events = cox_summary$nevent,
        model_type = model_type,
        covariates_requested = if (length(requested_covariates)) paste(requested_covariates, collapse = ",") else "none",
        covariates_retained = if (length(retained_covariates)) paste(retained_covariates, collapse = ",") else "none",
        stringsAsFactors = FALSE
      )
    })

    dplyr::bind_rows(rows)
  }
  
  # ============================================================================
  # ANALYSIS MODULES
  # ============================================================================
  
  # --- MODULE 1: INDIVIDUAL GENE UNIVARIATE COX ---
  if ("individual" %in% steps_to_run) {
    if (verbose) cat("\n=== Individual Gene Analysis (Univariate) ===\n")
    
    res_univ <- lapply(genes, function(g) {
      expr_g <- get_expr_aligned(g, grouped_gene_mat, merged_base, available_samples)
      df <- merged_base %>% dplyr::mutate(expr = expr_g)
      if (all(is.na(df$expr))) return(NULL)
      
      fit <- try(survival::coxph(survival::Surv(time, event) ~ expr, data = df), silent = TRUE)
      if (inherits(fit, "try-error")) return(NULL)
      
      s <- summary(fit)
      data.frame(
        Gene = g,
        HR = unname(s$conf.int[, "exp(coef)"]),
        CI_Lower = unname(s$conf.int[, "lower .95"]),
        CI_Upper = unname(s$conf.int[, "upper .95"]),
        P_Value = unname(s$coefficients[, "Pr(>|z|)"])
      )
    })
    
    res_univ <- do.call(rbind, res_univ)
    if (!is.null(res_univ) && nrow(res_univ)) {
      res_univ$FDR <- stats::p.adjust(res_univ$P_Value, method = p_adjust_method)
      res_univ <- res_univ[order(res_univ$FDR, res_univ$P_Value), ]
      results$results_univariate <- res_univ
      
      if (verbose) print(res_univ, row.names = FALSE)
      if (export_csv && !is.null(output_dir)) {
        write.csv(res_univ, file.path(output_dir, "results_univariate.csv"), row.names = FALSE)
      }
    }
  }
  
  # --- MODULE 2: INDIVIDUAL GENE COVARIATE-ADJUSTED COX ---
  if ("individual+covariates" %in% steps_to_run) {
    if (verbose) cat("\n=== Individual Gene Analysis (Covariate-Adjusted) ===\n")
    
    if (is.null(clinical_covariates) || length(clinical_covariates) == 0) {
      warning("No clinical covariates specified; skipping adjusted analysis")
    } else {
      # Validate covariates
      keep_covs <- intersect(clinical_covariates, colnames(merged_with_sig))
      
      valid_covs <- c()
      for (cov in keep_covs) {
        values <- merged_with_sig[[cov]]
        if (is.numeric(values)) {
          var_val <- var(values, na.rm = TRUE)
          n_missing <- sum(is.na(values))
          if (!is.na(var_val) && var_val > 1e-10 && n_missing < 0.8 * length(values)) {
            valid_covs <- c(valid_covs, cov)
          }
        } else if (is.factor(values) || is.character(values)) {
          n_unique <- length(unique(na.omit(values)))
          n_missing <- sum(is.na(values))
          if (n_unique > 1 && n_missing < 0.8 * length(values)) {
            valid_covs <- c(valid_covs, cov)
          }
        }
      }
      
      if (verbose) cat("Valid covariates:", paste(valid_covs, collapse = ", "), "\n")
      
      res_adj <- lapply(genes, function(g) {
        expr_g <- get_expr_aligned(g, gene_mat, merged_base, available_samples)
        df <- merged_with_sig %>% dplyr::mutate(expr = expr_g)
        if (all(is.na(df$expr))) return(NULL)
        
        if (length(valid_covs)) {
          model_cols <- c("time", "event", "expr", valid_covs)
          df_complete <- df[, model_cols, drop = FALSE]
          df_complete <- df_complete[complete.cases(df_complete), ]
          rhs <- paste(c("expr", valid_covs), collapse = " + ")
        } else {
          df_complete <- df[!is.na(df$expr) & !is.na(df$time) & !is.na(df$event), ]
          rhs <- "expr"
        }
        
        if (nrow(df_complete) < 10) return(NULL)
        
        fit <- try(survival::coxph(as.formula(paste("survival::Surv(time, event) ~", rhs)), 
                                  data = df_complete), silent = TRUE)
        if (inherits(fit, "try-error")) return(NULL)
        
        s <- summary(fit)
        ridx <- match("expr", rownames(s$coefficients))
        if (is.na(ridx)) return(NULL)
        
        data.frame(
          Gene = g,
          HR = unname(s$conf.int[ridx, "exp(coef)"]),
          CI_Lower = unname(s$conf.int[ridx, "lower .95"]),
          CI_Upper = unname(s$conf.int[ridx, "upper .95"]),
          P_Value = unname(s$coefficients[ridx, "Pr(>|z|)"]),
          N_samples = nrow(df_complete)
        )
      })
      
      res_adj <- do.call(rbind, res_adj)
      if (!is.null(res_adj) && nrow(res_adj)) {
        res_adj$FDR <- stats::p.adjust(res_adj$P_Value, method = p_adjust_method)
        res_adj <- res_adj[order(res_adj$FDR, res_adj$P_Value), ]
        results$results_adjusted <- res_adj
        
        if (verbose) print(res_adj, row.names = FALSE)
        if (export_csv && !is.null(output_dir)) {
          write.csv(res_adj, file.path(output_dir, "results_adjusted.csv"), row.names = FALSE)
        }
      }
    }
  }
  
  # --- MODULE 3: SIGNATURE ANALYSIS (GROUPED KM + COX) ---
  if ("signature" %in% steps_to_run || "signature+covariates" %in% steps_to_run) {
    if (verbose) cat("\n=== Gene Signature Analysis ===\n")
    
    # Build signature groups
    sig_vec <- merged_with_sig$signature_score  # aligned to merged_with_sig rows
    if (split_method == "median") {
      med <- median(sig_vec, na.rm = TRUE)
      labs <- ifelse(sig_vec < med, "A", ifelse(sig_vec >= med, "B", NA))
    } else if (split_method == "quartiles") {
      # Use cut() to avoid dropping samples at boundaries
      q25 <- quantile(sig_vec, 0.25, na.rm = TRUE)
      q75 <- quantile(sig_vec, 0.75, na.rm = TRUE)
      breaks <- c(-Inf, q25, q75, Inf)
      labs <- cut(sig_vec, breaks = breaks, labels = c("A", "B", "C"), 
                  include.lowest = TRUE, right = TRUE)
      labs <- as.character(labs)
    } else if (split_method == "tertiles") {
      # Use cut() to avoid dropping samples at boundaries
      t1 <- quantile(sig_vec, 1/3, na.rm = TRUE)
      t2 <- quantile(sig_vec, 2/3, na.rm = TRUE)
      breaks <- c(-Inf, t1, t2, Inf)
      labs <- cut(sig_vec, breaks = breaks, labels = c("A", "B", "C"), 
                  include.lowest = TRUE, right = TRUE)
      labs <- as.character(labs)
    } else if (split_method == "percentiles") {
      lo <- quantile(sig_vec, percentile_cutoffs[1]/100, na.rm = TRUE)
      hi <- quantile(sig_vec, percentile_cutoffs[2]/100, na.rm = TRUE)
      labs <- ifelse(sig_vec <= lo, "A",
                     ifelse(sig_vec >= hi, "C",
                            ifelse(sig_vec > lo & sig_vec < hi, "B", NA)))
    } else if (split_method == "detectable") {
      labs <- ifelse(is.na(sig_vec), NA, ifelse(sig_vec <= detection_threshold, "A", "C"))
    }
    
    surv_data <- merged_with_sig %>%
      dplyr::mutate(cluster = labs) %>%
      dplyr::filter(!is.na(cluster))
    
    # Keep only requested groups
    keep <- character(0)
    if ("low" %in% include_groups) keep <- c(keep, "A")
    if ("mid" %in% include_groups && !split_method %in% c("median", "detectable")) keep <- c(keep, "B")
    if ("high" %in% include_groups) keep <- c(keep, if (split_method == "median") "B" else "C")
    
    surv_data <- surv_data %>%
      dplyr::filter(cluster %in% keep) %>%
      dplyr::mutate(cluster = factor(cluster, levels = keep))
    
    if (nrow(surv_data) == 0 || nlevels(surv_data$cluster) < 2) {
      warning("Not enough strata for signature KM plot")
    } else {
      # Fit KM and Cox models
      fit <- survival::survfit(survival::Surv(time, event) ~ cluster, data = surv_data)
      cox_unadj <- survival::coxph(survival::Surv(time, event) ~ cluster, data = surv_data)
      
      # Adjusted Cox if covariates available
      cox_adj <- NULL
      valid_covs_local <- character(0)
      if ("signature+covariates" %in% steps_to_run && 
          !is.null(clinical_covariates) && length(clinical_covariates) > 0) {
        
        valid_covs_local <- intersect(clinical_covariates, 
                                      setdiff(colnames(surv_data), c("time", "event", "cluster")))
        if (length(valid_covs_local)) {
          f_adj <- as.formula(paste("survival::Surv(time, event) ~ cluster +", 
                                   paste(valid_covs_local, collapse = " + ")))
          surv_adj <- surv_data[complete.cases(surv_data[, all.vars(f_adj)]), ]
          if (nrow(surv_adj) >= 10) {
            cox_adj <- try(survival::coxph(f_adj, data = surv_adj), silent = TRUE)
            if (inherits(cox_adj, "try-error")) cox_adj <- NULL
          }
        }
      }

      # PH diagnostics and time-dependent modelling for grouped signature (Low/Mid/High)
      ph_test_signature <- td_model_signature <- NULL
      if (run_ph_checks) {
        ph_candidate <- try(survival::cox.zph(cox_unadj), silent = TRUE)
        if (!inherits(ph_candidate, "try-error")) {
          ph_test_signature <- ph_candidate
        }
      }

      if (run_time_dependent && !is.null(ph_test_signature) && !inherits(ph_test_signature, "try-error")) {
        ph_tab <- ph_test_signature$table
        p_global_sig <- if ("GLOBAL" %in% rownames(ph_tab)) ph_tab["GLOBAL", "p"] else NA_real_
        cluster_rows <- grep("^cluster", rownames(ph_tab))
        p_cluster_sig <- NA_real_
        if (length(cluster_rows) > 0) {
          high_rows <- cluster_rows[grep("clusterC$|clusterHigh$", rownames(ph_tab)[cluster_rows])]
          target_row <- if (length(high_rows) > 0) high_rows[1] else cluster_rows[1]
          p_cluster_sig <- ph_tab[target_row, "p"]
        }
        run_td_sig <- (is.finite(p_global_sig) && p_global_sig < ph_alpha) ||
                       (is.finite(p_cluster_sig) && p_cluster_sig < ph_alpha)

        if (run_td_sig) {
          keep_covs_td <- if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
            intersect(clinical_covariates, setdiff(colnames(surv_data), c("time", "event", "cluster")))
          } else character(0)

          if (length(keep_covs_td) > 0) {
            ok_td <- vapply(keep_covs_td, function(cv) {
              v <- surv_data[[cv]]
              if (is.numeric(v)) var(v, na.rm = TRUE) > 0
              else if (is.factor(v) || is.character(v)) length(unique(na.omit(v))) > 1
              else FALSE
            }, logical(1))
            keep_covs_td <- keep_covs_td[ok_td]
          }

          surv_sig_td <- surv_data
          if (!is.factor(surv_sig_td$cluster)) {
            surv_sig_td$cluster <- factor(surv_sig_td$cluster)
          }

          cluster_levels <- levels(surv_sig_td$cluster)
          if (length(cluster_levels) > 1) {
            surv_sig_td$time_early <- as.numeric(surv_sig_td$time <= cut_months)
            surv_sig_td$time_late  <- as.numeric(surv_sig_td$time >  cut_months)

            for (lev in cluster_levels[-1]) {
              surv_sig_td[[paste0("cluster_", lev, "_early")]] <-
                as.numeric(surv_sig_td$cluster == lev) * surv_sig_td$time_early
              surv_sig_td[[paste0("cluster_", lev, "_late")]] <-
                as.numeric(surv_sig_td$cluster == lev) * surv_sig_td$time_late
            }

            base_terms <- if (length(keep_covs_td) > 0) paste(keep_covs_td, collapse = " + ") else "1"
            td_terms <- paste0("cluster_", cluster_levels[-1], "_early", collapse = " + ")
            td_terms <- paste(td_terms, "+", paste0("cluster_", cluster_levels[-1], "_late", collapse = " + "))
            fml_td <- as.formula(paste("survival::Surv(time, event) ~", base_terms, "+", td_terms))

            td_fit_sig <- try(survival::coxph(fml_td, data = surv_sig_td), silent = TRUE)
            if (!inherits(td_fit_sig, "try-error")) {
              td_model_signature <- summary(td_fit_sig)
            }
          }
        }
      }
      
      # Continuous signature Cox models can use either the plotted groups or
      # every eligible sample. The latter is needed for analyses that retain
      # the middle quantile for continuous effect estimation.
      continuous_surv_data <- if (signature_continuous_cohort == "all") merged_with_sig else surv_data
      sig_cont_unadj_fit <- try(survival::coxph(survival::Surv(time, event) ~ signature_score, data = continuous_surv_data), silent = TRUE)
      sig_cont_unadj_sum <- if (!inherits(sig_cont_unadj_fit, "try-error")) summary(sig_cont_unadj_fit) else NULL
      sig_cont_adj_sum <- NULL
      if ("signature+covariates" %in% steps_to_run && !is.null(clinical_covariates) && length(clinical_covariates) > 0) {
          valid_covs_local2 <- intersect(clinical_covariates, setdiff(colnames(continuous_surv_data), c("time","event","cluster","signature_score")))
          if (length(valid_covs_local2)) {
            f_sig_cont_adj <- as.formula(paste("survival::Surv(time, event) ~ signature_score +", paste(valid_covs_local2, collapse = " + ")))
            sig_cont_adj_fit <- try(survival::coxph(f_sig_cont_adj, data = continuous_surv_data), silent = TRUE)
            if (!inherits(sig_cont_adj_fit, "try-error")) sig_cont_adj_sum <- summary(sig_cont_adj_fit)
          }
      }
      # Create plot
      # For signature plots, override to side panel if available
      sig_plot_params <- plot_params
  sig_plot_params$stats_panel_mode <- plot_params$stats_panel_mode_signature %||% "bottom_table"  # enforce bottom table for signature
      old_plot_params <- plot_params
      plot_params <<- sig_plot_params
      p <- draw_km_plot(
        fit, surv_data,
        cox_unadj = summary(cox_unadj),
        cox_adj = if (!is.null(cox_adj)) summary(cox_adj) else NULL,
        cont_unadj = sig_cont_unadj_sum,
        cont_adj = sig_cont_adj_sum,
        main_title = paste(plot_params$plot_title, "-", plot_params$signature_label %||% "Signature"),
        local_plot_params = sig_plot_params,
        attempted_adj = ("signature+covariates" %in% steps_to_run && !is.null(clinical_covariates) && length(clinical_covariates) > 0),
        adj_reason = if (is.null(cox_adj) && ("signature+covariates" %in% steps_to_run)) {
          if (length(valid_covs_local) == 0) {
            "no valid covariates"
          } else if (exists("surv_adj") && nrow(surv_adj) < 10) {
            paste0("insufficient complete cases (n=", nrow(surv_adj), ")")
          } else {
            "coxph error or convergence failure"
          }
        } else NULL,
        retained_covariates = valid_covs_local
      )
      plot_params <<- old_plot_params
      
      if (!is.null(ph_test_signature) && !inherits(ph_test_signature, "try-error")) {
        p$ph_test <- ph_test_signature

        sig_ph_tab <- ph_test_signature$table
        sig_p_global <- if ("GLOBAL" %in% rownames(sig_ph_tab)) sig_ph_tab["GLOBAL", "p"] else NA_real_
        cluster_rows <- grep("^cluster", rownames(sig_ph_tab))
        sig_p_signature <- NA_real_
        if (length(cluster_rows) > 0) {
          high_rows <- cluster_rows[grep("clusterC$|clusterHigh$", rownames(sig_ph_tab)[cluster_rows])]
          target_row <- if (length(high_rows) > 0) high_rows[1] else cluster_rows[1]
          sig_p_signature <- sig_ph_tab[target_row, "p"]
        }

        if (is.null(results$ph_diagnostics)) results$ph_diagnostics <- list()
        if (is.null(results$ph_diagnostics$grouped)) results$ph_diagnostics$grouped <- list()
        results$ph_diagnostics$p_global <- sig_p_global
        results$ph_diagnostics$p_signature <- sig_p_signature
        results$ph_diagnostics$ph_test <- ph_test_signature
        results$ph_diagnostics$grouped$ph_test <- ph_test_signature
        results$ph_diagnostics$grouped$p_global <- sig_p_global
        results$ph_diagnostics$grouped$p_signature <- sig_p_signature
      }

      if (!is.null(td_model_signature)) {
        p$td_model <- td_model_signature
        if (is.null(results$ph_diagnostics)) results$ph_diagnostics <- list()
        if (is.null(results$ph_diagnostics$grouped)) results$ph_diagnostics$grouped <- list()
        results$ph_diagnostics$time_dependent <- td_model_signature
        results$ph_diagnostics$grouped$time_dependent <- td_model_signature
      }

      results$plots$signature <- p
      
      if (verbose) {
        print(p)
        cat("\nCox Proportional Hazards (signature):\n")
        print(summary(cox_unadj))
        if (!is.null(cox_adj)) {
          cat("\nCox Proportional Hazards (signature, adjusted):\n")
          print(summary(cox_adj))
        }
      }
      
      if (save_plots) {
        ggplot2::ggsave(
          file.path(output_dir, "signature_km_plot.png"),
          plot = p$plot, width = 10, height = 8.5, dpi = 300, bg = "white"
        )
      }
      
      # Store signature Cox results
      results$results_signature <- list(
        unadjusted = summary(cox_unadj),
        adjusted = if (!is.null(cox_adj)) summary(cox_adj) else NULL
      )
      
      # --- GENERATE SIGNATURE EARLY/LATE WINDOW PLOTS (UNCONDITIONAL) ---
      if (run_time_dependent) {
        if (verbose) cat("\n=== Generating Signature Early/Late Window Plots ===\n")
        cat("Creating signature KM plots with cutpoint at", cut_months, "months\n")
        
        # Use merged_with_sig as base data (all samples with signature scores)
        window_data <- merged_with_sig[complete.cases(merged_with_sig[, c("time", "event", "signature_score")]), ]

        if (nrow(window_data) >= 20) {
          # Derive Low/High groups using the same split method as the primary analysis
          cluster_two <- make_clusters(
            window_data$signature_score,
            split_method = split_method,
            percentile_cutoffs = percentile_cutoffs,
            include_groups = c("low", "high"),
            detection_threshold = detection_threshold
          )
          cluster_chr <- as.character(cluster_two)
          cluster_chr[cluster_chr == "A"] <- "Low"
          cluster_chr[cluster_chr %in% c("B", "C")] <- "High"
          cluster_chr[!(cluster_chr %in% c("Low", "High"))] <- NA_character_
          window_data <- window_data %>%
            dplyr::mutate(
              cluster_low_high = factor(cluster_chr, levels = c("Low", "High"))
            )

          # Filter valid covariates for adjusted models (if available)
          valid_covs_window <- c()
          if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
            valid_covs <- intersect(clinical_covariates, colnames(window_data))
            for (cv in valid_covs) {
              values <- window_data[[cv]]
              if (is.numeric(values) && var(values, na.rm = TRUE) > 1e-10) {
                valid_covs_window <- c(valid_covs_window, cv)
              } else if ((is.factor(values) || is.character(values)) && 
                        length(unique(na.omit(values))) > 1) {
                valid_covs_window <- c(valid_covs_window, cv)
              }
            }
          }
          
          # Early window
          early_data <- window_data %>%
            dplyr::mutate(
              time = pmin(time, cut_months),
              event = ifelse(time <= cut_months & event == 1, 1, 0),
              cluster = cluster_low_high
            ) %>%
            dplyr::filter(!is.na(cluster))
          
          if (nrow(early_data) >= 10 && nlevels(early_data$cluster) == 2) {
            fit_early <- survival::survfit(survival::Surv(time, event) ~ cluster, 
                                           data = early_data)
            cox_early_unadj <- survival::coxph(survival::Surv(time, event) ~ cluster, 
                                               data = early_data)
            cox_early_cont_unadj <- try(survival::coxph(survival::Surv(time, event) ~ signature_score, 
                                                        data = early_data), silent = TRUE)
            if (inherits(cox_early_cont_unadj, "try-error")) cox_early_cont_unadj <- NULL
            
            # Adjusted models (if covariates available)
            valid_covs_early <- intersect(valid_covs_window, colnames(early_data))
            early_data_complete <- if (length(valid_covs_early) > 0) {
              early_data[complete.cases(early_data[, c("time", "event", "cluster", "signature_score", valid_covs_early)]), ]
            } else {
              early_data
            }
            
            cox_early_adj <- NULL
            if (length(valid_covs_early) > 0 && nrow(early_data_complete) >= 10) {
              fml_early <- as.formula(paste("survival::Surv(time, event) ~ cluster +", 
                                           paste(valid_covs_early, collapse = " + ")))
              cox_early_adj <- try(survival::coxph(fml_early, data = early_data_complete), silent = TRUE)
              if (inherits(cox_early_adj, "try-error")) cox_early_adj <- NULL
            }
            
            cox_early_cont_adj <- NULL
            if (length(valid_covs_early) > 0 && nrow(early_data_complete) >= 10) {
              fml_early_cont <- as.formula(paste("survival::Surv(time, event) ~ signature_score +", 
                                                 paste(valid_covs_early, collapse = " + ")))
              cox_early_cont_adj <- try(survival::coxph(fml_early_cont, data = early_data_complete), silent = TRUE)
              if (inherits(cox_early_cont_adj, "try-error")) cox_early_cont_adj <- NULL
            }
            
            early_params <- plot_params
            early_params$stats_panel_mode <- "bottom_table"
            early_params$title_size <- max(10, (plot_params$title_size %||% 15) * 0.85)
            early_title <- wrap_title(sprintf("%s - Early Window (≤%d months)", plot_params$plot_title, cut_months))
            p_early <- draw_km_plot(
              fit_early, early_data,
              cox_unadj = summary(cox_early_unadj),
              cox_adj = if (!is.null(cox_early_adj)) summary(cox_early_adj) else NULL,
              cont_unadj = if (!is.null(cox_early_cont_unadj)) summary(cox_early_cont_unadj) else NULL,
              cont_adj   = if (!is.null(cox_early_cont_adj)) summary(cox_early_cont_adj) else NULL,
              main_title = early_title,
              local_plot_params = early_params,
              attempted_adj = length(valid_covs_early) > 0,
              adj_reason = if (is.null(cox_early_adj) && length(valid_covs_early) > 0) "coxph error or insufficient data" else NULL
            )
            
            results$plots$signature_early <- p_early
            if (verbose) print(p_early)
            
            # Store early window Cox results
            results$results_signature_early <- list(
              unadjusted = summary(cox_early_unadj),
              adjusted = if (!is.null(cox_early_adj)) summary(cox_early_adj) else NULL,
              n_samples = nrow(early_data),
              n_complete = if (!is.null(cox_early_adj)) nrow(early_data_complete) else nrow(early_data)
            )
            
            if (save_plots) {
              ggplot2::ggsave(
                file.path(output_dir, sprintf("signature_early_%dm.png", cut_months)),
                plot = p_early$plot, width = 10, height = 8.5, dpi = 300, bg = "white"
              )
            }
          }
          
          # Late window
          late_data <- window_data %>%
            dplyr::filter(time > cut_months) %>%
            dplyr::mutate(
              cluster = cluster_low_high
            ) %>%
            dplyr::filter(!is.na(cluster))
          
          if (nrow(late_data) >= 10 && nlevels(late_data$cluster) == 2) {
            fit_late <- survival::survfit(survival::Surv(time, event) ~ cluster, 
                                          data = late_data)
            cox_late_unadj <- survival::coxph(survival::Surv(time, event) ~ cluster, 
                                              data = late_data)
            cox_late_cont_unadj <- try(survival::coxph(survival::Surv(time, event) ~ signature_score, 
                                                       data = late_data), silent = TRUE)
            if (inherits(cox_late_cont_unadj, "try-error")) cox_late_cont_unadj <- NULL
            
            valid_covs_late <- intersect(valid_covs_window, colnames(late_data))
            late_data_complete <- if (length(valid_covs_late) > 0) {
              late_data[complete.cases(late_data[, c("time", "event", "cluster", "signature_score", valid_covs_late)]), ]
            } else {
              late_data
            }
            
            cox_late_adj <- NULL
            if (length(valid_covs_late) > 0 && nrow(late_data_complete) >= 10) {
              fml_late <- as.formula(paste("survival::Surv(time, event) ~ cluster +", 
                                          paste(valid_covs_late, collapse = " + ")))
              cox_late_adj <- try(survival::coxph(fml_late, data = late_data_complete), silent = TRUE)
              if (inherits(cox_late_adj, "try-error")) cox_late_adj <- NULL
            }
            
            cox_late_cont_adj <- NULL
            if (length(valid_covs_late) > 0 && nrow(late_data_complete) >= 10) {
              fml_late_cont <- as.formula(paste("survival::Surv(time, event) ~ signature_score +", 
                                                paste(valid_covs_late, collapse = " + ")))
              cox_late_cont_adj <- try(survival::coxph(fml_late_cont, data = late_data_complete), silent = TRUE)
              if (inherits(cox_late_cont_adj, "try-error")) cox_late_cont_adj <- NULL
            }
            
            late_params <- plot_params
            late_params$stats_panel_mode <- "bottom_table"
            late_params$title_size <- max(10, (plot_params$title_size %||% 15) * 0.85)
            late_title <- wrap_title(sprintf("%s - Late Window (>%d months)", plot_params$plot_title, cut_months))
            p_late <- draw_km_plot(
              fit_late, late_data,
              cox_unadj = summary(cox_late_unadj),
              cox_adj = if (!is.null(cox_late_adj)) summary(cox_late_adj) else NULL,
              cont_unadj = if (!is.null(cox_late_cont_unadj)) summary(cox_late_cont_unadj) else NULL,
              cont_adj   = if (!is.null(cox_late_cont_adj)) summary(cox_late_cont_adj) else NULL,
              main_title = late_title,
              local_plot_params = late_params,
              attempted_adj = length(valid_covs_late) > 0,
              adj_reason = if (is.null(cox_late_adj) && length(valid_covs_late) > 0) "coxph error or insufficient data" else NULL
            )
            
            results$plots$signature_late <- p_late
            if (verbose) print(p_late)
            
            # Store late window Cox results
            results$results_signature_late <- list(
              unadjusted = summary(cox_late_unadj),
              adjusted = if (!is.null(cox_late_adj)) summary(cox_late_adj) else NULL,
              n_samples = nrow(late_data),
              n_complete = if (!is.null(cox_late_adj)) nrow(late_data_complete) else nrow(late_data)
            )
            
            if (save_plots) {
              ggplot2::ggsave(
                file.path(output_dir, sprintf("signature_late_%dm.png", cut_months)),
                plot = p_late$plot, width = 10, height = 8.5, dpi = 300, bg = "white"
              )
            }
          }
        }
      }
      
      # --- PROPORTIONAL HAZARDS CHECK AND TIME-DEPENDENT ANALYSIS ---
      if (run_ph_checks && "signature+covariates" %in% steps_to_run) {
        if (verbose) cat("\n=== Proportional Hazards Check ===\n")
        
        # Fit continuous signature + covariates model
        if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
          valid_covs <- intersect(clinical_covariates, colnames(merged_with_sig))
          
          # Filter valid covariates
          valid_covs_filtered <- c()
          for (cv in valid_covs) {
            values <- merged_with_sig[[cv]]
            if (is.numeric(values) && var(values, na.rm = TRUE) > 1e-10) {
              valid_covs_filtered <- c(valid_covs_filtered, cv)
            } else if ((is.factor(values) || is.character(values)) && 
                      length(unique(na.omit(values))) > 1) {
              valid_covs_filtered <- c(valid_covs_filtered, cv)
            }
          }
          
          if (length(valid_covs_filtered) > 0) {
            fml <- as.formula(paste("survival::Surv(time, event) ~ signature_score +", 
                                   paste(valid_covs_filtered, collapse = " + ")))
            
            model_data <- merged_with_sig[, c("time", "event", "signature_score", 
                                              valid_covs_filtered), drop = FALSE]
            model_data <- model_data[complete.cases(model_data), ]
            
            if (nrow(model_data) >= 20) {
              fit_sig_cov <- try(survival::coxph(fml, data = model_data), silent = TRUE)
              
              if (!inherits(fit_sig_cov, "try-error")) {
                ph_result <- try(survival::cox.zph(fit_sig_cov), silent = TRUE)
                
                if (!inherits(ph_result, "try-error")) {
                  if (verbose) {
                    cat("PH test global p =", signif(ph_result$table["GLOBAL", "p"], 3), "\n")
                    print(ph_result)
                  }
                  
                  # Check if PH assumption violated
                  p_sig <- if ("signature_score" %in% rownames(ph_result$table)) {
                    ph_result$table["signature_score", "p"]
                  } else NA_real_
                  p_glob <- ph_result$table["GLOBAL", "p"]
                  
                  run_tt <- (is.finite(p_sig) && p_sig < ph_alpha) || 
                           (is.finite(p_glob) && p_glob < ph_alpha)
                  
                  # Store PH diagnostics (even if no violation)
                  cont_p_global <- if ("GLOBAL" %in% rownames(ph_result$table)) ph_result$table["GLOBAL", "p"] else NA_real_
                  cont_p_signature <- if ("signature_score" %in% rownames(ph_result$table)) ph_result$table["signature_score", "p"] else NA_real_

                  if (is.null(results$ph_diagnostics)) results$ph_diagnostics <- list()
                  if (is.null(results$ph_diagnostics$continuous)) results$ph_diagnostics$continuous <- list()

                  results$ph_diagnostics$continuous$ph_test <- ph_result
                  results$ph_diagnostics$continuous$p_global <- cont_p_global
                  results$ph_diagnostics$continuous$p_signature <- cont_p_signature

                  if (is.null(results$ph_diagnostics$ph_test)) {
                    results$ph_diagnostics$ph_test <- ph_result
                  }
                  if (is.null(results$ph_diagnostics$p_global) || !is.finite(results$ph_diagnostics$p_global)) {
                    results$ph_diagnostics$p_global <- cont_p_global
                  }
                  if (is.null(results$ph_diagnostics$p_signature) || !is.finite(results$ph_diagnostics$p_signature)) {
                    results$ph_diagnostics$p_signature <- cont_p_signature
                  }

                  # Run time-dependent analysis if PH violation detected
                  if (run_tt && run_time_dependent) {
                    if (verbose) cat("\n=== Time-Dependent Analysis ===\n")
                    cat("PH violation detected (sig p =", signif(p_sig, 3), 
                        ", global p =", signif(p_glob, 3), ")\n")
                    cat("Fitting early/late model with cutpoint at", cut_months, "months\n")
                    
                    # Fit time-dependent Cox model
                    fml_tt <- as.formula(paste("survival::Surv(time, event) ~", 
                                              paste(c(valid_covs_filtered, "tt(signature_score)"), 
                                                   collapse = " + ")))
                    
                    fit_tt <- survival::coxph(
                      fml_tt,
                      data = model_data,
                      tt = function(x, t, ...) cbind(
                        sig_early = x * (t <= cut_months),
                        sig_late = x * (t > cut_months)
                      )
                    )
                    
                    if (verbose) {
                      cat("\nTime-dependent Cox model:\n")
                      print(summary(fit_tt))
                    }

                    # Augment stored diagnostics with time-dependent summary
                    if (is.null(results$ph_diagnostics)) results$ph_diagnostics <- list()
                    if (is.null(results$ph_diagnostics$continuous)) results$ph_diagnostics$continuous <- list()

                    results$ph_diagnostics$continuous$time_dependent <- summary(fit_tt)

                    if (is.null(results$ph_diagnostics$time_dependent)) {
                      results$ph_diagnostics$time_dependent <- summary(fit_tt)
                    }
                  } else {
                    if (verbose) cat("PH assumption holds (p =", signif(p_glob, 3), ")\n")
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  # --- MODULE 4: LASSO-COX SELECTION ---
  if (run_lasso && "lasso" %in% steps_to_run) {
    if (verbose) cat("\n=== LASSO-COX Selection ===\n")
    
    keep_idx <- with(merged_with_sig, !is.na(time) & !is.na(event) & 
                    is.finite(time) & time > 0)
    
    if (!any(keep_idx)) {
      warning("No samples with positive follow-up time for LASSO")
    } else {
      samps_keep <- merged_with_sig$sample_id[keep_idx]
      genes_present <- intersect(genes, rownames(gene_mat))
      
      if (length(genes_present) > 0) {
        X <- t(gene_mat[genes_present, samps_keep, drop = FALSE])
        
        # Add covariates if requested
        if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
          keep_covs <- intersect(clinical_covariates, colnames(merged_with_sig))
          
          if (length(keep_covs)) {
            Z_list <- lapply(keep_covs, function(cv) {
              col <- merged_with_sig[keep_idx, cv, drop = TRUE]
              if (is.factor(col)) col <- as.numeric(col)
              if (is.character(col)) col <- as.numeric(as.factor(col))
              col <- as.numeric(col)
              if (any(is.na(col))) {
                med_val <- median(col, na.rm = TRUE)
                if (is.na(med_val)) med_val <- 0
                col[is.na(col)] <- med_val
              }
              col
            })
            
            Z <- do.call(cbind, Z_list)
            colnames(Z) <- make.names(keep_covs, unique = TRUE)
            x_full <- cbind(X, Z)
            pen <- c(rep(1, ncol(X)), rep(0, ncol(Z)))
          } else {
            x_full <- X
            pen <- rep(1, ncol(X))
          }
        } else {
          x_full <- X
          pen <- rep(1, ncol(X))
        }
        
        y <- survival::Surv(merged_with_sig$time[keep_idx], 
                           merged_with_sig$event[keep_idx])
        
        set.seed(1)
        cvfit <- glmnet::cv.glmnet(
          x_full, y, family = "cox", alpha = lasso_alpha,
          nfolds = lasso_nfolds, penalty.factor = pen, standardize = FALSE
        )
        
        lambda_pick <- if (lasso_lambda == "min") cvfit$lambda.min else cvfit$lambda.1se
        cf <- as.matrix(glmnet::coef.glmnet(cvfit$glmnet.fit, s = lambda_pick))
        nz <- which(abs(cf[, 1]) > 0)
        
        if (length(nz)) {
          nz_names <- rownames(cf)[nz]
          selected_genes <- intersect(nz_names, colnames(X))
          results$lasso_selected <- selected_genes
          
          if (verbose) {
            cat("LASSO selected", length(selected_genes), "genes:\n")
            cat(paste(selected_genes, collapse = ", "), "\n")
          }
        } else {
          if (verbose) cat("LASSO selected no genes\n")
        }
      }
    }
  }
  
  # --- MODULE 5: INDIVIDUAL GENE KM PLOTS ---
  if ("km_individual" %in% steps_to_run) {
    if (verbose) cat("\n=== Individual Gene KM Plots ===\n")
  genes_to_plot <- genes
  genes_for_combination <- genes  # track full set actually plotted if we fallback
    grouped_gene_rows <- list()
    if (!is.null(results$results_adjusted) && nrow(results$results_adjusted) > 0) {
      sig_genes <- results$results_adjusted$Gene[results$results_adjusted$FDR < fdr_threshold]
      if (length(sig_genes) > 0) {
        genes_to_plot <- intersect(genes_to_plot, sig_genes)
        if (verbose) cat("Plotting", length(genes_to_plot), "significant genes\n")
      } else {
        if (verbose) cat("No genes passed FDR <", fdr_threshold, "; plotting all", length(genes_to_plot), "genes instead.\n")
        genes_for_combination <- genes_to_plot
      }
    }
    for (g in genes_to_plot) {
      if (!(g %in% rownames(gene_mat))) { if (verbose) message("Skipping ", g, " (not in expression matrix)"); next }
      expr_g <- get_expr_aligned(g, grouped_gene_mat, merged_base, available_samples)
      surv_gene <- merged_base %>% dplyr::mutate(cluster = make_clusters(expr_g, split_method, percentile_cutoffs, include_groups, detection_threshold)) %>% dplyr::filter(!is.na(cluster))
      if (nrow(surv_gene) == 0 || nlevels(surv_gene$cluster) < 2) { if (verbose) message("Skipping ", g, " (insufficient strata)"); next }
      fit_g <- survival::survfit(survival::Surv(time, event) ~ cluster, data = surv_gene)
      cox_unadj <- summary(survival::coxph(survival::Surv(time, event) ~ cluster, data = surv_gene))
      cox_adj <- NULL
      retained_covariates <- character(0)
      if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
        keep_covs_here <- intersect(clinical_covariates, setdiff(colnames(surv_gene), c("time", "event", "cluster")))
        if (length(keep_covs_here)) {
          ok <- vapply(keep_covs_here, function(cv) { v <- surv_gene[[cv]]; if (is.numeric(v)) { var(v, na.rm = TRUE) > 0 } else if (is.factor(v) || is.character(v)) { length(unique(na.omit(v))) > 1 } else FALSE }, logical(1))
          keep_covs_here <- keep_covs_here[ok]
          retained_covariates <- keep_covs_here
          if (length(keep_covs_here)) {
            fml_adj <- as.formula(paste("survival::Surv(time, event) ~ cluster +", paste(keep_covs_here, collapse = " + ")))
            surv_gene_adj <- surv_gene[, c("time", "event", "cluster", keep_covs_here), drop = FALSE]
            surv_gene_adj <- surv_gene_adj[complete.cases(surv_gene_adj), , drop = FALSE]
            cox_adj_try <- if (nrow(surv_gene_adj) >= 10) try(survival::coxph(fml_adj, data = surv_gene_adj), silent = TRUE) else NULL
            if (!inherits(cox_adj_try, "try-error") && !is.null(cox_adj_try)) {
              cox_adj <- summary(cox_adj_try)
            }
          }
        }
      }
      grouped_gene_rows[[length(grouped_gene_rows) + 1]] <- build_grouped_gene_export_rows(
        cox_summary = cox_unadj,
        feature_name = g,
        model_type = "grouped_unadjusted",
        split_method = split_method,
        requested_covariates = if (is.null(clinical_covariates)) character(0) else clinical_covariates,
        retained_covariates = character(0)
      )
      grouped_gene_rows[[length(grouped_gene_rows) + 1]] <- build_grouped_gene_export_rows(
        cox_summary = cox_adj,
        feature_name = g,
        model_type = "grouped_adjusted",
        split_method = split_method,
        requested_covariates = if (is.null(clinical_covariates)) character(0) else clinical_covariates,
        retained_covariates = retained_covariates
      )
      gene_params <- plot_params
      gene_params$stats_panel_mode <- "bottom_table"
      # Continuous univariate and adjusted Cox (per-unit expression) for display
      cont_unadj_fit <- try(survival::coxph(survival::Surv(time, event) ~ expr_g, data = surv_gene), silent = TRUE)
      cont_unadj_sum <- if (!inherits(cont_unadj_fit, "try-error")) summary(cont_unadj_fit) else NULL
      cont_adj_sum <- NULL
      keep_covs_here2 <- character(0)
      if (!is.null(cox_adj)) { # reuse adjusted covariates if available
        if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
          keep_covs_here2 <- intersect(clinical_covariates, setdiff(colnames(surv_gene), c("time","event")))
          if (length(keep_covs_here2)) {
            ok2 <- vapply(keep_covs_here2, function(cv) {
              v <- surv_gene[[cv]]; if (is.numeric(v)) var(v,na.rm=TRUE)>0 else if (is.factor(v)||is.character(v)) length(unique(na.omit(v)))>1 else FALSE
            }, logical(1))
            keep_covs_here2 <- keep_covs_here2[ok2]
            if (length(keep_covs_here2)) {
              f_cont_adj <- as.formula(paste("survival::Surv(time, event) ~ expr_val +", paste(keep_covs_here2, collapse=" + ")))
              surv_gene_cont_adj <- surv_gene[, c("time", "event", keep_covs_here2), drop = FALSE]
              surv_gene_cont_adj$expr_val <- expr_g
              surv_gene_cont_adj <- surv_gene_cont_adj[complete.cases(surv_gene_cont_adj), , drop = FALSE]
              cont_adj_fit <- if (nrow(surv_gene_cont_adj) >= 10) try(survival::coxph(f_cont_adj, data = surv_gene_cont_adj), silent = TRUE) else NULL
              if (!inherits(cont_adj_fit, "try-error") && !is.null(cont_adj_fit)) {
                cont_adj_sum <- summary(cont_adj_fit)
              }
            }
          }
        }
      }
      # Run PH tests and time-dependent models for individual genes
      ph_test_gene <- td_model_gene <- NULL
      if (run_ph_checks) {
        cox_fit_unadj <- survival::coxph(survival::Surv(time, event) ~ cluster, data = surv_gene)
        ph_test_gene <- try(survival::cox.zph(cox_fit_unadj), silent = TRUE)
      }
      
      # If PH violation detected and time-dependent analysis requested, run TD model
      if (run_time_dependent && !is.null(ph_test_gene) && !inherits(ph_test_gene, "try-error")) {
        ph_p_cluster <- if ("cluster" %in% rownames(ph_test_gene$table)) ph_test_gene$table["cluster", "p"] else NA_real_
        ph_p_global <- ph_test_gene$table["GLOBAL", "p"]
        
        # Run TD model if PH violation detected (p < 0.05)
        if (!is.na(ph_p_global) && ph_p_global < 0.05) {
          # Build time-dependent formula
          keep_covs_td <- if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
            intersect(clinical_covariates, setdiff(colnames(surv_gene), c("time", "event", "cluster")))
          } else character(0)
          
          if (length(keep_covs_td) > 0) {
            ok_td <- vapply(keep_covs_td, function(cv) {
              v <- surv_gene[[cv]]
              if (is.numeric(v)) var(v, na.rm = TRUE) > 0
              else if (is.factor(v) || is.character(v)) length(unique(na.omit(v))) > 1
              else FALSE
            }, logical(1))
            keep_covs_td <- keep_covs_td[ok_td]
          }
          
          base_terms <- if (length(keep_covs_td) > 0) paste(keep_covs_td, collapse = " + ") else "1"
          
          # Create time-dependent covariates directly in the data (60 months cutoff)
          cut_months <- 60
          surv_gene_td <- surv_gene
          # Create binary indicator for early vs late period
          surv_gene_td$time_early <- as.numeric(surv_gene_td$time <= cut_months)
          surv_gene_td$time_late <- as.numeric(surv_gene_td$time > cut_months)
          
          # Create interaction terms: cluster × time window
          # For factor cluster with levels A, B, C:
          # cluster_B_early = (cluster==B) * time_early
          # cluster_C_early = (cluster==C) * time_early  
          # cluster_B_late = (cluster==B) * time_late
          # cluster_C_late = (cluster==C) * time_late
          cluster_levels <- NULL
          if (is.factor(surv_gene_td$cluster)) {
            cluster_levels <- levels(surv_gene_td$cluster)
            # Create interactions for non-reference levels (B and C, assuming A is reference)
            for (lev in cluster_levels[-1]) {  # Skip first level (reference)
              surv_gene_td[[paste0("cluster_", lev, "_early")]] <- 
                as.numeric(surv_gene_td$cluster == lev) * surv_gene_td$time_early
              surv_gene_td[[paste0("cluster_", lev, "_late")]] <- 
                as.numeric(surv_gene_td$cluster == lev) * surv_gene_td$time_late
            }
            
            # Build formula with interaction terms
            td_terms <- paste0("cluster_", cluster_levels[-1], "_early", collapse = " + ")
            td_terms <- paste(td_terms, "+", paste0("cluster_", cluster_levels[-1], "_late", collapse = " + "))
            fml_td_text <- paste("Surv(time, event) ~", base_terms, "+", td_terms)
            fml_td <- as.formula(fml_td_text)
            
            td_fit <- try(survival::coxph(fml_td, data = surv_gene_td), silent = TRUE)
            if (!inherits(td_fit, "try-error")) {
              td_model_gene <- summary(td_fit)
            }
          }
        }
      }
      
      p <- draw_km_plot(
        fit = fit_g,
        surv_df = surv_gene,
        cox_unadj = cox_unadj,
        cox_adj = cox_adj,
        cont_unadj = cont_unadj_sum,
        cont_adj = cont_adj_sum,
        main_title = paste(plot_params$plot_title, "-", g),
        local_plot_params = gene_params,
        attempted_adj = (!is.null(clinical_covariates) && length(clinical_covariates) > 0),
        adj_reason = if (is.null(cox_adj) && (!is.null(clinical_covariates) && length(clinical_covariates) > 0)) "coxph error or convergence failure" else NULL
      )
      
      # Store PH test and TD model results with the plot
      if (!is.null(ph_test_gene) && !inherits(ph_test_gene, "try-error")) {
        p$ph_test <- ph_test_gene
      }
      if (!is.null(td_model_gene)) {
        p$td_model <- td_model_gene
      }
      
      results$plots[[g]] <- p
      if (verbose) print(p)
      
      if (save_plots) {
        ggplot2::ggsave(
          file.path(output_dir, paste0("KM_", g, ".png")),
          plot = p$plot, width = 10, height = 8.5, dpi = 300, bg = "white"
        )
      }
    }
    # Combine individual gene plots into multi-panel figures if requested
    if (combine_gene_plots) {
  # Use genes_for_combination if defined (may include all genes when none significant)
  gene_base_set <- if (exists("genes_for_combination")) genes_for_combination else genes_to_plot
  available_gene_plots <- Filter(function(nm) nm %in% gene_base_set, names(results$plots))
      grobs <- lapply(available_gene_plots, function(nm) results$plots[[nm]]$plot)
      if (length(grobs) > 0) {
        # Determine layout: cap columns, compute rows, paginate if needed
        total <- length(grobs)
        per_page <- min(max_plots_per_page, total)
        pages <- ceiling(total / per_page)
        if (verbose) cat("Combining", total, "gene KM plots across", pages, "page(s).\n")
        if (!requireNamespace("cowplot", quietly = TRUE) && !requireNamespace("patchwork", quietly = TRUE) && !requireNamespace("gridExtra", quietly = TRUE)) {
          warning("No layout package (cowplot/patchwork/gridExtra) available; skipping combined gene plot figure.")
        } else {
          for (pg in seq_len(pages)) {
            idx_start <- (pg - 1) * per_page + 1
            idx_end   <- min(pg * per_page, total)
            subset_grobs <- grobs[idx_start:idx_end]
            n_sub <- length(subset_grobs)
            ncol <- min(max_cols_combined, ceiling(sqrt(n_sub)))
            nrow <- ceiling(n_sub / ncol)
            # Adjust columns if too many rows (try to keep grid balanced)
            while ((nrow - ncol) > 1 && ncol < max_cols_combined) {
              ncol <- ncol + 1
              nrow <- ceiling(n_sub / ncol)
            }
            combined_plot <- NULL
            if (requireNamespace("patchwork", quietly = TRUE)) {
              # Use patchwork for simple assembly
              patch <- Reduce(function(a,b) a + b, subset_grobs)
              combined_plot <- patchwork::wrap_plots(subset_grobs, ncol = ncol)
            } else if (requireNamespace("cowplot", quietly = TRUE)) {
              combined_plot <- cowplot::plot_grid(plotlist = subset_grobs, ncol = ncol, align = "hv")
            } else if (requireNamespace("gridExtra", quietly = TRUE)) {
              combined_plot <- gridExtra::arrangeGrob(grobs = subset_grobs, ncol = ncol)
            }
            if (!is.null(combined_plot) && save_plots) {
              fname <- if (pages > 1) sprintf("KM_genes_combined_page%d.png", pg) else "KM_genes_combined.png"
              ggplot2::ggsave(file.path(output_dir, fname), plot = combined_plot, width = 12, height = 4 + 3 * nrow, dpi = 300, bg = "white")
            }
            if (!is.null(combined_plot)) {
              # Store first combined in results, subsequent pages with suffix
              key <- if (pages > 1) sprintf("genes_combined_%d", pg) else "genes_combined"
              results$plots[[key]] <- list(plot = combined_plot)
            }
          }
        }
      }
    }
    grouped_gene_rows <- Filter(Negate(is.null), grouped_gene_rows)
    if (length(grouped_gene_rows)) {
      results$results_grouped_gene <- dplyr::bind_rows(grouped_gene_rows)
    }

    # Optional: per-gene early/late window KM analyses (similar to signature time-dependent windows)
    if (gene_early_late && run_time_dependent) {
      if (verbose) cat("Generating per-gene early/late KM plots...\n")
      for (g in genes_to_plot) {
        # Build survival frame with expression for gene g
        if (!(g %in% rownames(gene_mat))) next
        expr_g <- get_expr_aligned(g, grouped_gene_mat, merged_with_sig, available_samples)
        gene_df <- merged_with_sig %>% dplyr::mutate(expr_val = expr_g)
        # Early window truncation / event logic mirroring signature implementation
        early_gene <- gene_df %>%
          dplyr::mutate(
            time = pmin(time, cut_months),
            event = ifelse(time <= cut_months & event == 1, 1, 0)
          )
        # Late window data (time > cut)
        late_gene <- gene_df %>% dplyr::filter(time > cut_months)

        # Helper to form clusters by expression (reuse make_clusters)
        assign_clusters <- function(df) {
          cl <- make_clusters(df$expr_val, split_method, percentile_cutoffs, include_groups, detection_threshold)
          df$cluster <- cl
          df <- df[!is.na(df$cluster), , drop = FALSE]
          df
        }
        early_gene <- assign_clusters(early_gene)
        late_gene  <- assign_clusters(late_gene)

        # Function to build and store plot if valid
        build_gene_window_plot <- function(df, window_label, suffix) {
          # Require a modest number of samples; allow down to 5 to be more permissive in windows
          if (nrow(df) < 5) return(NULL)
          cl_levels <- nlevels(factor(df$cluster))
          if (cl_levels < 2) {
            # Fallback 1: try high vs low using outer quartiles with cut()
            q25 <- stats::quantile(df$expr_val, 0.25, na.rm = TRUE)
            q75 <- stats::quantile(df$expr_val, 0.75, na.rm = TRUE)
            breaks <- c(-Inf, q25, q75, Inf)
            cl_q <- cut(df$expr_val, breaks = breaks, labels = c("A", "B", "C"),
                       include.lowest = TRUE, right = TRUE)
            cl_q <- as.character(cl_q)
            # Keep only low and high for this fallback
            cl_q[cl_q == "B"] <- NA
            if (nlevels(factor(cl_q)) < 2 || sum(!is.na(cl_q)) < 5) {
              # Fallback 2: median split
              med <- stats::median(df$expr_val, na.rm = TRUE)
              cl_med <- ifelse(df$expr_val <= med, "A", "C")
              if (nlevels(factor(cl_med)) < 2) return(NULL) else df$cluster <- cl_med
            } else {
              df$cluster <- cl_q
            }
          }
          df$cluster <- factor(df$cluster, levels = c("A","B","C"))
          fit <- survival::survfit(survival::Surv(time, event) ~ cluster, data = df)
          cox_unadj <- summary(survival::coxph(survival::Surv(time, event) ~ cluster, data = df))
          # Adjusted model (grouped)
          cox_adj <- NULL
            if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
              keep_covs <- intersect(clinical_covariates, setdiff(colnames(df), c("time","event","cluster","expr_val")))
              if (length(keep_covs)) {
                okc <- vapply(keep_covs, function(cv) { v<-df[[cv]]; if (is.numeric(v)) var(v,na.rm=TRUE)>0 else if (is.factor(v)||is.character(v)) length(unique(na.omit(v)))>1 else FALSE }, logical(1))
                keep_covs <- keep_covs[okc]
                if (length(keep_covs)) {
                  f_adj <- as.formula(paste("survival::Surv(time, event) ~ cluster +", paste(keep_covs, collapse=" + ")))
                  df_adj <- df[, c("time", "event", "cluster", keep_covs), drop = FALSE]
                  df_adj <- df_adj[complete.cases(df_adj), , drop = FALSE]
                  cox_adj_try <- if (nrow(df_adj) >= 10) try(survival::coxph(f_adj, data = df_adj), silent = TRUE) else NULL
                  if (!inherits(cox_adj_try, "try-error") && !is.null(cox_adj_try)) {
                    cox_adj <- summary(cox_adj_try)
                  }
                }
              }
            }
          # Continuous models for gene expression
          cont_unadj_fit <- try(survival::coxph(survival::Surv(time, event) ~ expr_val, data = df), silent = TRUE)
          cont_unadj <- if (!inherits(cont_unadj_fit, "try-error")) summary(cont_unadj_fit) else NULL
          cont_adj <- NULL
          if (!is.null(cox_adj)) {
            if (!is.null(clinical_covariates) && length(clinical_covariates) > 0) {
              keep_covs2 <- intersect(clinical_covariates, setdiff(colnames(df), c("time","event","expr_val")))
              if (length(keep_covs2)) {
                okc2 <- vapply(keep_covs2, function(cv) { v<-df[[cv]]; if (is.numeric(v)) var(v,na.rm=TRUE)>0 else if (is.factor(v)||is.character(v)) length(unique(na.omit(v)))>1 else FALSE }, logical(1))
                keep_covs2 <- keep_covs2[okc2]
                if (length(keep_covs2)) {
                  f_cont_adj <- as.formula(paste("survival::Surv(time, event) ~ expr_val +", paste(keep_covs2, collapse=" + ")))
                  df_cont_adj <- df[, c("time", "event", "expr_val", keep_covs2), drop = FALSE]
                  df_cont_adj <- df_cont_adj[complete.cases(df_cont_adj), , drop = FALSE]
                  cont_adj_fit <- if (nrow(df_cont_adj) >= 10) try(survival::coxph(f_cont_adj, data = df_cont_adj), silent = TRUE) else NULL
                  if (!inherits(cont_adj_fit, "try-error") && !is.null(cont_adj_fit)) {
                    cont_adj <- summary(cont_adj_fit)
                  }
                }
              }
            }
          }
          local_params <- plot_params
          local_params$stats_panel_mode <- "bottom_table"
          local_params$title_size <- max(10, (plot_params$title_size %||% 15) * 0.85)
          main_title <- wrap_title(paste0(window_label, " (", g, ")"))
          p_obj <- draw_km_plot(
            fit = fit,
            surv_df = df %>% dplyr::mutate(cluster = factor(cluster, levels = c("A","B","C"))),
            cox_unadj = cox_unadj,
            cox_adj = cox_adj,
            cont_unadj = cont_unadj,
            cont_adj = cont_adj,
            main_title = main_title,
            local_plot_params = local_params,
            attempted_adj = (!is.null(clinical_covariates) && length(clinical_covariates) > 0),
            adj_reason = if (is.null(cox_adj) && (!is.null(clinical_covariates) && length(clinical_covariates) > 0)) "coxph error or convergence failure" else NULL
          )
          plot_key <- paste0(g, suffix)
          results$plots[[plot_key]] <<- p_obj
          if (save_plots) {
            ggplot2::ggsave(file.path(output_dir, paste0("KM_", g, suffix, ".png")), plot = p_obj$plot, width = 8.5, height = 7.5, dpi = 300, bg = "white")
          }
        }
        build_gene_window_plot(early_gene, sprintf("%s - Early Window ≤%d m", plot_params$plot_title, cut_months), "_early")
        build_gene_window_plot(late_gene,  sprintf("%s - Late Window >%d m",  plot_params$plot_title, cut_months), "_late")
      }
      # After generating individual early/late plots, build combined panels
      if (combine_gene_plots) {
        if (!requireNamespace("cowplot", quietly = TRUE) && !requireNamespace("patchwork", quietly = TRUE) && !requireNamespace("gridExtra", quietly = TRUE)) {
          if (verbose) warning("No layout package for combined early/late gene panels.")
        } else {
          # Helper to assemble combined early or late plots
          assemble_window <- function(suffix, tag) {
            base_set <- if (exists("genes_for_combination")) genes_for_combination else genes_to_plot
            candidate_names <- paste0(base_set, suffix)
            existing <- intersect(candidate_names, names(results$plots))
            # Exclude any signature-level plots (now named signature_early / signature_late)
            existing <- setdiff(existing, c("signature_early", "signature_late"))
            if (!length(existing)) {
              # Fallback: look for already-saved PNG files (KM_<GENE><suffix>.png) if in-memory plots absent
              if (verbose) cat("No in-memory", tag, "plots found; searching files for fallback (suffix", suffix, ")\n")
              suffix_clean <- sub("^_", "", suffix)  # e.g. _early -> early
              # Build regex pattern to match any of the base_set genes
              gene_regex <- paste(base_set, collapse = "|")
              file_pattern <- paste0("^KM_(", gene_regex, ")", suffix, "\\.png$")
              png_files <- list.files(output_dir, pattern = file_pattern, full.names = TRUE)
              if (length(png_files) == 0) {
                if (verbose) cat("No", tag, "PNG files found for fallback combination.\n")
                return(NULL)
              }
              if (!requireNamespace("png", quietly = TRUE) || !requireNamespace("grid", quietly = TRUE)) {
                warning("Cannot combine file PNGs for ", tag, ": 'png' and 'grid' packages required.")
                return(NULL)
              }
              # Construct pseudo-ggplot objects from PNG rasters
              grobs <- lapply(png_files, function(fp) {
                img <- png::readPNG(fp)
                g <- grid::rasterGrob(img, interpolate = TRUE)
                # Derive a clean title from filename like KM_GENE_suffix.png -> GENE_suffix
                base_name <- basename(fp)
                title_txt <- sub("^KM_", "", base_name)
                title_txt <- sub("\\.png$", "", title_txt)
                ggplot2::ggplot() +
                  ggplot2::annotation_custom(g) +
                  ggplot2::theme_void() +
                  ggplot2::labs(title = title_txt)
              })
              # Proceed to layout & save combined from file-based grobs
              total <- length(grobs)
              per_page <- min(max_plots_per_page, total)
              pages <- ceiling(total / per_page)
              if (verbose) cat("Combining", total, tag, "gene KM PNG fallback plots across", pages, "page(s). (suffix=", suffix, ")\n", sep=" ")
              for (pg in seq_len(pages)) {
                idx_start <- (pg - 1) * per_page + 1
                idx_end <- min(pg * per_page, total)
                subset_grobs <- grobs[idx_start:idx_end]
                n_sub <- length(subset_grobs)
                ncol <- min(max_cols_combined, ceiling(sqrt(n_sub)))
                nrow <- ceiling(n_sub / ncol)
                while ((nrow - ncol) > 1 && ncol < max_cols_combined) {
                  ncol <- ncol + 1; nrow <- ceiling(n_sub / ncol)
                }
                combined_plot <- NULL
                if (requireNamespace("patchwork", quietly = TRUE)) {
                  combined_plot <- patchwork::wrap_plots(subset_grobs, ncol = ncol)
                } else if (requireNamespace("cowplot", quietly = TRUE)) {
                  combined_plot <- cowplot::plot_grid(plotlist = subset_grobs, ncol = ncol, align = "hv")
                } else if (requireNamespace("gridExtra", quietly = TRUE)) {
                  combined_plot <- gridExtra::arrangeGrob(grobs = subset_grobs, ncol = ncol)
                }
                if (!is.null(combined_plot) && save_plots) {
                  fname <- if (pages > 1) sprintf("KM_genes_%s_combined_page%d.png", tag, pg) else sprintf("KM_genes_%s_combined.png", tag)
                  ggplot2::ggsave(file.path(output_dir, fname), plot = combined_plot, width = 12, height = 4 + 3 * nrow, dpi = 300, bg = "white")
                  if (verbose) cat("Saved", fname, "(fallback) to", output_dir, "\n")
                }
                if (!is.null(combined_plot)) {
                  key <- if (pages > 1) sprintf("genes_%s_combined_%d", tag, pg) else sprintf("genes_%s_combined", tag)
                  results$plots[[key]] <- list(plot = combined_plot)
                }
              }
              return(invisible(NULL))
            }
            grobs <- lapply(existing, function(nm) results$plots[[nm]]$plot)
            total <- length(grobs)
            per_page <- min(max_plots_per_page, total)
            pages <- ceiling(total / per_page)
            if (verbose) cat("Combining", total, tag, "gene KM plots across", pages, "page(s). (suffix=", suffix, ")\n", sep=" ")
            for (pg in seq_len(pages)) {
              idx_start <- (pg - 1) * per_page + 1
              idx_end <- min(pg * per_page, total)
              subset_grobs <- grobs[idx_start:idx_end]
              n_sub <- length(subset_grobs)
              ncol <- min(max_cols_combined, ceiling(sqrt(n_sub)))
              nrow <- ceiling(n_sub / ncol)
              while ((nrow - ncol) > 1 && ncol < max_cols_combined) {
                ncol <- ncol + 1
                nrow <- ceiling(n_sub / ncol)
              }
              combined_plot <- NULL
              if (requireNamespace("patchwork", quietly = TRUE)) {
                combined_plot <- patchwork::wrap_plots(subset_grobs, ncol = ncol)
              } else if (requireNamespace("cowplot", quietly = TRUE)) {
                combined_plot <- cowplot::plot_grid(plotlist = subset_grobs, ncol = ncol, align = "hv")
              } else if (requireNamespace("gridExtra", quietly = TRUE)) {
                combined_plot <- gridExtra::arrangeGrob(grobs = subset_grobs, ncol = ncol)
              }
              if (!is.null(combined_plot) && save_plots) {
                fname <- if (pages > 1) sprintf("KM_genes_%s_combined_page%d.png", tag, pg) else sprintf("KM_genes_%s_combined.png", tag)
                ggplot2::ggsave(file.path(output_dir, fname), plot = combined_plot, width = 12, height = 4 + 3 * nrow, dpi = 300, bg = "white")
                if (verbose) cat("Saved", fname, "to", output_dir, "\n")
              }
              if (!is.null(combined_plot)) {
                key <- if (pages > 1) sprintf("genes_%s_combined_%d", tag, pg) else sprintf("genes_%s_combined", tag)
                results$plots[[key]] <- list(plot = combined_plot)
              }
            }
          }
          assemble_window("_early", "early")
          assemble_window("_late",  "late")
          if (verbose) cat("Finished assembling early and late combined gene KM panels.\n")
          # (Removed early+late single figure per user clarification; keep separate early and late combined only)
        }
      }
    }
  }
  
  # ============================================================================
  # BUILD SAMPLE GROUP MEMBERSHIPS
  # ============================================================================
  
  if (verbose) cat("\n=== Building Sample Group Memberships ===\n")
  
  group_list <- list()
  
  # Signature groups
  if (exists("surv_data", inherits = FALSE) && nrow(surv_data) > 0) {
    sig_members <- split(surv_data$sample_id, surv_data$cluster)
    for (grp in names(sig_members)) {
      col_name <- paste0("Signature_", 
                        if (grp == "A") "Low" else if (split_method == "median" && grp == "B") "High" else if (grp == "B") "Mid" else "High")
      group_list[[col_name]] <- sig_members[[grp]]
    }
  }
  
  # Per-gene groups
  for (g in genes) {
    expr_g <- get_expr_aligned(g, grouped_gene_mat, merged_base, available_samples)
    cl <- make_clusters(expr_g, split_method, percentile_cutoffs, include_groups, detection_threshold)
    if (!all(is.na(cl))) {
      sids_by_group <- split(merged_base$sample_id[!is.na(cl)], cl[!is.na(cl)])
      for (lev in names(sids_by_group)) {
        col_name <- paste0(
          g,
          "_",
          if (split_method == "detectable") {
            if (lev == "A") "Undetectable" else "Detectable"
          } else if (lev == "A") {
            "Low"
          } else if (split_method == "median" && lev == "B") {
            "High"
          } else if (lev == "B") {
            "Mid"
          } else {
            "High"
          }
        )
        group_list[[col_name]] <- sids_by_group[[lev]]
      }
    }
  }
  
  # Pad to equal length
  if (length(group_list) > 0) {
    maxlen <- max(vapply(group_list, length, 1L))
    group_list <- lapply(group_list, function(x) {
      c(x, rep(NA_character_, maxlen - length(x)))
    })
    results$sample_groups <- as.data.frame(group_list, stringsAsFactors = FALSE)
  }
  
  # ============================================================================
  # SUMMARY
  # ============================================================================
  
  if (verbose) {
    cat("\n", paste(rep("=", 80), collapse = ""), "\n")
    cat("ANALYSIS COMPLETE\n")
    cat(paste(rep("=", 80), collapse = ""), "\n")
    
    if (!is.null(results$results_adjusted)) {
      sig_genes <- results$results_adjusted$Gene[results$results_adjusted$FDR < fdr_threshold]
      cat(sprintf("Significant genes (FDR < %.2f): %d\n", fdr_threshold, length(sig_genes)))
      if (length(sig_genes) > 0) {
        cat("Genes:", paste(sig_genes, collapse = ", "), "\n")
      }
    }
    
    if (!is.null(results$lasso_selected) && length(results$lasso_selected) > 0) {
      cat(sprintf("LASSO selected: %d genes\n", length(results$lasso_selected)))
    }
    
    cat(sprintf("Total samples analyzed: %d\n", nrow(merged_base)))
    cat(sprintf("Genes in signature: %d\n", length(genes)))
  }
  
  return(results)
}
