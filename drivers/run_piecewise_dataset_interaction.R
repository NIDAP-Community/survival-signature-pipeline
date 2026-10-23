# run_piecewise_dataset_interaction.R
# Piecewise (0-60 vs >60 months) Cox models with interaction across two datasets.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
})

source("src/utilities/cli_helpers.R")
source("src/utilities/gene_identifier_resolution.R")
source("src/utilities/run_provenance_guardrails.R")

cli_defaults <- list(
  output_dir = default_run_output_dir("piecewise_dataset_interaction"),
  dataset1_expr_file = file.path("test_data", "tcga_expression.tsv"),
  dataset1_meta_file = file.path("test_data", "tcga_metadata.csv"),
  dataset2_expr_file = file.path("test_data", "metabric_expression.tsv"),
  dataset2_meta_file = file.path("test_data", "metabric_metadata.csv"),
  dataset1_label = "Dataset1",
  dataset2_label = "Dataset2",
  dataset1_sample_column = "sampnames",
  dataset2_sample_column = "sampnames",
  dataset1_subtype_column = "Subtype",
  dataset2_subtype_column = "Subtype",
  dataset1_stage_column = "ajcc_pathologic_tumor_stage",
  dataset2_stage_column = "stage_collapsed",
  dataset1_age_column = "Diagnosis_Age",
  dataset2_age_column = "Diagnosis_Age",
  gene_column = "Gene",
  gene_resolution = "none",
  signature_genes = "GENE_A,GENE_B,GENE_C,GENE_D",
  os_event_column = "OS",
  os_time_column = "OS_time",
  pfi_event_column = "PFI",
  pfi_time_column = "PFI_time",
  subtype = "ALL",
  adjust_covariates = "none",
  features = "GENE_A,GENE_B,GENE_C,GENE_D,Signature",
  endpoints = "OS,PFI",
  cut_months = "60"
)

print_usage <- function() {
  print_standard_cli_help(
    title = "Piecewise Dataset Interaction Driver",
    workflow_name = "piecewise-dataset-interaction",
    script_path = file.path("drivers", "run_piecewise_dataset_interaction.R"),
    defaults = cli_defaults,
    summary_lines = c(
      "This workflow expects two harmonized datasets and models dataset interaction terms across them.",
      "Endpoint columns are configurable through the OS/PFI event and time flags.",
      "Inputs are expected to be post-filtered upstream; the maintained driver does not apply an internal primary-only sample filter.",
      "gene_resolution=l2p opt-in updates user-supplied human gene symbols and requires a matched row in both datasets.",
      "dataset1_label and dataset2_label control the dataset titles written into pooled outputs."
    ),
    aliases = c(
      tcga_expr_file = "dataset1_expr_file",
      tcga_meta_file = "dataset1_meta_file",
      metabric_expr_file = "dataset2_expr_file",
      metabric_meta_file = "dataset2_meta_file",
      tcga_sample_column = "dataset1_sample_column",
      metabric_sample_column = "dataset2_sample_column",
      tcga_subtype_column = "dataset1_subtype_column",
      metabric_subtype_column = "dataset2_subtype_column",
      tcga_stage_column = "dataset1_stage_column",
      metabric_stage_column = "dataset2_stage_column",
      tcga_age_column = "dataset1_age_column",
      metabric_age_column = "dataset2_age_column"
    )
  )
}

opt <- parse_standard_cli_args(
  defaults = cli_defaults,
  workflow_name = "piecewise-dataset-interaction",
  help_printer = print_usage,
  aliases = c(
    tcga_expr_file = "dataset1_expr_file",
    tcga_meta_file = "dataset1_meta_file",
    metabric_expr_file = "dataset2_expr_file",
    metabric_meta_file = "dataset2_meta_file",
    tcga_sample_column = "dataset1_sample_column",
    metabric_sample_column = "dataset2_sample_column",
    tcga_subtype_column = "dataset1_subtype_column",
    metabric_subtype_column = "dataset2_subtype_column",
    tcga_stage_column = "dataset1_stage_column",
    metabric_stage_column = "dataset2_stage_column",
    tcga_age_column = "dataset1_age_column",
    metabric_age_column = "dataset2_age_column"
  )
)
opt$cut_months <- cli_parse_numeric_flag(opt$cut_months, "cut_months")
opt$adjust_covariates <- cli_parse_enum_flag(
  value = opt$adjust_covariates,
  allowed_values = c("none", "subtype_stage", "subtype_stage_age"),
  flag_name = "adjust_covariates",
  normalize = "lower"
)
opt$gene_resolution <- normalize_gene_resolution_mode(opt$gene_resolution)
feature_list <- cli_parse_csv_list(opt$features)
endpoint_list <- cli_parse_csv_list(opt$endpoints, uppercase = TRUE)
signature_genes <- cli_parse_csv_list(opt$signature_genes)

if (!nzchar(opt$dataset1_label) || !nzchar(opt$dataset2_label)) {
  stop("--dataset1-label and --dataset2-label must both be non-empty.")
}
if (identical(opt$dataset1_label, opt$dataset2_label)) {
  stop("--dataset1-label and --dataset2-label must be distinct.")
}

assert_output_guardrails(opt$output_dir)

cli_assert_existing_input_files(
  input_paths = c(
    dataset1_expr_file = opt$dataset1_expr_file,
    dataset1_meta_file = opt$dataset1_meta_file,
    dataset2_expr_file = opt$dataset2_expr_file,
    dataset2_meta_file = opt$dataset2_meta_file
  ),
  workflow_name = "piecewise-dataset-interaction"
)

parse_time_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x0 <- as.character(x)
  x0[x0 %in% c("", "NA", "NaN", "null", "_Not_Available_", "_Not_Applicable_", "[Not Available]", "[Not Applicable]")] <- NA
  x0 <- trimws(gsub(",", "", x0))
  x0 <- gsub("[^0-9eE+\\-\\.]", "", x0)
  suppressWarnings(as.numeric(x0))
}

parse_event_numeric <- function(x) {
  if (is.numeric(x)) return(x)
  x0 <- tolower(as.character(x))
  ifelse(x0 == "dead" | x0 == "deceased", 1,
         ifelse(x0 == "alive" | x0 == "living", 0,
                suppressWarnings(as.numeric(x0))))
}

make_clusters <- function(x_num) {
  tertile_1 <- stats::quantile(x_num, 1 / 3, na.rm = TRUE)
  tertile_2 <- stats::quantile(x_num, 2 / 3, na.rm = TRUE)
  breaks <- c(-Inf, tertile_1, tertile_2, Inf)
  cluster_values <- cut(x_num, breaks = breaks, labels = c("A", "B", "C"), include.lowest = TRUE, right = TRUE)
  cluster_values <- as.character(cluster_values)
  cluster_values[cluster_values == "A"] <- "Low"
  cluster_values[cluster_values == "C"] <- "High"
  cluster_values[!(cluster_values %in% c("Low", "High"))] <- NA_character_
  factor(cluster_values, levels = c("Low", "High"))
}

compute_feature_score <- function(expr, gene_column, sample_ids, feature_name, signature_genes) {
  if (tolower(feature_name) == "signature") {
    expression_subset <- expr %>% filter(.data[[gene_column]] %in% signature_genes)
    if (!nrow(expression_subset)) stop("None of the signature genes found in expression matrix.")
    matrix_data <- as.matrix(expression_subset[, setdiff(colnames(expression_subset), gene_column), drop = FALSE])
    rownames(matrix_data) <- expression_subset[[gene_column]]
    matrix_data <- matrix_data[, sample_ids, drop = FALSE]
    matrix_data <- apply(matrix_data, 2, as.numeric)
    if (is.null(dim(matrix_data))) matrix_data <- matrix(matrix_data, nrow = nrow(expression_subset))
    rownames(matrix_data) <- expression_subset[[gene_column]]
    colnames(matrix_data) <- sample_ids
    matrix_z <- t(scale(t(matrix_data)))
    score <- colSums(matrix_z, na.rm = TRUE)
    return(score)
  }

  expression_row <- expr %>% filter(.data[[gene_column]] == feature_name)
  if (!nrow(expression_row)) stop("Feature not found in expression matrix: ", feature_name)
  value_vector <- as.numeric(expression_row[1, setdiff(colnames(expression_row), gene_column), drop = TRUE])
  names(value_vector) <- setdiff(colnames(expression_row), gene_column)
  value_vector <- value_vector[sample_ids]
  as.numeric(scale(value_vector))
}

harmonize_stage <- function(x) {
  x_char <- tolower(trimws(as.character(x)))
  x_char[x_char %in% c("", "na", "nan", "[not available]", "[not applicable]", "not reported", "other/unknown")] <- NA

  early <- !is.na(x_char) & (
    grepl("early", x_char) |
      x_char %in% c("i", "ia", "ib", "ic") |
      grepl("^stage\\s*i([a-c])?$", x_char)
  )

  late <- !is.na(x_char) & (
    grepl("late", x_char) |
      x_char %in% c("ii", "iia", "iib", "iic", "iii", "iiia", "iiib", "iiic", "iv") |
      grepl("^stage\\s*ii([a-c])?$", x_char) |
      grepl("^stage\\s*iii([a-c])?$", x_char) |
      grepl("^stage\\s*iv([a-c])?$", x_char)
  )

  out <- rep(NA_character_, length(x_char))
  out[early] <- "Early"
  out[late] <- "Late"
  out
}

prepare_dataset <- function(dataset_label, expr, metadata, gene_column, sample_column, subtype_column, stage_column, age_column, subtype_filter, endpoint, feature_name, signature_genes, options) {
  if (!sample_column %in% names(metadata)) stop(dataset_label, " metadata missing sample column: ", sample_column)
  if (!subtype_column %in% names(metadata)) stop(dataset_label, " metadata missing subtype column: ", subtype_column)
  if (!stage_column %in% names(metadata)) stop(dataset_label, " metadata missing stage column: ", stage_column)
  if (!age_column %in% names(metadata)) stop(dataset_label, " metadata missing age column: ", age_column)

  endpoint_columns <- resolve_standard_endpoint_columns(endpoint, options, "piecewise-dataset-interaction")
  event_column <- endpoint_columns$event
  time_column <- endpoint_columns$time

  metadata_filtered <- metadata
  subtype_clean <- tolower(trimws(subtype_filter))
  if (!(subtype_clean %in% c("", "all", "overall", "any"))) {
    selected_subtypes <- cli_parse_csv_list(subtype_filter)
    metadata_filtered <- metadata_filtered %>% filter(.data[[subtype_column]] %in% selected_subtypes)
  }


  if (!all(c(event_column, time_column) %in% names(metadata_filtered))) {
    stop("Missing endpoint columns in ", dataset_label, ": ", event_column, ", ", time_column)
  }

  expression_samples <- setdiff(colnames(expr), gene_column)
  metadata_filtered <- metadata_filtered %>% filter(.data[[sample_column]] %in% expression_samples)
  metadata_filtered <- metadata_filtered %>%
    mutate(
      .time_raw = parse_time_numeric(.data[[time_column]]),
      event = parse_event_numeric(.data[[event_column]])
    ) %>%
    filter(!is.na(.time_raw) & !is.na(event))

  if (!nrow(metadata_filtered)) stop("No analyzable rows after filtering for ", dataset_label, ", ", feature_name, ", ", endpoint)

  median_time <- median(metadata_filtered$.time_raw, na.rm = TRUE)
  max_time <- max(metadata_filtered$.time_raw, na.rm = TRUE)
  in_days <- is.finite(max_time) && (max_time > 1200 || median_time > 365)
  metadata_filtered <- metadata_filtered %>% mutate(time = if (in_days) .time_raw / 30.44 else .time_raw)

  score <- compute_feature_score(expr, gene_column, metadata_filtered[[sample_column]], feature_name, signature_genes)
  analysis_df <- metadata_filtered %>%
    transmute(
      dataset = dataset_label,
      sample_id = .data[[sample_column]],
      time = time,
      event = event,
      subtype_cov = as.character(.data[[subtype_column]]),
      stage_cov = harmonize_stage(.data[[stage_column]]),
      age_cov = parse_time_numeric(.data[[age_column]]),
      feature_score = score
    ) %>%
    filter(!is.na(feature_score))

  analysis_df$cluster <- make_clusters(analysis_df$feature_score)
  analysis_df <- analysis_df %>% filter(!is.na(cluster))
  analysis_df$group_high <- ifelse(analysis_df$cluster == "High", 1, 0)

  analysis_df
}

extract_term <- function(fit, term_name) {
  fit_summary <- summary(fit)
  if (!term_name %in% rownames(fit_summary$coefficients)) {
    return(data.frame(term = term_name, beta = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_))
  }
  coefficient <- fit_summary$coefficients[term_name, "coef"]
  p_value <- fit_summary$coefficients[term_name, "Pr(>|z|)"]
  conf_int <- suppressWarnings(confint(fit))
  ci_low <- exp(conf_int[term_name, 1])
  ci_high <- exp(conf_int[term_name, 2])

  data.frame(
    term = term_name,
    beta = coefficient,
    HR = exp(coefficient),
    CI_low = ci_low,
    CI_high = ci_high,
    p_value = p_value
  )
}

extract_lrt_p <- function(anova_result) {
  p_column <- grep("Pr", colnames(anova_result), value = TRUE)
  if (!length(p_column)) return(NA_real_)
  as.numeric(anova_result[nrow(anova_result), p_column[1]])
}

extract_all_terms <- function(fit) {
  fit_summary <- summary(fit)
  conf_int <- suppressWarnings(confint(fit))
  term_names <- rownames(fit_summary$coefficients)

  out <- data.frame(
    term = term_names,
    beta = fit_summary$coefficients[, "coef"],
    HR = fit_summary$conf.int[, "exp(coef)"],
    CI_low = NA_real_,
    CI_high = NA_real_,
    p_value = fit_summary$coefficients[, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )

  common_terms <- intersect(term_names, rownames(conf_int))
  out$CI_low[match(common_terms, out$term)] <- exp(conf_int[common_terms, 1])
  out$CI_high[match(common_terms, out$term)] <- exp(conf_int[common_terms, 2])
  out
}

resolve_dataset_interaction_terms <- function(fit) {
  coefficient_names <- names(stats::coef(fit))
  early_term <- grep("^group_high:dataset[^:]+$", coefficient_names, value = TRUE)
  late_term <- grep("^group_high:dataset[^:]+:post60$|^group_high:post60:dataset[^:]+$", coefficient_names, value = TRUE)

  if (length(early_term) > 1 || length(late_term) > 1) {
    stop("Unable to resolve unique dataset interaction terms from fitted model coefficients.")
  }

  list(
    early = if (length(early_term)) early_term[[1]] else NA_character_,
    late = if (length(late_term)) late_term[[1]] else NA_character_
  )
}

extract_contrast <- function(fit, term_weights) {
  coefficients <- stats::coef(fit)
  variance_matrix <- stats::vcov(fit)
  required_terms <- names(term_weights)
  if (!all(required_terms %in% names(coefficients))) {
    return(data.frame(beta = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_))
  }

  beta_values <- coefficients[required_terms]
  variance_submatrix <- variance_matrix[required_terms, required_terms, drop = FALSE]
  weight_vector <- as.numeric(term_weights)

  if (any(is.na(beta_values)) || any(is.na(variance_submatrix))) {
    return(data.frame(beta = NA_real_, HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_))
  }

  beta_value <- sum(weight_vector * beta_values)
  se_value <- sqrt(as.numeric(t(weight_vector) %*% variance_submatrix %*% weight_vector))
  if (!is.finite(se_value) || se_value <= 0) {
    return(data.frame(beta = beta_value, HR = exp(beta_value), CI_low = NA_real_, CI_high = NA_real_, p_value = NA_real_))
  }

  z_value <- beta_value / se_value
  p_value <- 2 * stats::pnorm(abs(z_value), lower.tail = FALSE)

  data.frame(
    beta = beta_value,
    HR = exp(beta_value),
    CI_low = exp(beta_value - 1.96 * se_value),
    CI_high = exp(beta_value + 1.96 * se_value),
    p_value = p_value
  )
}

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
if (!dir.exists(file.path(opt$output_dir, "tables"))) dir.create(file.path(opt$output_dir, "tables"))
if (!dir.exists(file.path(opt$output_dir, "logs"))) dir.create(file.path(opt$output_dir, "logs"))

message("Loading expression and metadata files...")
dataset1_expr <- read_expression_table(opt$dataset1_expr_file, opt$gene_column)
dataset2_expr <- read_expression_table(opt$dataset2_expr_file, opt$gene_column)
dataset1_meta <- read_metadata_table(opt$dataset1_meta_file)
dataset2_meta <- read_metadata_table(opt$dataset2_meta_file)

individual_features <- feature_list[tolower(feature_list) != "signature"]
requested_genes <- unique(c(individual_features, if (any(tolower(feature_list) == "signature")) signature_genes else character(0)))
gene_resolution_table <- resolve_requested_gene_identifiers(
  requested_genes = requested_genes,
  expression_gene_ids = intersect(dataset1_expr[[opt$gene_column]], dataset2_expr[[opt$gene_column]]),
  mode = opt$gene_resolution
)
assert_gene_resolution_is_unambiguous(gene_resolution_table)
feature_list[tolower(feature_list) != "signature"] <- replace_requested_gene_identifiers(individual_features, gene_resolution_table)
if (any(tolower(feature_list) == "signature")) {
  signature_genes <- replace_requested_gene_identifiers(signature_genes, gene_resolution_table)
}
write_gene_resolution_report(gene_resolution_table, opt$output_dir)

result_rows <- list()
term_rows <- list()
implied_rows <- list()
all_term_rows <- list()

for (feature_name in feature_list) {
  for (endpoint_name in endpoint_list) {
    message("Fitting piecewise model for ", feature_name, " (", endpoint_name, ")")

    dataset1_df <- prepare_dataset(
      dataset_label = opt$dataset1_label,
      expr = dataset1_expr,
      metadata = dataset1_meta,
      gene_column = opt$gene_column,
      sample_column = opt$dataset1_sample_column,
      subtype_column = opt$dataset1_subtype_column,
      stage_column = opt$dataset1_stage_column,
      age_column = opt$dataset1_age_column,
      subtype_filter = opt$subtype,
      endpoint = endpoint_name,
      feature_name = feature_name,
      signature_genes = signature_genes,
      options = opt
    )

    dataset2_df <- prepare_dataset(
      dataset_label = opt$dataset2_label,
      expr = dataset2_expr,
      metadata = dataset2_meta,
      gene_column = opt$gene_column,
      sample_column = opt$dataset2_sample_column,
      subtype_column = opt$dataset2_subtype_column,
      stage_column = opt$dataset2_stage_column,
      age_column = opt$dataset2_age_column,
      subtype_filter = opt$subtype,
      endpoint = endpoint_name,
      feature_name = feature_name,
      signature_genes = signature_genes,
      options = opt
    )

    combined_df <- bind_rows(dataset1_df, dataset2_df)
    combined_df$dataset <- factor(combined_df$dataset, levels = c(opt$dataset1_label, opt$dataset2_label))

    early_df <- combined_df %>%
      transmute(
        dataset = dataset,
        sample_id = sample_id,
        group_high = group_high,
        subtype_cov = subtype_cov,
        stage_cov = stage_cov,
        age_cov = age_cov,
        tstart = 0,
        tstop = pmin(time, opt$cut_months),
        event = ifelse(time <= opt$cut_months & event == 1, 1, 0),
        post60 = 0
      )

    late_df <- combined_df %>%
      filter(time > opt$cut_months) %>%
      transmute(
        dataset = dataset,
        sample_id = sample_id,
        group_high = group_high,
        subtype_cov = subtype_cov,
        stage_cov = stage_cov,
        age_cov = age_cov,
        tstart = opt$cut_months,
        tstop = time,
        event = event,
        post60 = 1
      )

    split_df <- bind_rows(early_df, late_df)

    covariate_mode <- opt$adjust_covariates
    if (covariate_mode == "subtype_stage") {
      split_df <- split_df %>%
        filter(!is.na(subtype_cov), !is.na(stage_cov)) %>%
        mutate(
          subtype_cov = factor(subtype_cov),
          stage_cov = factor(stage_cov, levels = c("Early", "Late"))
        )
      if (!nrow(split_df)) stop("No rows left after subtype/stage complete-case filtering for ", feature_name, " ", endpoint_name)
    } else if (covariate_mode == "subtype_stage_age") {
      split_df <- split_df %>%
        filter(!is.na(subtype_cov), !is.na(stage_cov), !is.na(age_cov)) %>%
        mutate(
          subtype_cov = factor(subtype_cov),
          stage_cov = factor(stage_cov, levels = c("Early", "Late")),
          age_cov = as.numeric(age_cov)
        )
      if (!nrow(split_df)) stop("No rows left after subtype/stage/age complete-case filtering for ", feature_name, " ", endpoint_name)
    }

    formula_full <- if (covariate_mode == "subtype_stage") {
      survival::Surv(tstart, tstop, event) ~ group_high * dataset * post60 + subtype_cov + stage_cov
    } else if (covariate_mode == "subtype_stage_age") {
      survival::Surv(tstart, tstop, event) ~ group_high * dataset * post60 + subtype_cov + stage_cov + age_cov
    } else {
      survival::Surv(tstart, tstop, event) ~ group_high * dataset * post60
    }

    formula_no_dataset_interaction <- if (covariate_mode == "subtype_stage") {
      survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:post60 + dataset:post60 + subtype_cov + stage_cov
    } else if (covariate_mode == "subtype_stage_age") {
      survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:post60 + dataset:post60 + subtype_cov + stage_cov + age_cov
    } else {
      survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:post60 + dataset:post60
    }

    formula_no_threeway <- if (covariate_mode == "subtype_stage") {
      survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:dataset + group_high:post60 + dataset:post60 + subtype_cov + stage_cov
    } else if (covariate_mode == "subtype_stage_age") {
      survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:dataset + group_high:post60 + dataset:post60 + subtype_cov + stage_cov + age_cov
    } else {
      survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:dataset + group_high:post60 + dataset:post60
    }

    fit_full <- survival::coxph(
      formula_full,
      data = split_df
    )

    fit_no_dataset_interaction <- survival::coxph(
      formula_no_dataset_interaction,
      data = split_df
    )

    fit_no_threeway <- survival::coxph(
      formula_no_threeway,
      data = split_df
    )

    lrt_dataset <- anova(fit_no_dataset_interaction, fit_full, test = "LRT")
    lrt_threeway <- anova(fit_no_threeway, fit_full, test = "LRT")
    dataset_terms <- resolve_dataset_interaction_terms(fit_full)

    term_table <- bind_rows(
      extract_term(fit_full, "group_high"),
      extract_term(fit_full, dataset_terms$early),
      extract_term(fit_full, "group_high:post60"),
      extract_term(fit_full, dataset_terms$late)
    ) %>%
      mutate(feature = feature_name, endpoint = endpoint_name)

    term_rows[[length(term_rows) + 1]] <- term_table
    all_term_rows[[length(all_term_rows) + 1]] <- extract_all_terms(fit_full) %>%
      mutate(feature = feature_name, endpoint = endpoint_name, covariates = covariate_mode)

    dataset1_counts <- dataset1_df %>% group_by(cluster) %>% summarise(n = n(), events = sum(event == 1), .groups = "drop")
    dataset2_counts <- dataset2_df %>% group_by(cluster) %>% summarise(n = n(), events = sum(event == 1), .groups = "drop")

    pull_count <- function(df_counts, cluster_name, value_name) {
      value <- df_counts %>% filter(cluster == cluster_name) %>% pull(.data[[value_name]])
      if (!length(value)) return(NA_real_)
      value[1]
    }

    result_rows[[length(result_rows) + 1]] <- data.frame(
      feature = feature_name,
      endpoint = endpoint_name,
      model_covariates = covariate_mode,
      cut_months = opt$cut_months,
      model_rows = nrow(split_df),
      model_events = sum(split_df$event == 1, na.rm = TRUE),
      dataset1_low_n = pull_count(dataset1_counts, "Low", "n"),
      dataset1_high_n = pull_count(dataset1_counts, "High", "n"),
      dataset1_low_events = pull_count(dataset1_counts, "Low", "events"),
      dataset1_high_events = pull_count(dataset1_counts, "High", "events"),
      dataset2_low_n = pull_count(dataset2_counts, "Low", "n"),
      dataset2_high_n = pull_count(dataset2_counts, "High", "n"),
      dataset2_low_events = pull_count(dataset2_counts, "Low", "events"),
      dataset2_high_events = pull_count(dataset2_counts, "High", "events"),
      hr_high_dataset1_early = extract_term(fit_full, "group_high")$HR,
      p_high_dataset1_early = extract_term(fit_full, "group_high")$p_value,
      hr_dataset_interaction_early = extract_term(fit_full, dataset_terms$early)$HR,
      p_dataset_interaction_early = extract_term(fit_full, dataset_terms$early)$p_value,
      hr_dataset1_late_shift = extract_term(fit_full, "group_high:post60")$HR,
      p_dataset1_late_shift = extract_term(fit_full, "group_high:post60")$p_value,
      hr_dataset_late_shift_interaction = extract_term(fit_full, dataset_terms$late)$HR,
      p_dataset_late_shift_interaction = extract_term(fit_full, dataset_terms$late)$p_value,
      lrt_dataset_interaction_p = extract_lrt_p(lrt_dataset),
      lrt_threeway_p = extract_lrt_p(lrt_threeway)
    )

    dataset1_early <- extract_contrast(fit_full, c("group_high" = 1))
    dataset1_late <- extract_contrast(fit_full, c("group_high" = 1, "group_high:post60" = 1))
    dataset2_early <- extract_contrast(fit_full, setNames(c(1, 1), c("group_high", dataset_terms$early)))
    dataset2_late <- extract_contrast(
      fit_full,
      setNames(c(1, 1, 1, 1), c("group_high", dataset_terms$early, "group_high:post60", dataset_terms$late))
    )

    implied_rows[[length(implied_rows) + 1]] <- bind_rows(
      dataset1_early %>% mutate(feature = feature_name, endpoint = endpoint_name, dataset = opt$dataset1_label, window = "0-60m"),
      dataset1_late %>% mutate(feature = feature_name, endpoint = endpoint_name, dataset = opt$dataset1_label, window = ">60m"),
      dataset2_early %>% mutate(feature = feature_name, endpoint = endpoint_name, dataset = opt$dataset2_label, window = "0-60m"),
      dataset2_late %>% mutate(feature = feature_name, endpoint = endpoint_name, dataset = opt$dataset2_label, window = ">60m")
    )
  }
}

results_df <- bind_rows(result_rows)
terms_df <- bind_rows(term_rows)
implied_df <- bind_rows(implied_rows)
all_terms_df <- bind_rows(all_term_rows)

results_df <- results_df %>%
  group_by(endpoint) %>%
  mutate(
    q_dataset_interaction_early_bh = p.adjust(p_dataset_interaction_early, method = "BH"),
    q_dataset_late_shift_interaction_bh = p.adjust(p_dataset_late_shift_interaction, method = "BH"),
    q_lrt_dataset_interaction_bh = p.adjust(lrt_dataset_interaction_p, method = "BH"),
    q_lrt_threeway_bh = p.adjust(lrt_threeway_p, method = "BH")
  ) %>%
  ungroup()

terms_df <- terms_df %>%
  group_by(endpoint, term) %>%
  mutate(p_value_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup()

all_terms_df <- all_terms_df %>%
  group_by(endpoint, term) %>%
  mutate(p_value_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup()

implied_df <- implied_df %>%
  group_by(endpoint, window) %>%
  mutate(p_value_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup()

readr::write_csv(results_df, file.path(opt$output_dir, "tables", "piecewise_dataset_interaction_summary.csv"))
readr::write_csv(terms_df, file.path(opt$output_dir, "tables", "piecewise_dataset_interaction_terms.csv"))
readr::write_csv(all_terms_df, file.path(opt$output_dir, "tables", "piecewise_dataset_interaction_all_terms.csv"))
readr::write_csv(implied_df, file.path(opt$output_dir, "tables", "implied_high_vs_low_by_dataset_window.csv"))
readr::write_csv(
  implied_df %>% filter(feature == "Signature"),
  file.path(opt$output_dir, "tables", "signature_implied_high_vs_low_by_dataset_window.csv")
)

format_num <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = digits))

md_header <- c(
  "| Feature | Endpoint | Covariates | Model rows | Model events | Dataset1 early HR (High vs Low) | Dataset1 early p | Dataset interaction HR (early) | Dataset interaction p (early) | Dataset interaction q (BH) | Dataset1 late-shift HR | Dataset1 late-shift p | Dataset late-shift interaction HR | Dataset late-shift interaction p | Dataset late-shift interaction q (BH) | LRT dataset interaction p | LRT dataset interaction q (BH) | LRT three-way p | LRT three-way q (BH) |",
  "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
)

md_rows <- apply(results_df, 1, function(row_values) {
  paste0(
    "| ", row_values[["feature"]],
    " | ", row_values[["endpoint"]],
    " | ", row_values[["model_covariates"]],
    " | ", format_num(as.numeric(row_values[["model_rows"]]), 4),
    " | ", format_num(as.numeric(row_values[["model_events"]]), 4),
    " | ", format_num(as.numeric(row_values[["hr_high_dataset1_early"]]), 4),
    " | ", format_num(as.numeric(row_values[["p_high_dataset1_early"]]), 3),
    " | ", format_num(as.numeric(row_values[["hr_dataset_interaction_early"]]), 4),
    " | ", format_num(as.numeric(row_values[["p_dataset_interaction_early"]]), 3),
    " | ", format_num(as.numeric(row_values[["q_dataset_interaction_early_bh"]]), 3),
    " | ", format_num(as.numeric(row_values[["hr_dataset1_late_shift"]]), 4),
    " | ", format_num(as.numeric(row_values[["p_dataset1_late_shift"]]), 3),
    " | ", format_num(as.numeric(row_values[["hr_dataset_late_shift_interaction"]]), 4),
    " | ", format_num(as.numeric(row_values[["p_dataset_late_shift_interaction"]]), 3),
    " | ", format_num(as.numeric(row_values[["q_dataset_late_shift_interaction_bh"]]), 3),
    " | ", format_num(as.numeric(row_values[["lrt_dataset_interaction_p"]]), 3),
    " | ", format_num(as.numeric(row_values[["q_lrt_dataset_interaction_bh"]]), 3),
    " | ", format_num(as.numeric(row_values[["lrt_threeway_p"]]), 3),
    " | ", format_num(as.numeric(row_values[["q_lrt_threeway_bh"]]), 3),
    " |"
  )
})

writeLines(c(md_header, md_rows), con = file.path(opt$output_dir, "tables", "piecewise_dataset_interaction_summary.md"))

implied_md <- implied_df %>%
  mutate(
    `HR (High vs Low)` = format_num(HR, 4),
    `95% CI` = paste0("[", format_num(CI_low, 4), ", ", format_num(CI_high, 4), "]"),
    `p` = format_num(p_value, 3),
    `q (BH)` = format_num(p_value_bh, 3)
  ) %>%
  select(feature, endpoint, dataset, window, `HR (High vs Low)`, `95% CI`, p, `q (BH)`)

implied_md_header <- c(
  "| Feature | Endpoint | Dataset | Window | HR (High vs Low) | 95% CI | p | q (BH) |",
  "| --- | --- | --- | --- | --- | --- | --- | --- |"
)
implied_md_rows <- apply(implied_md, 1, function(row_values) paste0("| ", paste(row_values, collapse = " | "), " |"))
writeLines(c(implied_md_header, implied_md_rows), con = file.path(opt$output_dir, "tables", "implied_high_vs_low_by_dataset_window.md"))

signature_md <- implied_md %>% filter(feature == "Signature")
signature_md_rows <- apply(signature_md, 1, function(row_values) paste0("| ", paste(row_values, collapse = " | "), " |"))
writeLines(c(implied_md_header, signature_md_rows), con = file.path(opt$output_dir, "tables", "signature_implied_high_vs_low_by_dataset_window.md"))

write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Fit piecewise Cox models with interaction across two harmonized datasets for overall survival endpoints.",
  driver = "drivers/run_piecewise_dataset_interaction.R",
  inputs = list(
    dataset1_expr_file = opt$dataset1_expr_file,
    dataset1_meta_file = opt$dataset1_meta_file,
    dataset2_expr_file = opt$dataset2_expr_file,
    dataset2_meta_file = opt$dataset2_meta_file
  ),
  parameters = list(
    gene_column = opt$gene_column,
    gene_resolution = opt$gene_resolution,
    gene_resolution_report = file.path("tables", "gene_identifier_resolution.csv"),
    signature_genes = paste(signature_genes, collapse = ","),
    subtype = opt$subtype,
    adjust_covariates = opt$adjust_covariates,
    dataset1_stage_column = opt$dataset1_stage_column,
    dataset2_stage_column = opt$dataset2_stage_column,
    dataset1_age_column = opt$dataset1_age_column,
    dataset2_age_column = opt$dataset2_age_column,
    dataset1_label = opt$dataset1_label,
    dataset2_label = opt$dataset2_label,
    features = paste(feature_list, collapse = ","),
    endpoints = paste(endpoint_list, collapse = ","),
    cut_months = opt$cut_months
  ),
  scope = "exploratory"
)

message("Piecewise dataset-interaction summary written to: ", opt$output_dir)
