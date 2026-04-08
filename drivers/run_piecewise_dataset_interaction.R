# run_piecewise_dataset_interaction.R
# Piecewise (0-60 vs >60 months) Cox models with dataset interaction (TCGA vs METABRIC).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
})

source("src/utilities/run_provenance_guardrails.R")

cli_defaults <- list(
  output_dir = default_run_output_dir("piecewise_dataset_interaction"),
  tcga_expr_file = file.path("test_data", "tcga_expression.tsv"),
  tcga_meta_file = file.path("test_data", "tcga_metadata.csv"),
  metabric_expr_file = file.path("test_data", "metabric_expression.tsv"),
  metabric_meta_file = file.path("test_data", "metabric_metadata.csv"),
  tcga_sample_column = "sampnames",
  metabric_sample_column = "sampnames",
  tcga_subtype_column = "Subtype",
  metabric_subtype_column = "Subtype",
  tcga_stage_column = "ajcc_pathologic_tumor_stage",
  metabric_stage_column = "stage_collapsed",
  tcga_age_column = "Diagnosis_Age",
  metabric_age_column = "Diagnosis_Age",
  gene_column = "Gene",
  signature_genes = "SAT1,SRM,ODC1,AMD1",
  os_event_column = "OS",
  os_time_column = "OS_time",
  pfi_event_column = "PFI",
  pfi_time_column = "PFI_time",
  tcga_primary_pattern = "_01",
  subtype = "ALL",
  tcga_primary_filter = "tcga_primary",
  adjust_covariates = "none",
  features = "AMD1,SRM,SAT1,ODC1,Signature",
  endpoints = "OS,PFI",
  cut_months = "60"
)

print_usage <- function(defaults) {
  cat(
    paste(
      "Piecewise Dataset Interaction Driver",
      "",
      "Usage:",
      "  Rscript drivers/run_piecewise_dataset_interaction.R [--flag value]...",
      "  Rscript entrypoint.R piecewise-dataset-interaction [--flag value]...",
      "",
      "Key notes:",
      "  - This workflow expects two datasets and models dataset interaction terms across TCGA and METABRIC inputs.",
      "  - Endpoint columns are configurable through the OS/PFI event and time flags.",
      "  - tcga_primary_filter uses tcga_primary_pattern to subset TCGA sample IDs when enabled.",
      "",
      "Important flags:",
      paste(vapply(names(defaults), function(name) {
        paste0("  --", gsub("_", "-", name), " (default: ", defaults[[name]], ")")
      }, character(1)), collapse = "\n"),
      sep = "\n"
    ),
    "\n"
  )
}

parse_cli <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || any(args %in% c("--help", "-h", "help"))) {
    print_usage(defaults)
    quit(save = "no", status = 0)
  }
  if (!length(args)) return(defaults)
  index <- 1
  while (index <= length(args)) {
    key <- args[[index]]
    if (!startsWith(key, "--")) stop("Arguments must be supplied as --flag value pairs. Offending token: ", key)
    key_clean <- gsub("-", "_", substring(key, 3))
    if (!key_clean %in% names(defaults)) stop("Unknown argument: ", key)
    if (index == length(args)) stop("No value supplied for argument ", key)
    defaults[[key_clean]] <- args[[index + 1]]
    index <- index + 2
  }
  defaults
}

opt <- parse_cli(cli_defaults)
opt$cut_months <- as.numeric(opt$cut_months)
opt$adjust_covariates <- tolower(opt$adjust_covariates)
if (!opt$adjust_covariates %in% c("none", "subtype_stage", "subtype_stage_age")) {
  stop("--adjust-covariates must be one of: none, subtype_stage, subtype_stage_age")
}
feature_list <- trimws(unlist(strsplit(opt$features, ",")))
endpoint_list <- toupper(trimws(unlist(strsplit(opt$endpoints, ","))))
signature_genes <- trimws(unlist(strsplit(opt$signature_genes, ",")))

assert_output_guardrails(opt$output_dir)

resolve_endpoint_columns <- function(endpoint_name, options) {
  endpoint_name <- toupper(endpoint_name)
  if (endpoint_name == "OS") {
    return(list(event = options$os_event_column, time = options$os_time_column))
  }
  if (endpoint_name == "PFI") {
    return(list(event = options$pfi_event_column, time = options$pfi_time_column))
  }
  stop("Unsupported endpoint: ", endpoint_name)
}

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

read_expression <- function(path, gene_column) {
  expression_data <- readr::read_tsv(path, show_col_types = FALSE)
  if (!gene_column %in% colnames(expression_data)) colnames(expression_data)[1] <- gene_column
  expression_data
}

read_metadata <- function(path) {
  metadata <- readr::read_csv(path, show_col_types = FALSE)
  colnames(metadata) <- gsub("\\n", "", colnames(metadata))
  metadata
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

prepare_dataset <- function(dataset_label, expr, metadata, gene_column, sample_column, subtype_column, stage_column, age_column, subtype_filter, endpoint, primary_filter, primary_pattern, feature_name, signature_genes, options) {
  if (!sample_column %in% names(metadata)) stop(dataset_label, " metadata missing sample column: ", sample_column)
  if (!subtype_column %in% names(metadata)) stop(dataset_label, " metadata missing subtype column: ", subtype_column)
  if (!stage_column %in% names(metadata)) stop(dataset_label, " metadata missing stage column: ", stage_column)
  if (!age_column %in% names(metadata)) stop(dataset_label, " metadata missing age column: ", age_column)

  endpoint_columns <- resolve_endpoint_columns(endpoint, options)
  event_column <- endpoint_columns$event
  time_column <- endpoint_columns$time

  metadata_filtered <- metadata
  subtype_clean <- tolower(trimws(subtype_filter))
  if (!(subtype_clean %in% c("", "all", "overall", "any"))) {
    selected_subtypes <- trimws(unlist(strsplit(subtype_filter, ",")))
    metadata_filtered <- metadata_filtered %>% filter(.data[[subtype_column]] %in% selected_subtypes)
  }

  if (identical(primary_filter, "tcga_primary")) {
    metadata_filtered <- metadata_filtered %>% filter(grepl(primary_pattern, .data[[sample_column]]))
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
tcga_expr <- read_expression(opt$tcga_expr_file, opt$gene_column)
metabric_expr <- read_expression(opt$metabric_expr_file, opt$gene_column)
tcga_meta <- read_metadata(opt$tcga_meta_file)
metabric_meta <- read_metadata(opt$metabric_meta_file)

result_rows <- list()
term_rows <- list()
implied_rows <- list()
all_term_rows <- list()

for (feature_name in feature_list) {
  for (endpoint_name in endpoint_list) {
    message("Fitting piecewise model for ", feature_name, " (", endpoint_name, ")")

    tcga_df <- prepare_dataset(
      dataset_label = "TCGA",
      expr = tcga_expr,
      metadata = tcga_meta,
      gene_column = opt$gene_column,
      sample_column = opt$tcga_sample_column,
      subtype_column = opt$tcga_subtype_column,
      stage_column = opt$tcga_stage_column,
      age_column = opt$tcga_age_column,
      subtype_filter = opt$subtype,
      endpoint = endpoint_name,
      primary_filter = opt$tcga_primary_filter,
      primary_pattern = opt$tcga_primary_pattern,
      feature_name = feature_name,
      signature_genes = signature_genes,
      options = opt
    )

    metabric_df <- prepare_dataset(
      dataset_label = "METABRIC",
      expr = metabric_expr,
      metadata = metabric_meta,
      gene_column = opt$gene_column,
      sample_column = opt$metabric_sample_column,
      subtype_column = opt$metabric_subtype_column,
      stage_column = opt$metabric_stage_column,
      age_column = opt$metabric_age_column,
      subtype_filter = opt$subtype,
      endpoint = endpoint_name,
      primary_filter = "none",
      primary_pattern = opt$tcga_primary_pattern,
      feature_name = feature_name,
      signature_genes = signature_genes,
      options = opt
    )

    combined_df <- bind_rows(tcga_df, metabric_df)
    combined_df$dataset <- factor(combined_df$dataset, levels = c("TCGA", "METABRIC"))

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

    term_table <- bind_rows(
      extract_term(fit_full, "group_high"),
      extract_term(fit_full, "group_high:datasetMETABRIC"),
      extract_term(fit_full, "group_high:post60"),
      extract_term(fit_full, "group_high:datasetMETABRIC:post60")
    ) %>%
      mutate(feature = feature_name, endpoint = endpoint_name)

    term_rows[[length(term_rows) + 1]] <- term_table
    all_term_rows[[length(all_term_rows) + 1]] <- extract_all_terms(fit_full) %>%
      mutate(feature = feature_name, endpoint = endpoint_name, covariates = covariate_mode)

    tcga_counts <- tcga_df %>% group_by(cluster) %>% summarise(n = n(), events = sum(event == 1), .groups = "drop")
    metabric_counts <- metabric_df %>% group_by(cluster) %>% summarise(n = n(), events = sum(event == 1), .groups = "drop")

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
      tcga_low_n = pull_count(tcga_counts, "Low", "n"),
      tcga_high_n = pull_count(tcga_counts, "High", "n"),
      tcga_low_events = pull_count(tcga_counts, "Low", "events"),
      tcga_high_events = pull_count(tcga_counts, "High", "events"),
      metabric_low_n = pull_count(metabric_counts, "Low", "n"),
      metabric_high_n = pull_count(metabric_counts, "High", "n"),
      metabric_low_events = pull_count(metabric_counts, "Low", "events"),
      metabric_high_events = pull_count(metabric_counts, "High", "events"),
      hr_high_tcga_early = extract_term(fit_full, "group_high")$HR,
      p_high_tcga_early = extract_term(fit_full, "group_high")$p_value,
      hr_dataset_interaction_early = extract_term(fit_full, "group_high:datasetMETABRIC")$HR,
      p_dataset_interaction_early = extract_term(fit_full, "group_high:datasetMETABRIC")$p_value,
      hr_tcga_late_shift = extract_term(fit_full, "group_high:post60")$HR,
      p_tcga_late_shift = extract_term(fit_full, "group_high:post60")$p_value,
      hr_dataset_late_shift_interaction = extract_term(fit_full, "group_high:datasetMETABRIC:post60")$HR,
      p_dataset_late_shift_interaction = extract_term(fit_full, "group_high:datasetMETABRIC:post60")$p_value,
      lrt_dataset_interaction_p = extract_lrt_p(lrt_dataset),
      lrt_threeway_p = extract_lrt_p(lrt_threeway)
    )

    tcga_early <- extract_contrast(fit_full, c("group_high" = 1))
    tcga_late <- extract_contrast(fit_full, c("group_high" = 1, "group_high:post60" = 1))
    metabric_early <- extract_contrast(fit_full, c("group_high" = 1, "group_high:datasetMETABRIC" = 1))
    metabric_late <- extract_contrast(
      fit_full,
      c(
        "group_high" = 1,
        "group_high:datasetMETABRIC" = 1,
        "group_high:post60" = 1,
        "group_high:datasetMETABRIC:post60" = 1
      )
    )

    implied_rows[[length(implied_rows) + 1]] <- bind_rows(
      tcga_early %>% mutate(feature = feature_name, endpoint = endpoint_name, dataset = "TCGA", window = "0-60m"),
      tcga_late %>% mutate(feature = feature_name, endpoint = endpoint_name, dataset = "TCGA", window = ">60m"),
      metabric_early %>% mutate(feature = feature_name, endpoint = endpoint_name, dataset = "METABRIC", window = "0-60m"),
      metabric_late %>% mutate(feature = feature_name, endpoint = endpoint_name, dataset = "METABRIC", window = ">60m")
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
  "| Feature | Endpoint | Covariates | Model rows | Model events | TCGA early HR (High vs Low) | TCGA early p | Dataset interaction HR (early) | Dataset interaction p (early) | Dataset interaction q (BH) | TCGA late-shift HR | TCGA late-shift p | Dataset late-shift interaction HR | Dataset late-shift interaction p | Dataset late-shift interaction q (BH) | LRT dataset interaction p | LRT dataset interaction q (BH) | LRT three-way p | LRT three-way q (BH) |",
  "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
)

md_rows <- apply(results_df, 1, function(row_values) {
  paste0(
    "| ", row_values[["feature"]],
    " | ", row_values[["endpoint"]],
    " | ", row_values[["model_covariates"]],
    " | ", format_num(as.numeric(row_values[["model_rows"]]), 4),
    " | ", format_num(as.numeric(row_values[["model_events"]]), 4),
    " | ", format_num(as.numeric(row_values[["hr_high_tcga_early"]]), 4),
    " | ", format_num(as.numeric(row_values[["p_high_tcga_early"]]), 3),
    " | ", format_num(as.numeric(row_values[["hr_dataset_interaction_early"]]), 4),
    " | ", format_num(as.numeric(row_values[["p_dataset_interaction_early"]]), 3),
    " | ", format_num(as.numeric(row_values[["q_dataset_interaction_early_bh"]]), 3),
    " | ", format_num(as.numeric(row_values[["hr_tcga_late_shift"]]), 4),
    " | ", format_num(as.numeric(row_values[["p_tcga_late_shift"]]), 3),
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
  purpose = "Fit piecewise Cox models with TCGA vs METABRIC interaction for overall survival endpoints.",
  driver = "drivers/run_piecewise_dataset_interaction.R",
  inputs = list(
    tcga_expr_file = opt$tcga_expr_file,
    tcga_meta_file = opt$tcga_meta_file,
    metabric_expr_file = opt$metabric_expr_file,
    metabric_meta_file = opt$metabric_meta_file
  ),
  parameters = list(
    gene_column = opt$gene_column,
    signature_genes = paste(signature_genes, collapse = ","),
    subtype = opt$subtype,
    adjust_covariates = opt$adjust_covariates,
    tcga_stage_column = opt$tcga_stage_column,
    metabric_stage_column = opt$metabric_stage_column,
    tcga_age_column = opt$tcga_age_column,
    metabric_age_column = opt$metabric_age_column,
    features = paste(feature_list, collapse = ","),
    endpoints = paste(endpoint_list, collapse = ","),
    cut_months = opt$cut_months,
    tcga_primary_pattern = opt$tcga_primary_pattern
  ),
  scope = "exploratory"
)

message("Piecewise dataset-interaction summary written to: ", opt$output_dir)
