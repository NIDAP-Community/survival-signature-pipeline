# run_piecewise_subtype_specific.R
# Piecewise (0-60 vs >60 months) Cox models run within each subtype.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
})

source("src/utilities/cli_helpers.R")
source("src/utilities/run_provenance_guardrails.R")

cli_defaults <- list(
  output_dir = default_run_output_dir("piecewise_subtype_specific"),
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
  tcga_primary_filter = "tcga_primary",
  features = "Signature",
  endpoints = "OS,PFI",
  subtype_list = "BRCA_Basal,BRCA_Her2,BRCA_LumA,BRCA_LumB,BRCA_Normal",
  adjust_covariates = "stage_age",
  cut_months = "60"
)

print_usage <- function() {
  print_standard_cli_help(
    title = "Piecewise Subtype-Specific Driver",
    workflow_name = "piecewise-subtype-specific",
    script_path = file.path("drivers", "run_piecewise_subtype_specific.R"),
    defaults = cli_defaults,
    summary_lines = c(
      "This workflow expects two datasets and runs pooled piecewise models within each subtype value in subtype_list.",
      "subtype_list controls which subtypes are evaluated; both datasets are still required.",
      "adjust_covariates supports none or stage_age in the current maintained version.",
      "It emits subgroup_effect_summary.csv using the shared subgroup-summary contract for downstream heatmaps."
    )
  )
}

opt <- parse_standard_cli_args(
  defaults = cli_defaults,
  workflow_name = "piecewise-subtype-specific",
  help_printer = print_usage
)
opt$adjust_covariates <- cli_parse_enum_flag(
  value = opt$adjust_covariates,
  allowed_values = c("none", "stage_age"),
  flag_name = "adjust_covariates",
  normalize = "lower"
)
opt$cut_months <- cli_parse_numeric_flag(opt$cut_months, "cut_months")
feature_list <- cli_parse_csv_list(opt$features)
endpoint_list <- cli_parse_csv_list(opt$endpoints, uppercase = TRUE)
subtype_list <- cli_parse_csv_list(opt$subtype_list)
signature_genes <- cli_parse_csv_list(opt$signature_genes)

assert_output_guardrails(opt$output_dir)

cli_assert_existing_input_files(
  input_paths = c(
    tcga_expr_file = opt$tcga_expr_file,
    tcga_meta_file = opt$tcga_meta_file,
    metabric_expr_file = opt$metabric_expr_file,
    metabric_meta_file = opt$metabric_meta_file
  ),
  workflow_name = "piecewise-subtype-specific"
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

make_clusters <- function(x_num) {
  tertile_1 <- stats::quantile(x_num, 1 / 3, na.rm = TRUE)
  tertile_2 <- stats::quantile(x_num, 2 / 3, na.rm = TRUE)
  cluster_values <- cut(x_num, breaks = c(-Inf, tertile_1, tertile_2, Inf), labels = c("A", "B", "C"), include.lowest = TRUE, right = TRUE)
  cluster_values <- as.character(cluster_values)
  cluster_values[cluster_values == "A"] <- "Low"
  cluster_values[cluster_values == "C"] <- "High"
  cluster_values[!(cluster_values %in% c("Low", "High"))] <- NA_character_
  factor(cluster_values, levels = c("Low", "High"))
}

compute_feature_score <- function(expr, gene_column, sample_ids, feature_name, signature_genes) {
  if (tolower(feature_name) == "signature") {
    expr_sub <- expr %>% filter(.data[[gene_column]] %in% signature_genes)
    if (!nrow(expr_sub)) stop("None of signature genes found in expression matrix.")
    mat <- as.matrix(expr_sub[, setdiff(colnames(expr_sub), gene_column), drop = FALSE])
    rownames(mat) <- expr_sub[[gene_column]]
    mat <- mat[, sample_ids, drop = FALSE]
    mat <- apply(mat, 2, as.numeric)
    if (is.null(dim(mat))) mat <- matrix(mat, nrow = nrow(expr_sub))
    rownames(mat) <- expr_sub[[gene_column]]
    colnames(mat) <- sample_ids
    mat_z <- t(scale(t(mat)))
    return(colSums(mat_z, na.rm = TRUE))
  }

  expr_row <- expr %>% filter(.data[[gene_column]] == feature_name)
  if (!nrow(expr_row)) stop("Feature not found in expression matrix: ", feature_name)
  value_vector <- as.numeric(expr_row[1, setdiff(colnames(expr_row), gene_column), drop = TRUE])
  names(value_vector) <- setdiff(colnames(expr_row), gene_column)
  value_vector <- value_vector[sample_ids]
  as.numeric(scale(value_vector))
}

prepare_dataset <- function(dataset_label, expr, metadata, gene_column, sample_column, subtype_column, stage_column, age_column, endpoint, primary_filter, primary_pattern, feature_name, subtype_value, signature_genes, options) {
  if (!sample_column %in% names(metadata)) stop(dataset_label, " metadata missing sample column: ", sample_column)
  if (!subtype_column %in% names(metadata)) stop(dataset_label, " metadata missing subtype column: ", subtype_column)
  if (!stage_column %in% names(metadata)) stop(dataset_label, " metadata missing stage column: ", stage_column)
  if (!age_column %in% names(metadata)) stop(dataset_label, " metadata missing age column: ", age_column)

  endpoint_columns <- resolve_standard_endpoint_columns(endpoint, options, "piecewise-subtype-specific")
  event_column <- endpoint_columns$event
  time_column <- endpoint_columns$time

  meta_filtered <- metadata %>% filter(.data[[subtype_column]] == subtype_value)
  if (identical(primary_filter, "tcga_primary")) {
    meta_filtered <- meta_filtered %>% filter(grepl(primary_pattern, .data[[sample_column]]))
  }

  if (!all(c(event_column, time_column) %in% names(meta_filtered))) {
    stop("Missing endpoint columns in ", dataset_label, ": ", event_column, ", ", time_column)
  }

  expr_samples <- setdiff(colnames(expr), gene_column)
  meta_filtered <- meta_filtered %>% filter(.data[[sample_column]] %in% expr_samples)
  meta_filtered <- meta_filtered %>%
    mutate(
      .time_raw = parse_time_numeric(.data[[time_column]]),
      event = parse_event_numeric(.data[[event_column]])
    ) %>%
    filter(!is.na(.time_raw) & !is.na(event))

  if (!nrow(meta_filtered)) return(NULL)

  median_time <- median(meta_filtered$.time_raw, na.rm = TRUE)
  max_time <- max(meta_filtered$.time_raw, na.rm = TRUE)
  in_days <- is.finite(max_time) && (max_time > 1200 || median_time > 365)
  meta_filtered <- meta_filtered %>% mutate(time = if (in_days) .time_raw / 30.44 else .time_raw)

  score <- compute_feature_score(expr, gene_column, meta_filtered[[sample_column]], feature_name, signature_genes)

  analysis_df <- meta_filtered %>%
    transmute(
      dataset = dataset_label,
      sample_id = .data[[sample_column]],
      time = time,
      event = event,
      stage_cov = harmonize_stage(.data[[stage_column]]),
      age_cov = parse_time_numeric(.data[[age_column]]),
      feature_score = score
    ) %>%
    filter(!is.na(feature_score))

  if (!nrow(analysis_df)) return(NULL)

  analysis_df$cluster <- make_clusters(analysis_df$feature_score)
  analysis_df <- analysis_df %>% filter(!is.na(cluster))
  if (!nrow(analysis_df)) return(NULL)

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

  data.frame(
    term = term_name,
    beta = coefficient,
    HR = exp(coefficient),
    CI_low = exp(conf_int[term_name, 1]),
    CI_high = exp(conf_int[term_name, 2]),
    p_value = p_value
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

extract_lrt_p <- function(anova_result) {
  p_column <- grep("Pr", colnames(anova_result), value = TRUE)
  if (!length(p_column)) return(NA_real_)
  as.numeric(anova_result[nrow(anova_result), p_column[1]])
}

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
if (!dir.exists(file.path(opt$output_dir, "tables"))) dir.create(file.path(opt$output_dir, "tables"))
if (!dir.exists(file.path(opt$output_dir, "logs"))) dir.create(file.path(opt$output_dir, "logs"))

message("Loading expression and metadata files...")
tcga_expr <- read_expression_table(opt$tcga_expr_file, opt$gene_column)
metabric_expr <- read_expression_table(opt$metabric_expr_file, opt$gene_column)
tcga_meta <- read_metadata_table(opt$tcga_meta_file)
metabric_meta <- read_metadata_table(opt$metabric_meta_file)

summary_rows <- list()
term_rows <- list()
implied_rows <- list()

for (feature_name in feature_list) {
  for (endpoint_name in endpoint_list) {
    for (subtype_value in subtype_list) {
      message("Fitting subtype-specific model for ", feature_name, " (", endpoint_name, ", ", subtype_value, ")")

      tcga_df <- prepare_dataset(
        dataset_label = "TCGA",
        expr = tcga_expr,
        metadata = tcga_meta,
        gene_column = opt$gene_column,
        sample_column = opt$tcga_sample_column,
        subtype_column = opt$tcga_subtype_column,
        stage_column = opt$tcga_stage_column,
        age_column = opt$tcga_age_column,
        endpoint = endpoint_name,
        primary_filter = opt$tcga_primary_filter,
        primary_pattern = opt$tcga_primary_pattern,
        feature_name = feature_name,
        subtype_value = subtype_value,
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
        endpoint = endpoint_name,
        primary_filter = "none",
        primary_pattern = opt$tcga_primary_pattern,
        feature_name = feature_name,
        subtype_value = subtype_value,
        signature_genes = signature_genes,
        options = opt
      )

      if (is.null(tcga_df) || is.null(metabric_df)) next

      combined_df <- bind_rows(tcga_df, metabric_df)
      combined_df$dataset <- factor(combined_df$dataset, levels = c("TCGA", "METABRIC"))
      if (nlevels(droplevels(combined_df$dataset)) < 2) next

      early_df <- combined_df %>%
        transmute(
          dataset = dataset,
          sample_id = sample_id,
          group_high = group_high,
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
          stage_cov = stage_cov,
          age_cov = age_cov,
          tstart = opt$cut_months,
          tstop = time,
          event = event,
          post60 = 1
        )

      split_df <- bind_rows(early_df, late_df)
      if (opt$adjust_covariates == "stage_age") {
        split_df <- split_df %>%
          filter(!is.na(stage_cov), !is.na(age_cov)) %>%
          mutate(
            stage_cov = factor(stage_cov, levels = c("Early", "Late")),
            age_cov = as.numeric(age_cov)
          )
      }
      if (!nrow(split_df)) next

      formula_full <- if (opt$adjust_covariates == "stage_age") {
        survival::Surv(tstart, tstop, event) ~ group_high * dataset * post60 + stage_cov + age_cov
      } else {
        survival::Surv(tstart, tstop, event) ~ group_high * dataset * post60
      }

      formula_no_dataset_interaction <- if (opt$adjust_covariates == "stage_age") {
        survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:post60 + dataset:post60 + stage_cov + age_cov
      } else {
        survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:post60 + dataset:post60
      }

      formula_no_threeway <- if (opt$adjust_covariates == "stage_age") {
        survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:dataset + group_high:post60 + dataset:post60 + stage_cov + age_cov
      } else {
        survival::Surv(tstart, tstop, event) ~ group_high + dataset + post60 + group_high:dataset + group_high:post60 + dataset:post60
      }

      fit_full <- try(survival::coxph(formula_full, data = split_df), silent = TRUE)
      if (inherits(fit_full, "try-error")) next
      fit_no_dataset_interaction <- try(survival::coxph(formula_no_dataset_interaction, data = split_df), silent = TRUE)
      fit_no_threeway <- try(survival::coxph(formula_no_threeway, data = split_df), silent = TRUE)
      if (inherits(fit_no_dataset_interaction, "try-error") || inherits(fit_no_threeway, "try-error")) next

      lrt_dataset <- anova(fit_no_dataset_interaction, fit_full, test = "LRT")
      lrt_threeway <- anova(fit_no_threeway, fit_full, test = "LRT")

      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        subtype = subtype_value,
        feature = feature_name,
        endpoint = endpoint_name,
        covariates = opt$adjust_covariates,
        model_rows = nrow(split_df),
        model_events = sum(split_df$event == 1, na.rm = TRUE),
        p_dataset_interaction_early = extract_term(fit_full, "group_high:datasetMETABRIC")$p_value,
        p_dataset_late_shift_interaction = extract_term(fit_full, "group_high:datasetMETABRIC:post60")$p_value,
        lrt_dataset_interaction_p = extract_lrt_p(lrt_dataset),
        lrt_threeway_p = extract_lrt_p(lrt_threeway)
      )

      term_rows[[length(term_rows) + 1]] <- bind_rows(
        extract_term(fit_full, "group_high"),
        extract_term(fit_full, "group_high:datasetMETABRIC"),
        extract_term(fit_full, "group_high:post60"),
        extract_term(fit_full, "group_high:datasetMETABRIC:post60")
      ) %>% mutate(subtype = subtype_value, feature = feature_name, endpoint = endpoint_name)

      implied_rows[[length(implied_rows) + 1]] <- bind_rows(
        extract_contrast(fit_full, c("group_high" = 1)) %>% mutate(subtype = subtype_value, feature = feature_name, endpoint = endpoint_name, dataset = "TCGA", window = "0-60m"),
        extract_contrast(fit_full, c("group_high" = 1, "group_high:post60" = 1)) %>% mutate(subtype = subtype_value, feature = feature_name, endpoint = endpoint_name, dataset = "TCGA", window = ">60m"),
        extract_contrast(fit_full, c("group_high" = 1, "group_high:datasetMETABRIC" = 1)) %>% mutate(subtype = subtype_value, feature = feature_name, endpoint = endpoint_name, dataset = "METABRIC", window = "0-60m"),
        extract_contrast(fit_full, c("group_high" = 1, "group_high:datasetMETABRIC" = 1, "group_high:post60" = 1, "group_high:datasetMETABRIC:post60" = 1)) %>% mutate(subtype = subtype_value, feature = feature_name, endpoint = endpoint_name, dataset = "METABRIC", window = ">60m")
      )
    }
  }
}

if (!length(summary_rows) || !length(implied_rows)) {
  stop("No subtype-specific models converged. Check subtype sizes and covariate completeness.")
}

summary_df <- bind_rows(summary_rows) %>%
  group_by(endpoint) %>%
  mutate(
    q_dataset_interaction_early_bh = p.adjust(p_dataset_interaction_early, method = "BH"),
    q_dataset_late_shift_interaction_bh = p.adjust(p_dataset_late_shift_interaction, method = "BH"),
    q_lrt_dataset_interaction_bh = p.adjust(lrt_dataset_interaction_p, method = "BH"),
    q_lrt_threeway_bh = p.adjust(lrt_threeway_p, method = "BH")
  ) %>%
  ungroup()

terms_df <- bind_rows(term_rows) %>%
  group_by(endpoint, term) %>%
  mutate(p_value_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup()

implied_df <- bind_rows(implied_rows) %>%
  group_by(endpoint, window, subtype) %>%
  mutate(p_value_bh = p.adjust(p_value, method = "BH")) %>%
  ungroup()

subgroup_effect_df <- implied_df %>%
  transmute(
    subgroup = subtype,
    subgroup_type = "subtype",
    feature = feature,
    endpoint = endpoint,
    comparison = paste(dataset, window, sep = " | "),
    dataset = dataset,
    window = window,
    HR = HR,
    CI_low = CI_low,
    CI_high = CI_high,
    p_value = p_value,
    p_value_bh = p_value_bh,
    model = "piecewise_high_vs_low"
  )

readr::write_csv(summary_df, file.path(opt$output_dir, "tables", "subtype_piecewise_interaction_summary.csv"))
readr::write_csv(terms_df, file.path(opt$output_dir, "tables", "subtype_piecewise_interaction_terms.csv"))
readr::write_csv(implied_df, file.path(opt$output_dir, "tables", "subtype_implied_high_vs_low_by_dataset_window.csv"))
readr::write_csv(implied_df %>% filter(feature == "Signature"), file.path(opt$output_dir, "tables", "subtype_signature_implied_high_vs_low_by_dataset_window.csv"))
readr::write_csv(subgroup_effect_df, file.path(opt$output_dir, "tables", "subgroup_effect_summary.csv"))

format_num <- function(x, digits = 3) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = digits))

summary_md <- summary_df %>%
  mutate(
    `Dataset interaction p (early)` = format_num(p_dataset_interaction_early, 3),
    `Dataset interaction q (early)` = format_num(q_dataset_interaction_early_bh, 3),
    `LRT dataset interaction p` = format_num(lrt_dataset_interaction_p, 3),
    `LRT dataset interaction q` = format_num(q_lrt_dataset_interaction_bh, 3)
  ) %>%
  select(subtype, feature, endpoint, covariates, model_rows, model_events,
         `Dataset interaction p (early)`, `Dataset interaction q (early)`,
         `LRT dataset interaction p`, `LRT dataset interaction q`)

summary_header <- c(
  "| Subtype | Feature | Endpoint | Covariates | Model rows | Model events | Dataset interaction p (early) | Dataset interaction q (early) | LRT dataset interaction p | LRT dataset interaction q |",
  "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
)
summary_rows_md <- apply(summary_md, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
writeLines(c(summary_header, summary_rows_md), con = file.path(opt$output_dir, "tables", "subtype_piecewise_interaction_summary.md"))

implied_md <- implied_df %>%
  mutate(
    `HR (High vs Low)` = format_num(HR, 4),
    `95% CI` = paste0("[", format_num(CI_low, 4), ", ", format_num(CI_high, 4), "]"),
    p = format_num(p_value, 3),
    `q (BH)` = format_num(p_value_bh, 3)
  ) %>%
  select(subtype, feature, endpoint, dataset, window, `HR (High vs Low)`, `95% CI`, p, `q (BH)`)

implied_header <- c(
  "| Subtype | Feature | Endpoint | Dataset | Window | HR (High vs Low) | 95% CI | p | q (BH) |",
  "| --- | --- | --- | --- | --- | --- | --- | --- | --- |"
)
implied_rows_md <- apply(implied_md, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
writeLines(c(implied_header, implied_rows_md), con = file.path(opt$output_dir, "tables", "subtype_implied_high_vs_low_by_dataset_window.md"))

signature_md <- implied_md %>% filter(feature == "Signature")
signature_rows_md <- apply(signature_md, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
writeLines(c(implied_header, signature_rows_md), con = file.path(opt$output_dir, "tables", "subtype_signature_implied_high_vs_low_by_dataset_window.md"))

write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Run subtype-specific pooled piecewise Cox models with dataset interaction.",
  driver = "drivers/run_piecewise_subtype_specific.R",
  inputs = list(
    tcga_expr_file = opt$tcga_expr_file,
    tcga_meta_file = opt$tcga_meta_file,
    metabric_expr_file = opt$metabric_expr_file,
    metabric_meta_file = opt$metabric_meta_file
  ),
  parameters = list(
    gene_column = opt$gene_column,
    signature_genes = paste(signature_genes, collapse = ","),
    features = paste(feature_list, collapse = ","),
    endpoints = paste(endpoint_list, collapse = ","),
    subtype_list = paste(subtype_list, collapse = ","),
    adjust_covariates = opt$adjust_covariates,
    cut_months = opt$cut_months,
    tcga_primary_pattern = opt$tcga_primary_pattern
  ),
  scope = "exploratory"
)

message("Subtype-specific piecewise outputs written to: ", opt$output_dir)
