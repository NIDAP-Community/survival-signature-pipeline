#!/usr/bin/env Rscript

# Execute a declarative set of continuous-feature Cox models for one cohort.
# The plan is deliberately cohort-agnostic: it describes model intent and
# metadata columns, but it does not contain disease- or project-specific logic.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(survival)
})

source("src/utilities/cli_helpers.R")
source("src/utilities/gene_identifier_resolution.R")
source("src/utilities/run_provenance_guardrails.R")

cli_defaults <- list(
  output_dir = default_run_output_dir("single_cohort_analysis_plan"),
  expr_file = file.path("test_data", "tcga_expression.tsv"),
  meta_file = file.path("test_data", "tcga_metadata.csv"),
  plan_file = file.path("test_data", "single_cohort_analysis_plan.csv"),
  sample_column = "sampnames",
  gene_column = "Gene",
  gene_resolution = "none",
  os_event_column = "OS",
  os_time_column = "OS_time",
  pfi_event_column = "PFI",
  pfi_time_column = "PFI_time",
  verbose = "true"
)

print_usage <- function() {
  print_standard_cli_help(
    title = "Single-Cohort Analysis-Plan Driver",
    workflow_name = "single-cohort-analysis-plan",
    script_path = file.path("drivers", "run_single_cohort_analysis_plan.R"),
    defaults = cli_defaults,
    summary_lines = c(
      "Runs generic overall, filtered-subgroup, and feature-by-group interaction Cox models described in a CSV plan.",
      "The plan, not the driver, defines cohort-specific metadata columns, filters, covariates, and interaction axes.",
      "Each successful model emits a shared subgroup-summary row; interaction rows are marked subgroup_type=interaction_model.",
      "gene_resolution=l2p opt-in updates user-supplied human gene symbols and writes an auditable mapping table without changing expression identifiers.",
      "For endpoints other than OS or PFI, provide event_column and time_column in the plan."
    )
  )
}

is_truthy <- function(value) {
  tolower(as.character(value)) %in% c("true", "t", "1", "yes", "y")
}

value_or_default <- function(value, default) {
  value <- trimws(as.character(value))
  value[is.na(value) | !nzchar(value)] <- default
  value
}

quote_formula_name <- function(value) {
  paste0("`", gsub("`", "", value, fixed = TRUE), "`")
}

parse_plan_list <- function(value) {
  values <- cli_parse_csv_list(value_or_default(value, ""), allow_empty = TRUE)
  values[tolower(values) != "none"]
}

read_analysis_plan <- function(path) {
  plan <- readr::read_csv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    name_repair = "minimal",
    show_col_types = FALSE
  )
  names(plan) <- trimws(names(plan))

  required <- c("analysis_id", "analysis_type", "feature", "endpoint")
  missing <- setdiff(required, names(plan))
  if (length(missing)) {
    stop("Analysis plan is missing required columns: ", paste(missing, collapse = ", "))
  }
  if (!nrow(plan)) stop("Analysis plan contains no model rows.")

  optional_defaults <- c(
    feature_type = "gene",
    signature_genes = "",
    covariates = "",
    filter_column = "",
    filter_values = "",
    subgroup_type = "subgroup",
    subgroup_label = "",
    interaction_column = "",
    interaction_reference = "",
    interaction_label = "",
    event_column = "",
    time_column = ""
  )
  for (column_name in names(optional_defaults)) {
    if (!column_name %in% names(plan)) plan[[column_name]] <- optional_defaults[[column_name]]
    plan[[column_name]] <- value_or_default(plan[[column_name]], optional_defaults[[column_name]])
  }

  plan <- plan %>% mutate(
    analysis_id = trimws(.data$analysis_id),
    analysis_type = tolower(trimws(.data$analysis_type)),
    feature = trimws(.data$feature),
    feature_type = tolower(trimws(.data$feature_type)),
    endpoint = trimws(.data$endpoint)
  )

  if (any(!nzchar(plan$analysis_id)) || anyDuplicated(plan$analysis_id)) {
    stop("analysis_id must be non-empty and unique in the analysis plan.")
  }
  if (any(!nzchar(plan$feature))) stop("feature must be non-empty for every analysis-plan row.")
  if (any(!nzchar(plan$endpoint))) stop("endpoint must be non-empty for every analysis-plan row.")

  valid_types <- c("overall", "subgroup", "interaction")
  invalid_types <- setdiff(unique(plan$analysis_type), valid_types)
  if (length(invalid_types)) {
    stop("analysis_type must be one of: ", paste(valid_types, collapse = ", "), ". Invalid values: ", paste(invalid_types, collapse = ", "))
  }
  invalid_feature_types <- setdiff(unique(plan$feature_type), c("gene", "signature"))
  if (length(invalid_feature_types)) {
    stop("feature_type must be gene or signature. Invalid values: ", paste(invalid_feature_types, collapse = ", "))
  }
  if (any(plan$feature_type == "signature" & !nzchar(plan$signature_genes))) {
    stop("signature_genes is required for feature_type=signature.")
  }
  if (any(xor(nzchar(plan$filter_column), nzchar(plan$filter_values)))) {
    stop("filter_column and filter_values must either both be supplied or both be blank.")
  }
  if (any(plan$analysis_type == "interaction" & !nzchar(plan$interaction_column))) {
    stop("interaction_column is required for analysis_type=interaction.")
  }
  if (any(xor(nzchar(plan$event_column), nzchar(plan$time_column)))) {
    stop("event_column and time_column must either both be supplied or both be blank.")
  }

  plan
}

resolve_endpoint_columns <- function(plan_row, options) {
  if (nzchar(plan_row$event_column) && nzchar(plan_row$time_column)) {
    return(list(event = plan_row$event_column, time = plan_row$time_column))
  }
  resolve_standard_endpoint_columns(plan_row$endpoint, options, "single-cohort-analysis-plan")
}

build_feature_values <- function(plan_row, expression_data, gene_column) {
  if (plan_row$feature_type == "gene") {
    hits <- which(expression_data[[gene_column]] == plan_row$feature)
    if (length(hits) != 1L) {
      stop("Feature '", plan_row$feature, "' must match exactly one expression row; found ", length(hits), ".")
    }
    values <- suppressWarnings(as.numeric(expression_data[hits, setdiff(names(expression_data), gene_column), drop = TRUE]))
    names(values) <- setdiff(names(expression_data), gene_column)
    return(values)
  }

  signature_genes <- parse_plan_list(plan_row$signature_genes)
  signature_rows <- match(signature_genes, expression_data[[gene_column]])
  missing_genes <- signature_genes[is.na(signature_rows)]
  if (length(missing_genes)) {
    stop("Signature '", plan_row$feature, "' is missing genes: ", paste(missing_genes, collapse = ", "))
  }

  sample_columns <- setdiff(names(expression_data), gene_column)
  signature_matrix <- as.matrix(expression_data[signature_rows, sample_columns, drop = FALSE])
  storage.mode(signature_matrix) <- "numeric"
  values <- colMeans(signature_matrix, na.rm = FALSE)
  names(values) <- sample_columns
  values
}

resolve_plan_gene_identifiers <- function(plan, expression_gene_ids, mode) {
  gene_inputs <- c(
    plan$feature[plan$feature_type == "gene"],
    unlist(lapply(plan$signature_genes[plan$feature_type == "signature"], parse_plan_list), use.names = FALSE)
  )
  gene_inputs <- unique(gene_inputs)
  resolution_table <- resolve_requested_gene_identifiers(gene_inputs, expression_gene_ids, mode)
  assert_gene_resolution_is_unambiguous(resolution_table)

  plan$feature_requested <- plan$feature
  plan$signature_genes_requested <- plan$signature_genes
  gene_rows <- plan$feature_type == "gene"
  plan$feature[gene_rows] <- replace_requested_gene_identifiers(plan$feature[gene_rows], resolution_table)
  signature_rows <- plan$feature_type == "signature"
  plan$signature_genes[signature_rows] <- vapply(plan$signature_genes[signature_rows], function(value) {
    paste(replace_requested_gene_identifiers(parse_plan_list(value), resolution_table), collapse = ",")
  }, character(1))

  list(plan = plan, resolution_table = resolution_table)
}

prepare_model_data <- function(plan_row, expression_data, metadata, options) {
  endpoint_columns <- resolve_endpoint_columns(plan_row, options)
  required_metadata <- unique(c(
    options$sample_column,
    endpoint_columns$event,
    endpoint_columns$time,
    parse_plan_list(plan_row$covariates),
    if (nzchar(plan_row$filter_column)) plan_row$filter_column else character(0),
    if (plan_row$analysis_type == "interaction") plan_row$interaction_column else character(0)
  ))
  missing_metadata <- setdiff(required_metadata, names(metadata))
  if (length(missing_metadata)) {
    stop("Metadata is missing required column(s) for analysis_id '", plan_row$analysis_id, "': ", paste(missing_metadata, collapse = ", "))
  }

  model_data <- metadata
  if (nzchar(plan_row$filter_column)) {
    requested_values <- parse_plan_list(plan_row$filter_values)
    model_data <- model_data %>% filter(.data[[plan_row$filter_column]] %in% requested_values)
  }
  if (!nrow(model_data)) stop("No metadata rows remain after the configured filter for analysis_id '", plan_row$analysis_id, "'.")

  feature_values <- build_feature_values(plan_row, expression_data, options$gene_column)
  model_data$feature_value <- unname(feature_values[as.character(model_data[[options$sample_column]])])
  model_data$time_value <- suppressWarnings(as.numeric(model_data[[endpoint_columns$time]]))
  model_data$event_value <- suppressWarnings(as.numeric(model_data[[endpoint_columns$event]]))

  covariates <- parse_plan_list(plan_row$covariates)
  if (any(covariates %in% c("feature_value", "time_value", "event_value", "interaction_factor"))) {
    stop("Plan covariates may not use reserved model-data column names for analysis_id '", plan_row$analysis_id, "'.")
  }

  if (plan_row$analysis_type == "interaction") {
    if (plan_row$interaction_column %in% covariates) {
      stop("interaction_column must not also be listed in covariates for analysis_id '", plan_row$analysis_id, "'.")
    }
    model_data$interaction_factor <- as.factor(model_data[[plan_row$interaction_column]])
    if (nzchar(plan_row$interaction_reference)) {
      if (!plan_row$interaction_reference %in% levels(model_data$interaction_factor)) {
        stop("interaction_reference is not an observed level for analysis_id '", plan_row$analysis_id, "'.")
      }
      model_data$interaction_factor <- stats::relevel(model_data$interaction_factor, ref = plan_row$interaction_reference)
    }
  }

  complete_columns <- c("feature_value", "time_value", "event_value", covariates)
  if (plan_row$analysis_type == "interaction") complete_columns <- c(complete_columns, "interaction_factor")
  model_data <- model_data[stats::complete.cases(model_data[, complete_columns, drop = FALSE]), , drop = FALSE]
  model_data <- model_data[is.finite(model_data$feature_value) & is.finite(model_data$time_value) & is.finite(model_data$event_value), , drop = FALSE]
  if (!nrow(model_data)) stop("No complete, finite observations remain for analysis_id '", plan_row$analysis_id, "'.")
  if (!all(model_data$event_value %in% c(0, 1))) stop("Event values must be coded 0 or 1 for analysis_id '", plan_row$analysis_id, "'.")
  if (any(model_data$time_value < 0)) stop("Survival times must be non-negative for analysis_id '", plan_row$analysis_id, "'.")
  if (plan_row$analysis_type == "interaction" && nlevels(model_data$interaction_factor) < 2L) {
    stop("Interaction model requires at least two observed interaction levels for analysis_id '", plan_row$analysis_id, "'.")
  }

  list(data = model_data, covariates = covariates, endpoint_columns = endpoint_columns)
}

fit_plan_model <- function(plan_row, prepared) {
  model_terms <- if (plan_row$analysis_type == "interaction") {
    c("feature_value * interaction_factor", vapply(prepared$covariates, quote_formula_name, character(1)))
  } else {
    c("feature_value", vapply(prepared$covariates, quote_formula_name, character(1)))
  }
  formula_text <- paste("survival::Surv(time_value, event_value) ~", paste(model_terms, collapse = " + "))
  fit <- survival::coxph(stats::as.formula(formula_text), data = prepared$data, x = TRUE)
  list(fit = fit, formula_text = formula_text)
}

coefficient_table <- function(fit, plan_row, formula_text) {
  model_summary <- summary(fit)
  if (is.null(model_summary$coefficients) || !nrow(model_summary$coefficients)) return(data.frame())
  data.frame(
    analysis_id = plan_row$analysis_id,
    analysis_type = plan_row$analysis_type,
    feature = plan_row$feature,
    endpoint = plan_row$endpoint,
    term = rownames(model_summary$coefficients),
    coef = unname(model_summary$coefficients[, "coef"]),
    HR = unname(model_summary$conf.int[, "exp(coef)"]),
    CI_low = unname(model_summary$conf.int[, "lower .95"]),
    CI_high = unname(model_summary$conf.int[, "upper .95"]),
    p_value = unname(model_summary$coefficients[, "Pr(>|z|)"]),
    N = model_summary$n,
    events = model_summary$nevent,
    model_formula = formula_text,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

make_effect_rows <- function(plan_row, coefficient_rows, prepared) {
  if (!nrow(coefficient_rows)) return(data.frame())
  if (plan_row$analysis_type == "interaction") {
    interaction_rows <- coefficient_rows[grepl("feature_value:interaction_factor", coefficient_rows$term, fixed = TRUE), , drop = FALSE]
    if (!nrow(interaction_rows)) return(data.frame())
    levels <- levels(prepared$data$interaction_factor)
    reference <- levels[[1]]
    suffixes <- sub("^feature_value:interaction_factor", "", interaction_rows$term)
    level_index <- match(suffixes, make.names(levels))
    compared_level <- ifelse(is.na(level_index), suffixes, levels[level_index])
    interaction_label <- value_or_default(plan_row$interaction_label, plan_row$interaction_column)
    return(interaction_rows %>% transmute(
      subgroup = paste0(plan_row$analysis_id, "__", sanitize_path_component(compared_level)),
      subgroup_label = paste0(interaction_label, ": ", compared_level, " vs ", reference),
      subgroup_type = "interaction_model",
      feature = .data$feature,
      endpoint = .data$endpoint,
      comparison = paste0("Interaction: feature x ", interaction_label, " (", compared_level, " vs ", reference, ")"),
      HR = .data$HR,
      CI_low = .data$CI_low,
      CI_high = .data$CI_high,
      p_value = .data$p_value,
      N = .data$N,
      events = .data$events,
      analysis_id = .data$analysis_id,
      analysis_type = .data$analysis_type,
      model_formula = .data$model_formula,
      covariates = plan_row$covariates,
      filter_column = plan_row$filter_column,
      filter_values = plan_row$filter_values,
      interaction_column = plan_row$interaction_column,
      interaction_reference = reference
    ))
  }

  feature_rows <- coefficient_rows[coefficient_rows$term == "feature_value", , drop = FALSE]
  if (!nrow(feature_rows)) return(data.frame())
  subgroup_type <- if (plan_row$analysis_type == "overall") "overall" else value_or_default(plan_row$subgroup_type, "subgroup")
  subgroup_label <- value_or_default(plan_row$subgroup_label, plan_row$analysis_id)
  feature_rows %>% transmute(
    subgroup = plan_row$analysis_id,
    subgroup_label = subgroup_label,
    subgroup_type = subgroup_type,
    feature = .data$feature,
    endpoint = .data$endpoint,
    comparison = "Adjusted continuous Cox",
    HR = .data$HR,
    CI_low = .data$CI_low,
    CI_high = .data$CI_high,
    p_value = .data$p_value,
    N = .data$N,
    events = .data$events,
    analysis_id = .data$analysis_id,
    analysis_type = .data$analysis_type,
    model_formula = .data$model_formula,
    covariates = plan_row$covariates,
    filter_column = plan_row$filter_column,
    filter_values = plan_row$filter_values,
    interaction_column = "",
    interaction_reference = ""
  )
}

opt <- parse_standard_cli_args(
  defaults = cli_defaults,
  workflow_name = "single-cohort-analysis-plan",
  help_printer = print_usage
)
opt$gene_resolution <- normalize_gene_resolution_mode(opt$gene_resolution)
assert_output_guardrails(opt$output_dir)
cli_assert_existing_input_files(
  input_paths = c(expr_file = opt$expr_file, meta_file = opt$meta_file, plan_file = opt$plan_file),
  workflow_name = "single-cohort-analysis-plan"
)

expression_data <- read_expression_table(opt$expr_file, opt$gene_column)
metadata <- read_metadata_table(opt$meta_file)
if (!opt$sample_column %in% names(metadata)) stop("Metadata missing sample column: ", opt$sample_column)
if (anyDuplicated(metadata[[opt$sample_column]])) stop("Metadata sample column contains duplicate sample identifiers.")
if (!opt$gene_column %in% names(expression_data)) stop("Expression data missing gene column: ", opt$gene_column)
if (anyDuplicated(expression_data[[opt$gene_column]])) stop("Expression gene column contains duplicate feature identifiers.")

plan <- read_analysis_plan(opt$plan_file)
resolved_plan <- resolve_plan_gene_identifiers(
  plan = plan,
  expression_gene_ids = expression_data[[opt$gene_column]],
  mode = opt$gene_resolution
)
plan <- resolved_plan$plan
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
readr::write_csv(plan, file.path(opt$output_dir, "tables", "analysis_plan_resolved.csv"))
write_gene_resolution_report(resolved_plan$resolution_table, opt$output_dir)

effect_rows <- list()
coefficient_rows <- list()
status_rows <- list()

for (row_index in seq_len(nrow(plan))) {
  plan_row <- as.list(plan[row_index, , drop = FALSE])
  names(plan_row) <- names(plan)
  plan_row <- lapply(plan_row, function(value) as.character(value[[1]]))
  analysis_id <- plan_row$analysis_id

  attempt <- tryCatch({
    prepared <- prepare_model_data(plan_row, expression_data, metadata, opt)
    fitted <- fit_plan_model(plan_row, prepared)
    coefficients <- coefficient_table(fitted$fit, plan_row, fitted$formula_text)
    effects <- make_effect_rows(plan_row, coefficients, prepared)
    list(status = "completed", message = "", coefficients = coefficients, effects = effects)
  }, error = function(error) {
    list(status = "failed", message = conditionMessage(error), coefficients = data.frame(), effects = data.frame())
  })

  status_rows[[length(status_rows) + 1L]] <- data.frame(
    analysis_id = analysis_id,
    analysis_type = plan_row$analysis_type,
    status = attempt$status,
    message = attempt$message,
    stringsAsFactors = FALSE
  )
  if (nrow(attempt$coefficients)) coefficient_rows[[length(coefficient_rows) + 1L]] <- attempt$coefficients
  if (nrow(attempt$effects)) effect_rows[[length(effect_rows) + 1L]] <- attempt$effects
  if (is_truthy(opt$verbose)) message("[", attempt$status, "] ", analysis_id, if (nzchar(attempt$message)) paste0(": ", attempt$message) else "")
}

status_df <- bind_rows(status_rows)
coefficients_df <- bind_rows(coefficient_rows)
effects_df <- bind_rows(effect_rows)
if (nrow(effects_df)) {
  effects_df <- effects_df %>% group_by(.data$feature, .data$endpoint, .data$analysis_type) %>%
    mutate(p_value_bh = stats::p.adjust(.data$p_value, method = "BH")) %>%
    ungroup() %>%
    select(subgroup, subgroup_label, subgroup_type, feature, endpoint, comparison, HR, CI_low, CI_high, p_value, p_value_bh, everything())
}

readr::write_csv(status_df, file.path(opt$output_dir, "tables", "analysis_plan_status.csv"))
readr::write_csv(coefficients_df, file.path(opt$output_dir, "tables", "model_coefficients.csv"))
readr::write_csv(effects_df, file.path(opt$output_dir, "tables", "subgroup_effect_summary.csv"))

write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Run a declarative single-cohort survival analysis plan and emit shared subgroup-summary results.",
  driver = "drivers/run_single_cohort_analysis_plan.R",
  inputs = list(expr_file = opt$expr_file, meta_file = opt$meta_file, plan_file = opt$plan_file),
  parameters = list(
    sample_column = opt$sample_column,
    gene_column = opt$gene_column,
    gene_resolution = opt$gene_resolution,
    gene_resolution_report = file.path("tables", "gene_identifier_resolution.csv"),
    plan_rows = nrow(plan)
  ),
  notes = c(
    "Plan rows are independent model specifications; failed rows are recorded in tables/analysis_plan_status.csv.",
    "p_value_bh is adjusted within feature, endpoint, and analysis_type among emitted primary effect rows."
  ),
  scope = "exploratory"
)

if (any(status_df$status == "failed")) {
  message("Analysis plan completed with one or more failed rows; inspect tables/analysis_plan_status.csv.")
} else {
  message("Analysis plan completed successfully.")
}
message("Single-cohort analysis-plan outputs written to: ", opt$output_dir)
