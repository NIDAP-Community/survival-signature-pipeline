# run_signature_consistency_summary.R
# Build a compact comparison table: overall adjusted Signature vs subtype-specific adjusted Signature.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("src/utilities/run_provenance_guardrails.R")

cli_defaults <- list(
  output_dir = default_run_output_dir("signature_consistency_summary"),
  overall_implied_csv = "",
  overall_summary_csv = "",
  subtype_implied_csv = "",
  subtype_summary_csv = "",
  subtype_terms_csv = "",
  signature_feature_name = "Signature",
  subtype_term_pattern = "^subtype_cov"
)

print_usage <- function(defaults) {
  cat(
    paste(
      "Signature Consistency Summary Driver",
      "",
      "Usage:",
      "  Rscript drivers/run_signature_consistency_summary.R [--flag value]...",
      "  Rscript entrypoint.R signature-consistency-summary [--flag value]...",
      "",
      "Key notes:",
      "  - This workflow is downstream-only and requires explicit upstream result tables.",
      "  - No raw expression or metadata inputs are used here.",
      "  - signature_feature_name and subtype_term_pattern control which summary rows are extracted.",
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
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) stop("Arguments must be supplied as --flag value pairs. Offending token: ", key)
    key_clean <- gsub("-", "_", substring(key, 3))
    if (!key_clean %in% names(defaults)) stop("Unknown argument: ", key)
    if (i == length(args)) stop("No value supplied for argument ", key)
    defaults[[key_clean]] <- args[[i + 1]]
    i <- i + 2
  }
  defaults
}

opt <- parse_cli(cli_defaults)

assert_output_guardrails(opt$output_dir)

required_inputs <- c(
  overall_implied_csv = opt$overall_implied_csv,
  overall_summary_csv = opt$overall_summary_csv,
  subtype_implied_csv = opt$subtype_implied_csv,
  subtype_summary_csv = opt$subtype_summary_csv,
  subtype_terms_csv = opt$subtype_terms_csv
)

missing_flags <- names(required_inputs)[!nzchar(required_inputs)]
if (length(missing_flags)) {
  stop("Missing required input argument(s): ", paste(paste0("--", gsub("_", "-", missing_flags)), collapse = ", "))
}

if (!dir.exists(opt$output_dir)) dir.create(opt$output_dir, recursive = TRUE)
if (!dir.exists(file.path(opt$output_dir, "tables"))) dir.create(file.path(opt$output_dir, "tables"))
if (!dir.exists(file.path(opt$output_dir, "logs"))) dir.create(file.path(opt$output_dir, "logs"))

overall_implied <- readr::read_csv(opt$overall_implied_csv, show_col_types = FALSE) %>%
  mutate(level = "Overall", subtype = "ALL")

overall_summary <- readr::read_csv(opt$overall_summary_csv, show_col_types = FALSE) %>%
  select(feature, endpoint, model_rows, model_events)

overall_implied <- overall_implied %>%
  left_join(overall_summary, by = c("feature", "endpoint"))

subtype_implied <- readr::read_csv(opt$subtype_implied_csv, show_col_types = FALSE) %>%
  mutate(level = "Subtype")

subtype_summary <- readr::read_csv(opt$subtype_summary_csv, show_col_types = FALSE) %>%
  select(subtype, feature, endpoint, model_rows, model_events)

subtype_implied <- subtype_implied %>%
  left_join(subtype_summary, by = c("subtype", "feature", "endpoint"))

consistency <- bind_rows(
  overall_implied %>%
    transmute(level, subtype, feature, endpoint, dataset, window, HR, CI_low, CI_high, p_value, p_value_bh, model_rows, model_events),
  subtype_implied %>%
    transmute(level, subtype, feature, endpoint, dataset, window, HR, CI_low, CI_high, p_value, p_value_bh, model_rows, model_events)
) %>%
  arrange(endpoint, window, dataset, level, subtype)

readr::write_csv(consistency, file.path(opt$output_dir, "tables", "signature_overall_vs_subtype_consistency.csv"))

fmt <- function(x, d = 3) ifelse(is.na(x), "NA", formatC(x, format = "g", digits = d))
md <- consistency %>%
  mutate(
    `HR (High vs Low)` = fmt(HR, 4),
    `95% CI` = paste0("[", fmt(CI_low, 4), ", ", fmt(CI_high, 4), "]"),
    `p` = fmt(p_value, 3),
    `q (BH)` = fmt(p_value_bh, 3),
    `Rows` = as.character(model_rows),
    `Events` = as.character(model_events)
  ) %>%
  select(level, subtype, endpoint, dataset, window, `Rows`, `Events`, `HR (High vs Low)`, `95% CI`, `p`, `q (BH)`)

md_header <- c(
  "| Level | Subtype | Endpoint | Dataset | Window | Rows | Events | HR (High vs Low) | 95% CI | p | q (BH) |",
  "| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |"
)
md_rows <- apply(md, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
writeLines(c(md_header, md_rows), con = file.path(opt$output_dir, "tables", "signature_overall_vs_subtype_consistency.md"))

subtype_terms <- readr::read_csv(opt$subtype_terms_csv, show_col_types = FALSE) %>%
  filter(feature == opt$signature_feature_name, grepl(opt$subtype_term_pattern, term)) %>%
  mutate(reference = "Reference subtype is first factor level in model data") %>%
  arrange(endpoint, term)

readr::write_csv(subtype_terms, file.path(opt$output_dir, "tables", "overall_model_subtype_covariate_effects.csv"))

subtype_terms_md <- subtype_terms %>%
  mutate(
    HR = fmt(HR, 4),
    `95% CI` = paste0("[", fmt(CI_low, 4), ", ", fmt(CI_high, 4), "]"),
    p = fmt(p_value, 3),
    `q (BH)` = fmt(p_value_bh, 3)
  ) %>%
  select(endpoint, term, HR, `95% CI`, p, `q (BH)`, covariates)

terms_header <- c(
  "| Endpoint | Subtype term (vs reference) | HR | 95% CI | p | q (BH) | Covariates |",
  "| --- | --- | --- | --- | --- | --- | --- |"
)
terms_rows <- apply(subtype_terms_md, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
writeLines(c(terms_header, terms_rows), con = file.path(opt$output_dir, "tables", "overall_model_subtype_covariate_effects.md"))

write_run_metadata_guardrail(
  output_dir = opt$output_dir,
  purpose = "Compare overall vs subtype-specific signature piecewise results and extract subtype covariate effects from the overall model.",
  driver = "drivers/run_signature_consistency_summary.R",
  inputs = as.list(required_inputs),
  parameters = list(
    signature_feature_name = opt$signature_feature_name,
    subtype_term_pattern = opt$subtype_term_pattern
  ),
  scope = "exploratory"
)

message("Signature consistency summary written to: ", opt$output_dir)
