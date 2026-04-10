# Run_OS_PFI_Quartile_Survival_BySubtype_HTML.R
# Updated subtype-stratified survival analysis (quartiles) with full feature parity:
# - Individual & adjusted Cox (genes + signature)
# - LASSO
# - Per-gene KM + combined multi-gene panels (standard + early/late)
# - Time-dependent PH handling (early/late) for signature & genes
# - Markdown + standalone HTML reports per subtype & endpoint
# - Aggregate cross-subtype summary (CSV + MD + HTML)
#
# Requirements: GeneExpressionBasedSurvivalAnalys.r (updated master function)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(survival)
  library(survminer)
  library(ggplot2)
})

source("src/core/GeneExpressionBasedSurvivalAnalys.r")
source("src/utilities/run_provenance_guardrails.R")

analyze_survival_expression <- get("analyze_survival_expression", mode = "function")
assert_output_guardrails <- get("assert_output_guardrails", mode = "function")
write_run_metadata_guardrail <- get("write_run_metadata_guardrail", mode = "function")
normalize_path_safe <- get("normalize_path_safe", mode = "function")
is_legacy_output_root <- get("is_legacy_output_root", mode = "function")

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------
cli_defaults <- list(
  expr_file = "data/IODC-0001-TCGA-BRCA.normalized_counts.tsv",
  meta_file = "data/Survival_subtype_joined_metadata.csv",
  output_root = "outputs/quartile_analysis/by_subtype/Output_Quartile_Survival_Subtypes",
  primary_filter = "tcga_primary",
  subtype_column = "Subtype",
  allow_legacy_output = "false"
)

parse_cli <- function(defaults) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) {
    return(defaults)
  }
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop("Arguments must be supplied as --flag value pairs. Offending token: ", key)
    }
    key_clean <- gsub("-", "_", substring(key, 3))
    if (!key_clean %in% names(defaults)) {
      stop("Unknown argument: ", key)
    }
    if (i == length(args)) {
      stop("No value supplied for argument ", key)
    }
    defaults[[key_clean]] <- args[[i + 1]]
    i <- i + 2
  }
  defaults
}

opt <- parse_cli(cli_defaults)
opt$primary_filter <- tolower(opt$primary_filter)
if (!opt$primary_filter %in% c("tcga_primary", "none")) {
  stop("--primary-filter must be one of: tcga_primary, none")
}

expr_file <- opt$expr_file
meta_file <- opt$meta_file
output_root <- opt$output_root
allow_legacy_output <- tolower(as.character(opt$allow_legacy_output)) %in% c("1", "true", "yes", "y") || identical(Sys.getenv("SURVIVAL_ALLOW_LEGACY_OUTPUT"), "1")
assert_output_guardrails(output_root, allow_legacy_output = allow_legacy_output)
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

signature_genes <- c("GENE_A", "GENE_B", "GENE_C", "GENE_D")
base_covariates <- c("Diagnosis_Age", "ajcc_pathologic_tumor_stage", "race_collapsed") # exclude Subtype inside strata

min_samples <- 30
min_events <- 5
cut_months <- 60

write_run_metadata <- function(output_dir, expr_path, meta_path, n_before, n_after, covariates_used, subtype_column) {
  write_run_metadata_guardrail(
    output_dir = output_dir,
    purpose = "Run OS and PFI subtype-stratified survival analyses using a placeholder 4-gene signature with quartile grouping.",
    driver = "src/drivers/Run_OS_PFI_Quartile_Survival_BySubtype_HTML.R",
    inputs = list(Expression_TSV = expr_path, Metadata_CSV = meta_path),
    parameters = list(
      primary_filter = opt$primary_filter,
      metadata_rows_pre_filter = n_before,
      metadata_rows_post_filter = n_after,
      subtype_column = subtype_column,
      covariates_used = ifelse(length(covariates_used), paste(covariates_used, collapse = ", "), "<none>"),
      min_samples = min_samples,
      min_events = min_events,
      cut_months = cut_months,
      output_root = normalize_path_safe(output_dir),
      allow_legacy_output = allow_legacy_output
    ),
    notes = if (is_legacy_output_root(output_dir)) "Legacy Output_* root used via explicit override (SURVIVAL_ALLOW_LEGACY_OUTPUT=1 or --allow-legacy-output true)." else character(0)
  )
}

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
read_expression <- function(path) {
  expr <- readr::read_tsv(path, show_col_types = FALSE)
  if (!"Gene" %in% colnames(expr)) colnames(expr)[1] <- "Gene"
  expr
}

read_metadata <- function(path) {
  meta <- readr::read_csv(path, show_col_types = FALSE)
  colnames(meta) <- gsub("\n", "", colnames(meta))
  meta
}

md_table <- function(df, max_rows = NULL, digits = 4) {
  if (is.null(df) || length(df) == 0) {
    return(character(0))
  }
  if (!is.data.frame(df)) df <- as.data.frame(df)
  if (nrow(df) == 0) {
    return(character(0))
  }
  if (!is.null(max_rows)) df <- head(df, max_rows)
  df_fmt <- df
  num_cols <- vapply(df_fmt, is.numeric, logical(1))
  df_fmt[num_cols] <- lapply(df_fmt[num_cols], function(x) format(signif(x, digits), trim = TRUE, scientific = FALSE))
  header <- paste(colnames(df_fmt), collapse = " | ")
  sep <- paste(vapply(colnames(df_fmt), function(x) "---", character(1)), collapse = " | ")
  rows <- apply(df_fmt, 1, function(r) paste(r, collapse = " | "))
  c(paste0("| ", header, " |"), paste0("| ", sep, " |"), paste0("| ", rows, " |"), "")
}

# Lightweight MD -> HTML (table + heading + paragraph). Avoid regex for initial pipe detection.
md_to_html <- function(lines, page_title) {
  html <- c()
  i <- 1L
  n <- length(lines)
  esc <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  # Process markdown formatting: **bold**, *italic*, `code`
  fmt <- function(x) {
    x <- esc(x)
    x <- gsub("\\*\\*([^*]+)\\*\\*", "<strong>\\1</strong>", x) # **bold**
    x <- gsub("\\*([^*]+)\\*", "<em>\\1</em>", x) # *italic*
    x <- gsub("`([^`]+)`", "<code>\\1</code>", x) # `code`
    x
  }
  while (i <= n) {
    ln <- lines[i]
    if (nzchar(ln) && substr(ln, 1, 1) == "|") {
      block <- c()
      while (i <= n && nzchar(lines[i]) && substr(lines[i], 1, 1) == "|") {
        block <- c(block, lines[i])
        i <- i + 1L
      }
      # remove alignment row if present (| --- | --- |)
      if (length(block) >= 2 && grepl("^\\| *-+", block[2])) {
        header_line <- block[1]
        body <- block[-c(1, 2)]
      } else {
        header_line <- block[1]
        body <- block[-1]
      }
      sf <- function(line) {
        line <- sub("^\\|", "", line)
        line <- sub("\\n*$", "", line)
        line <- sub("\\r*$", "", line)
        line <- sub("\\| *$", "", line)
        trimws(strsplit(line, "\\|")[[1]])
      }
      heads <- sf(header_line)
      html <- c(html, '<table class="report-table">', "<thead><tr>", paste0("<th>", esc(heads), "</th>", collapse = ""), "</tr></thead>", "<tbody>")
      if (length(body)) {
        for (bl in body) {
          if (!nzchar(gsub("[\n[:space:]]", "", bl))) next
          cells <- sf(bl)
          if (length(cells) == length(heads)) html <- c(html, "<tr>", paste0("<td>", esc(cells), "</td>", collapse = ""), "</tr>")
        }
      }
      html <- c(html, "</tbody></table>")
      next
    }
    if (grepl("^### ", ln)) {
      html <- c(html, paste0("<h3>", fmt(sub("^### ", "", ln)), "</h3>"))
    } else if (grepl("^## ", ln)) {
      html <- c(html, paste0("<h2>", fmt(sub("^## ", "", ln)), "</h2>"))
    } else if (grepl("^# ", ln)) {
      html <- c(html, paste0("<h1>", fmt(sub("^# ", "", ln)), "</h1>"))
    } else if (grepl("^- ", ln)) {
      html <- c(html, paste0("<li>", fmt(sub("^- ", "", ln)), "</li>"))
    } # List items
    else if (nzchar(ln)) html <- c(html, paste0("<p>", fmt(ln), "</p>"))
    i <- i + 1L
  }
  doc <- c(
    "<!DOCTYPE html>", '<html lang="en"><head><meta charset="UTF-8"/>',
    sprintf("<title>%s</title>", esc(page_title)),
    "<style>body{font-family:Arial,Helvetica,sans-serif;margin:1.5rem;line-height:1.35;font-size:14px;} h1,h2,h3{font-weight:600;margin-top:1.7em;} table.report-table{border-collapse:collapse;margin:1em 0;font-size:13px;} table.report-table th,table.report-table td{border:1px solid #999;padding:4px 6px;text-align:left;} table.report-table th{background:#f2f2f2;} p{margin:0.45em 0;} .note{color:#555;font-size:12px;margin-top:3em;border-top:1px solid #ddd;padding-top:0.5em;} code{background:#f5f5f5;padding:2px 4px;border-radius:3px;}</style>",
    "</head><body>", html,
    sprintf('<div class="note">Generated on %s</div>', format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "</body></html>"
  )
  doc
}

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------
expr <- read_expression(expr_file)
meta <- read_metadata(meta_file)
stopifnot("sampnames" %in% colnames(meta))
stopifnot(all(c("OS", "OS_time", "PFI", "PFI_time") %in% colnames(meta)))
if (!opt$subtype_column %in% colnames(meta)) stop("Subtype column missing in metadata: ", opt$subtype_column)

meta_primary <- meta
if (identical(opt$primary_filter, "tcga_primary")) {
  meta_primary <- meta_primary %>% filter(grepl("_01", sampnames))
}
if (!nrow(meta_primary)) stop("No metadata rows remain after applying primary_filter=", opt$primary_filter, ".")

covariates_used <- base_covariates[base_covariates %in% colnames(meta_primary)]
write_run_metadata(
  output_dir = output_root,
  expr_path = expr_file,
  meta_path = meta_file,
  n_before = nrow(meta),
  n_after = nrow(meta_primary),
  covariates_used = covariates_used,
  subtype_column = opt$subtype_column
)

subtypes <- sort(unique(na.omit(meta_primary[[opt$subtype_column]])))
if (!length(subtypes)) stop("No subtype levels found")
message("Subtypes detected: ", paste(subtypes, collapse = ", "))

# Aggregate rows across subtype & endpoint
agg_rows <- list()
gene_agg_rows <- list() # For individual gene summaries

# Helper function to format numbers for markdown tables
format_table_number <- function(x) {
  if (is.na(x) || is.infinite(x)) {
    return(as.character(x))
  }
  if (abs(x) > 1000) {
    # Use scientific notation for large numbers
    return(sprintf("%.2e", x))
  }
  return(as.character(round(x, 4)))
}

# -----------------------------------------------------------------------------
# Core per-subtype runner
# -----------------------------------------------------------------------------
run_endpoint_subtype <- function(meta_sub, endpoint_event, endpoint_time, endpoint_label, subtype_label) {
  safe_sub <- gsub("[^A-Za-z0-9]+", "_", subtype_label)
  out_dir <- file.path(output_root, safe_sub, endpoint_label)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  # Feasibility counts
  events_vec <- meta_sub[[endpoint_event]]
  time_vec <- meta_sub[[endpoint_time]]
  events_num <- suppressWarnings(as.numeric(events_vec))
  if (all(is.na(events_num))) {
    events_num <- ifelse(tolower(as.character(events_vec)) == "dead", 1,
      ifelse(tolower(as.character(events_vec)) == "alive", 0, NA)
    )
  }
  time_num <- suppressWarnings(as.numeric(time_vec))
  n_events <- sum(events_num == 1 & !is.na(time_num) & !is.na(events_num))
  n_samples <- sum(!is.na(time_num) & !is.na(events_num))

  if (n_samples < min_samples || n_events < min_events) {
    msg <- sprintf("Skipping %s / %s (n=%d, events=%d) below thresholds (min_samples=%d, min_events=%d)", subtype_label, endpoint_label, n_samples, n_events, min_samples, min_events)
    message(msg)
    writeLines(msg, con = file.path(out_dir, paste0(endpoint_label, "_Summary.txt")))
    return(NULL)
  }

  lasso_ok <- requireNamespace("glmnet", quietly = TRUE)

  res <- analyze_survival_expression(
    expression_data = expr,
    metadata = meta_sub,
    sample_col = "sampnames",
    gene_col = "Gene",
    genes = signature_genes,
    event_col = endpoint_event,
    time_col = endpoint_time,
    clinical_covariates = covariates_used,
    split_method = "quartiles",
    include_groups = c("low", "mid", "high"),
    steps_to_run = c("individual", "individual+covariates", "signature", "signature+covariates", "lasso", "km_individual"),
    run_lasso = lasso_ok,
    run_ph_checks = TRUE,
    run_time_dependent = TRUE,
    cut_months = cut_months,
    save_plots = TRUE,
    export_csv = TRUE,
    output_dir = out_dir,
    plot_params = list(plot_title = paste0(endpoint_label, " (", subtype_label, ") Survival")),
    verbose = TRUE
  )

  saveRDS(res, file = file.path(out_dir, paste0(endpoint_label, "_results.rds")))

  # ---------------- Plain Text Summary ----------------
  summary_lines <- c(
    paste0("Subtype: ", subtype_label),
    paste0("Endpoint: ", endpoint_label),
    paste0("Samples analyzed: ", ifelse(!is.null(res$sample_groups), length(unique(na.omit(unlist(res$sample_groups)))), NA))
  )
  if (!is.null(res$results_univariate)) summary_lines <- c(summary_lines, "\nUnivariate Cox (first 10 rows):", utils::capture.output(print(head(res$results_univariate, 10))))
  if (!is.null(res$results_adjusted)) summary_lines <- c(summary_lines, "\nAdjusted Cox (first 10 rows):", utils::capture.output(print(head(res$results_adjusted, 10))))
  if (!is.null(res$results_signature)) {
    summary_lines <- c(summary_lines, "\nSignature Cox (unadjusted):", utils::capture.output(print(res$results_signature$unadjusted)))
    if (!is.null(res$results_signature$adjusted)) summary_lines <- c(summary_lines, "\nSignature Cox (adjusted):", utils::capture.output(print(res$results_signature$adjusted)))
  }
  if (!is.null(res$ph_diagnostics)) {
    summary_lines <- c(summary_lines, "\nProportional Hazards Test:", utils::capture.output(print(res$ph_diagnostics$ph_test)))
    if (!is.null(res$ph_diagnostics$time_dependent)) summary_lines <- c(summary_lines, "\nTime-Dependent Model Summary:", utils::capture.output(print(res$ph_diagnostics$time_dependent)))
  }
  writeLines(summary_lines, con = file.path(out_dir, paste0(endpoint_label, "_Summary.txt")))

  # ---------------- Markdown + HTML Report ----------------
  md <- c(paste0("# ", endpoint_label, " Survival (", subtype_label, ") Quartile Analysis"), "")
  md <- c(md, "## Samples", "", paste0("Subtype: ", subtype_label), paste0("Total samples (metadata): ", nrow(meta_sub)), paste0("Analyzed samples (with survival & expression): ", ifelse(!is.null(res$sample_groups), length(unique(na.omit(unlist(res$sample_groups)))), NA)), "")
  if (!is.null(res$results_univariate)) md <- c(md, "## Univariate Cox (top 10)", "", md_table(res$results_univariate, 10))
  if (!is.null(res$results_adjusted)) md <- c(md, "## Adjusted Cox (top 10)", "", md_table(res$results_adjusted, 10))
  if (!is.null(res$results_signature)) {
    # Check if we have any valid signature results before adding section header
    has_unadj <- !is.null(res$results_signature$unadjusted) && inherits(res$results_signature$unadjusted, "summary.coxph")
    has_adj <- !is.null(res$results_signature$adjusted) && inherits(res$results_signature$adjusted, "summary.coxph")

    if (has_unadj || has_adj) {
      md <- c(md, "## Signature Cox", "")

      if (has_unadj) {
        su <- res$results_signature$unadjusted
        co <- as.data.frame(su$coefficients)
        ci <- as.data.frame(su$conf.int)
        tab <- data.frame(
          Term = rownames(co),
          HR = round(ci[, "exp(coef)"], 4),
          CI_Lower = round(ci[, "lower .95"], 4),
          CI_Upper = round(ci[, "upper .95"], 4),
          P_Value = round(co[, "Pr(>|z|)"], 4),
          stringsAsFactors = FALSE
        )
        md <- c(md, "### Unadjusted", "", md_table(tab))
      } else {
        md <- c(md, "### Unadjusted", "", "No unadjusted signature results available.", "")
      }

      if (has_adj) {
        sa <- res$results_signature$adjusted
        co <- as.data.frame(sa$coefficients)
        ci <- as.data.frame(sa$conf.int)
        tab2 <- data.frame(
          Term = rownames(co),
          HR = sapply(ci[, "exp(coef)"], format_table_number),
          CI_Lower = sapply(ci[, "lower .95"], format_table_number),
          CI_Upper = sapply(ci[, "upper .95"], format_table_number),
          P_Value = sapply(co[, "Pr(>|z|)"], format_table_number),
          stringsAsFactors = FALSE
        )
        md <- c(md, "### Adjusted", "", md_table(tab2))
      } else {
        md <- c(md, "### Adjusted", "", "No adjusted signature results available (insufficient covariates or model failure).", "")
      }
    }
  }

  # Add early window Cox results
  if (!is.null(res$results_signature_early)) {
    has_early_unadj <- !is.null(res$results_signature_early$unadjusted) && inherits(res$results_signature_early$unadjusted, "summary.coxph")
    has_early_adj <- !is.null(res$results_signature_early$adjusted) && inherits(res$results_signature_early$adjusted, "summary.coxph")

    if (has_early_unadj || has_early_adj) {
      md <- c(md, sprintf("## Signature Cox - Early Window (≤60 months, n=%d)", res$results_signature_early$n_samples), "")

      if (has_early_unadj) {
        su <- res$results_signature_early$unadjusted
        co <- as.data.frame(su$coefficients)
        ci <- as.data.frame(su$conf.int)
        tab <- data.frame(
          Term = rownames(co),
          HR = sapply(ci[, "exp(coef)"], format_table_number),
          CI_Lower = sapply(ci[, "lower .95"], format_table_number),
          CI_Upper = sapply(ci[, "upper .95"], format_table_number),
          P_Value = sapply(co[, "Pr(>|z|)"], format_table_number),
          stringsAsFactors = FALSE
        )
        md <- c(md, "### Unadjusted", "", md_table(tab))
      }

      if (has_early_adj) {
        sa <- res$results_signature_early$adjusted
        co <- as.data.frame(sa$coefficients)
        ci <- as.data.frame(sa$conf.int)
        tab2 <- data.frame(
          Term = rownames(co),
          HR = sapply(ci[, "exp(coef)"], format_table_number),
          CI_Lower = sapply(ci[, "lower .95"], format_table_number),
          CI_Upper = sapply(ci[, "upper .95"], format_table_number),
          P_Value = sapply(co[, "Pr(>|z|)"], format_table_number),
          stringsAsFactors = FALSE
        )
        md <- c(md, "### Adjusted", "", md_table(tab2))
      }
    }
  }

  # Add late window Cox results
  if (!is.null(res$results_signature_late)) {
    has_late_unadj <- !is.null(res$results_signature_late$unadjusted) && inherits(res$results_signature_late$unadjusted, "summary.coxph")
    has_late_adj <- !is.null(res$results_signature_late$adjusted) && inherits(res$results_signature_late$adjusted, "summary.coxph")

    if (has_late_unadj || has_late_adj) {
      md <- c(md, sprintf("## Signature Cox - Late Window (>60 months, n=%d)", res$results_signature_late$n_samples), "")

      if (has_late_unadj) {
        su <- res$results_signature_late$unadjusted
        co <- as.data.frame(su$coefficients)
        ci <- as.data.frame(su$conf.int)
        tab <- data.frame(
          Term = rownames(co),
          HR = sapply(ci[, "exp(coef)"], format_table_number),
          CI_Lower = sapply(ci[, "lower .95"], format_table_number),
          CI_Upper = sapply(ci[, "upper .95"], format_table_number),
          P_Value = sapply(co[, "Pr(>|z|)"], format_table_number),
          stringsAsFactors = FALSE
        )
        md <- c(md, "### Unadjusted", "", md_table(tab))
      }

      if (has_late_adj) {
        sa <- res$results_signature_late$adjusted
        co <- as.data.frame(sa$coefficients)
        ci <- as.data.frame(sa$conf.int)
        tab2 <- data.frame(
          Term = rownames(co),
          HR = sapply(ci[, "exp(coef)"], format_table_number),
          CI_Lower = sapply(ci[, "lower .95"], format_table_number),
          CI_Upper = sapply(ci[, "upper .95"], format_table_number),
          P_Value = sapply(co[, "Pr(>|z|)"], format_table_number),
          stringsAsFactors = FALSE
        )
        md <- c(md, "### Adjusted", "", md_table(tab2))
      }
    }
  }

  if (!is.null(res$ph_diagnostics)) {
    md <- c(md, "## Proportional Hazards Test", "")
    if (!is.null(res$ph_diagnostics$ph_test) && inherits(res$ph_diagnostics$ph_test, "cox.zph")) {
      ph <- res$ph_diagnostics$ph_test
      ph_tab <- as.data.frame(ph$table)
      ph_tab$Term <- rownames(ph_tab)
      ph_tab <- ph_tab[, c("Term", setdiff(colnames(ph_tab), "Term"))]
      md <- c(md, md_table(ph_tab))
    }
    if (!is.null(res$ph_diagnostics$time_dependent) && inherits(res$ph_diagnostics$time_dependent, "summary.coxph")) {
      td <- res$ph_diagnostics$time_dependent
      co <- as.data.frame(td$coefficients)
      ci <- as.data.frame(td$conf.int)
      td_tab <- cbind(Term = rownames(co), HR = ci$`exp(coef)`, CI_Lower = ci$`lower .95`, CI_Upper = ci$`upper .95`, P_Value = co$`Pr(>|z|)`)
      md <- c(md, "### Time-Dependent Model Coefficients", "", md_table(td_tab))
    }
  }
  if (!is.null(res$lasso_selected)) {
    md <- c(md, "## LASSO Selection", "")
    if (length(res$lasso_selected) == 0) md <- c(md, "No genes selected.", "") else md <- c(md, paste0("Selected genes: ", paste(res$lasso_selected, collapse = ", ")), "")
  }
  if (!is.null(res$sample_groups)) {
    grp_counts <- vapply(res$sample_groups, function(col) sum(!is.na(col)), integer(1))
    grp_df <- data.frame(Group = names(grp_counts), N = as.integer(grp_counts))
    md <- c(md, "## Group Sizes", "", md_table(grp_df))
  }

  md_path <- file.path(out_dir, paste0(endpoint_label, "_Summary.md"))
  writeLines(md, md_path)
  html_doc <- md_to_html(md, page_title = paste0(endpoint_label, " (", subtype_label, ") Survival Quartile Analysis"))
  writeLines(html_doc, file.path(out_dir, paste0(endpoint_label, "_Summary.html")))

  # Aggregate metrics across subtypes
  sig_unadj <- res$results_signature$unadjusted
  hr_high <- p_high <- ci_low <- ci_high <- NA_real_
  adj_hr_high <- adj_p_high <- adj_ci_low <- adj_ci_high <- NA_real_
  events_total <- NA_real_
  if (!is.null(sig_unadj) && "clusterC" %in% rownames(sig_unadj$coefficients)) {
    idx <- which(rownames(sig_unadj$coefficients) == "clusterC")[1]
    hr_high <- sig_unadj$conf.int[idx, "exp(coef)"]
    ci_low <- sig_unadj$conf.int[idx, "lower .95"]
    ci_high <- sig_unadj$conf.int[idx, "upper .95"]
    p_high <- sig_unadj$coefficients[idx, "Pr(>|z|)"]
    events_total <- sig_unadj$nevent
  }
  sig_adj <- res$results_signature$adjusted
  if (!is.null(sig_adj) && "clusterC" %in% rownames(sig_adj$coefficients)) {
    idx2 <- which(rownames(sig_adj$coefficients) == "clusterC")[1]
    adj_hr_high <- sig_adj$conf.int[idx2, "exp(coef)"]
    adj_ci_low <- sig_adj$conf.int[idx2, "lower .95"]
    adj_ci_high <- sig_adj$conf.int[idx2, "upper .95"]
    adj_p_high <- sig_adj$coefficients[idx2, "Pr(>|z|)"]
  }
  ph_global_p <- if (!is.null(res$ph_diagnostics) && !is.null(res$ph_diagnostics$p_global)) res$ph_diagnostics$p_global else NA_real_
  ph_signature_p <- if (!is.null(res$ph_diagnostics) && !is.null(res$ph_diagnostics$p_signature)) res$ph_diagnostics$p_signature else NA_real_

  # Extract time-dependent coefficients for early and late windows
  td_early_coef <- td_early_p <- td_late_coef <- td_late_p <- NA_real_
  if (!is.null(res$ph_diagnostics)) {
    td_sum <- NULL
    td_source <- NULL
    if (!is.null(res$ph_diagnostics$grouped) && !is.null(res$ph_diagnostics$grouped$time_dependent)) {
      td_sum <- res$ph_diagnostics$grouped$time_dependent
      td_source <- "grouped"
    } else if (!is.null(res$ph_diagnostics$time_dependent)) {
      td_sum <- res$ph_diagnostics$time_dependent
      td_source <- "continuous"
    }

    if (!is.null(td_sum) && inherits(td_sum, "summary.coxph")) {
      rn <- rownames(td_sum$coefficients)
      if (td_source == "grouped") {
        pick_td_row <- function(phase) {
          candidates <- c(
            paste0("cluster_C_", phase),
            paste0("clusterHigh_", phase),
            paste0("cluster_B_", phase)
          )
          hit <- candidates[candidates %in% rn]
          if (!length(hit)) {
            alt <- rn[grepl(paste0("_", phase, "$"), rn) & !grepl("cluster_A_", rn)]
            if (length(alt)) hit <- alt[1]
          }
          if (length(hit)) hit[1] else NULL
        }
        row_early <- pick_td_row("early")
        row_late <- pick_td_row("late")
        if (!is.null(row_early)) {
          td_early_coef <- td_sum$coefficients[row_early, "coef"]
          td_early_p <- td_sum$coefficients[row_early, "Pr(>|z|)"]
        }
        if (!is.null(row_late)) {
          td_late_coef <- td_sum$coefficients[row_late, "coef"]
          td_late_p <- td_sum$coefficients[row_late, "Pr(>|z|)"]
        }
      } else {
        wearly <- which(grepl("sig_early", rn))
        wlate <- which(grepl("sig_late", rn))
        if (length(wearly)) {
          td_early_coef <- td_sum$coefficients[wearly[1], "coef"]
          td_early_p <- td_sum$coefficients[wearly[1], "Pr(>|z|)"]
        }
        if (length(wlate)) {
          td_late_coef <- td_sum$coefficients[wlate[1], "coef"]
          td_late_p <- td_sum$coefficients[wlate[1], "Pr(>|z|)"]
        }
      }
    }
  }

  # Extract early and late window Cox models (High vs Low)
  early_hr <- early_p <- late_hr <- late_p <- NA_real_
  if (!is.null(res$results_signature_early) && !is.null(res$results_signature_early$unadjusted)) {
    early_sum <- res$results_signature_early$unadjusted
    if ("clusterHigh" %in% rownames(early_sum$coefficients)) {
      idx_early <- which(rownames(early_sum$coefficients) == "clusterHigh")[1]
      early_hr <- early_sum$conf.int[idx_early, "exp(coef)"]
      early_p <- early_sum$coefficients[idx_early, "Pr(>|z|)"]
    }
  }
  if (!is.null(res$results_signature_late) && !is.null(res$results_signature_late$unadjusted)) {
    late_sum <- res$results_signature_late$unadjusted
    if ("clusterHigh" %in% rownames(late_sum$coefficients)) {
      idx_late <- which(rownames(late_sum$coefficients) == "clusterHigh")[1]
      late_hr <- late_sum$conf.int[idx_late, "exp(coef)"]
      late_p <- late_sum$coefficients[idx_late, "Pr(>|z|)"]
    }
  }

  # Round all numeric values to 4 decimal places
  agg_rows[[length(agg_rows) + 1]] <<- data.frame(
    Subtype = subtype_label,
    Endpoint = endpoint_label,
    N = n_samples,
    Events = events_total,
    HR_High_vs_Low = round(hr_high, 4),
    CI_Lower = round(ci_low, 4),
    CI_Upper = round(ci_high, 4),
    P_High_vs_Low = round(p_high, 4),
    Adj_HR_High_vs_Low = round(adj_hr_high, 4),
    Adj_CI_Lower = round(adj_ci_low, 4),
    Adj_CI_Upper = round(adj_ci_high, 4),
    Adj_P_High_vs_Low = round(adj_p_high, 4),
    PH_Global_p = round(ph_global_p, 4),
    PH_Signature_p = round(ph_signature_p, 4),
    Early_HR = round(early_hr, 4),
    Early_p = round(early_p, 4),
    Late_HR = round(late_hr, 4),
    Late_p = round(late_p, 4),
    TD_Early_coef = round(td_early_coef, 4),
    TD_Early_p = round(td_early_p, 4),
    TD_Late_coef = round(td_late_coef, 4),
    TD_Late_p = round(td_late_p, 4),
    stringsAsFactors = FALSE
  )

  # Extract individual gene statistics from the analysis results (matching plots exactly)
  for (gene in signature_genes) {
    # Check if plot results exist for this gene
    if (!(gene %in% names(res$plots))) next

    plot_obj <- res$plots[[gene]]
    if (is.null(plot_obj)) next

    # Extract unadjusted Cox results (quartile grouping)
    hr_high <- ci_low <- ci_high <- p_high <- NA_real_
    if (!is.null(plot_obj$cox_unadj) && inherits(plot_obj$cox_unadj, "summary.coxph")) {
      cox_unadj <- plot_obj$cox_unadj
      # Look for clusterC (High vs Low with Low as reference)
      idx_high <- which(rownames(cox_unadj$coefficients) == "clusterC")
      if (length(idx_high) > 0) {
        hr_high <- cox_unadj$conf.int[idx_high[1], "exp(coef)"]
        ci_low <- cox_unadj$conf.int[idx_high[1], "lower .95"]
        ci_high <- cox_unadj$conf.int[idx_high[1], "upper .95"]
        p_high <- cox_unadj$coefficients[idx_high[1], "Pr(>|z|)"]
      }
    }

    # Extract adjusted Cox results
    adj_hr_high <- adj_ci_low <- adj_ci_high <- adj_p_high <- NA_real_
    if (!is.null(plot_obj$cox_adj) && inherits(plot_obj$cox_adj, "summary.coxph")) {
      cox_adj <- plot_obj$cox_adj
      idx_high_adj <- which(rownames(cox_adj$coefficients) == "clusterC")
      if (length(idx_high_adj) > 0) {
        adj_hr_high <- cox_adj$conf.int[idx_high_adj[1], "exp(coef)"]
        adj_ci_low <- cox_adj$conf.int[idx_high_adj[1], "lower .95"]
        adj_ci_high <- cox_adj$conf.int[idx_high_adj[1], "upper .95"]
        adj_p_high <- cox_adj$coefficients[idx_high_adj[1], "Pr(>|z|)"]
      }
    }

    # Get sample size and event count from the unadjusted model
    n_samples <- events_count <- NA_integer_
    if (!is.null(plot_obj$cox_unadj)) {
      n_samples <- plot_obj$cox_unadj$n
      events_count <- plot_obj$cox_unadj$nevent
    }

    # PH tests - extract from plot object if available
    ph_global_p <- ph_gene_p <- NA_real_
    if (!is.null(plot_obj$ph_test) && !inherits(plot_obj$ph_test, "try-error")) {
      ph_global_p <- plot_obj$ph_test$table["GLOBAL", "p"]
      if ("cluster" %in% rownames(plot_obj$ph_test$table)) {
        ph_gene_p <- plot_obj$ph_test$table["cluster", "p"]
      }
    }

    # Time-dependent model coefficients
    td_early_coef <- td_early_p <- td_late_coef <- td_late_p <- NA_real_
    if (!is.null(plot_obj$td_model) && inherits(plot_obj$td_model, "summary.coxph")) {
      coef_names <- rownames(plot_obj$td_model$coefficients)
      # Look for cluster_C_early and cluster_C_late (High vs Low comparison)
      # These are created as interaction terms: cluster level × time window
      idx_early <- grep("^cluster_C_early$", coef_names)
      idx_late <- grep("^cluster_C_late$", coef_names)

      if (length(idx_early) > 0) {
        td_early_coef <- plot_obj$td_model$coefficients[idx_early[1], "coef"]
        td_early_p <- plot_obj$td_model$coefficients[idx_early[1], "Pr(>|z|)"]
      }
      if (length(idx_late) > 0) {
        td_late_coef <- plot_obj$td_model$coefficients[idx_late[1], "coef"]
        td_late_p <- plot_obj$td_model$coefficients[idx_late[1], "Pr(>|z|)"]
      }
    }

    # Early window results (from gene_early plots if they exist) - using ADJUSTED values
    early_hr <- early_p <- NA_real_
    early_key <- paste0(gene, "_early")
    if (early_key %in% names(res$plots) && !is.null(res$plots[[early_key]])) {
      plot_early <- res$plots[[early_key]]
      # Try adjusted first, fall back to unadjusted if not available
      cox_early <- plot_early$cox_adj
      if (is.null(cox_early) || !inherits(cox_early, "summary.coxph")) {
        cox_early <- plot_early$cox_unadj
      }
      if (!is.null(cox_early) && inherits(cox_early, "summary.coxph")) {
        idx_e <- which(rownames(cox_early$coefficients) == "clusterC")
        if (length(idx_e) > 0) {
          early_hr <- cox_early$conf.int[idx_e[1], "exp(coef)"]
          early_p <- cox_early$coefficients[idx_e[1], "Pr(>|z|)"]
        }
      }
    }

    # Late window results (from gene_late plots if they exist) - using ADJUSTED values
    late_hr <- late_p <- NA_real_
    late_key <- paste0(gene, "_late")
    if (late_key %in% names(res$plots) && !is.null(res$plots[[late_key]])) {
      plot_late <- res$plots[[late_key]]
      # Try adjusted first, fall back to unadjusted if not available
      cox_late <- plot_late$cox_adj
      if (is.null(cox_late) || !inherits(cox_late, "summary.coxph")) {
        cox_late <- plot_late$cox_unadj
      }
      if (!is.null(cox_late) && inherits(cox_late, "summary.coxph")) {
        idx_l <- which(rownames(cox_late$coefficients) == "clusterC")
        if (length(idx_l) > 0) {
          late_hr <- cox_late$conf.int[idx_l[1], "exp(coef)"]
          late_p <- cox_late$coefficients[idx_l[1], "Pr(>|z|)"]
        }
      }
    }

    # Store comprehensive gene results
    gene_agg_rows[[length(gene_agg_rows) + 1]] <<- data.frame(
      Gene = gene,
      Subtype = subtype_label,
      Endpoint = endpoint_label,
      N = n_samples,
      Events = events_count,
      HR_High_vs_Low = round(hr_high, 4),
      CI_Lower = round(ci_low, 4),
      CI_Upper = round(ci_high, 4),
      P_High_vs_Low = round(p_high, 4),
      Adj_HR_High_vs_Low = round(adj_hr_high, 4),
      Adj_CI_Lower = round(adj_ci_low, 4),
      Adj_CI_Upper = round(adj_ci_high, 4),
      Adj_P_High_vs_Low = round(adj_p_high, 4),
      PH_Global_p = round(ph_global_p, 4),
      PH_Gene_p = round(ph_gene_p, 4),
      Early_HR = round(early_hr, 4),
      Early_p = round(early_p, 4),
      Late_HR = round(late_hr, 4),
      Late_p = round(late_p, 4),
      TD_Early_coef = round(td_early_coef, 4),
      TD_Early_p = round(td_early_p, 4),
      TD_Late_coef = round(td_late_coef, 4),
      TD_Late_p = round(td_late_p, 4),
      stringsAsFactors = FALSE
    )
  }

  invisible(res)
}

# -----------------------------------------------------------------------------
# Iterate subtypes
# -----------------------------------------------------------------------------
for (st in subtypes) {
  message("\n==============================")
  message("Processing Subtype: ", st)
  message("==============================")
  meta_sub <- meta_primary %>% filter(Subtype == st)
  run_endpoint_subtype(meta_sub, "OS", "OS_time", "OS", st)
  run_endpoint_subtype(meta_sub, "PFI", "PFI_time", "PFI", st)
}

# -----------------------------------------------------------------------------
# Aggregate summary across all subtype analyses
# -----------------------------------------------------------------------------
if (length(agg_rows)) {
  agg_df <- dplyr::bind_rows(agg_rows)

  # Format large numbers to scientific notation
  format_number <- function(x) {
    if (is.na(x) || is.infinite(x)) {
      return(x)
    }
    if (abs(x) > 1000) {
      return(sprintf("%.2e", x)) # Force scientific notation for large numbers
    }
    return(round(x, 4))
  }

  # Apply formatting to numeric columns
  numeric_cols <- c(
    "HR_High_vs_Low", "CI_Lower", "CI_Upper", "P_High_vs_Low",
    "Adj_HR_High_vs_Low", "Adj_CI_Lower", "Adj_CI_Upper", "Adj_P_High_vs_Low",
    "PH_Global_p", "PH_Signature_p", "Early_HR", "Early_p", "Late_HR", "Late_p",
    "TD_Early_coef", "TD_Early_p", "TD_Late_coef", "TD_Late_p"
  )

  for (col in numeric_cols) {
    if (col %in% names(agg_df)) {
      agg_df[[col]] <- sapply(agg_df[[col]], format_number)
    }
  }

  out_csv <- file.path(output_root, "subtype_signature_summary.csv")
  write.csv(agg_df, out_csv, row.names = FALSE)
  message("Aggregated summary written to ", out_csv)

  # Markdown aggregate with legend
  legend <- c(
    "", "## Column Definitions", "",
    "- **Subtype**: BRCA molecular subtype",
    "- **Endpoint**: Overall Survival (OS) or Progression-Free Interval (PFI)",
    "- **N**: Total number of samples analyzed",
    "- **Events**: Number of events (deaths for OS, progressions for PFI)",
    "- **HR_High_vs_Low**: Hazard ratio comparing highest tertile/quartile vs lowest (unadjusted)",
    "- **CI_Lower/CI_Upper**: 95% confidence interval bounds for HR_High_vs_Low",
    "- **P_High_vs_Low**: P-value for High vs Low comparison (unadjusted)",
    "- **Adj_HR_High_vs_Low**: Adjusted hazard ratio (controlling for age, stage, race)",
    "- **Adj_CI_Lower/Adj_CI_Upper**: 95% CI for adjusted HR",
    "- **Adj_P_High_vs_Low**: P-value for adjusted comparison",
    "- **PH_Global_p**: Global proportional hazards test p-value",
    "- **PH_Signature_p**: Signature-specific proportional hazards test p-value",
    "- **Early_HR**: Hazard ratio (High vs Low) in early window (≤60 months)",
    "- **Early_p**: P-value for early window HR",
    "- **Late_HR**: Hazard ratio (High vs Low) in late window (>60 months)",
    "- **Late_p**: P-value for late window HR",
    "- **TD_Early_coef**: Time-dependent coefficient for early period (from PH violation model)",
    "- **TD_Early_p**: P-value for time-dependent early coefficient",
    "- **TD_Late_coef**: Time-dependent coefficient for late period (from PH violation model)",
    "- **TD_Late_p**: P-value for time-dependent late coefficient",
    "",
    "**Note**: NA values indicate the test was not performed (e.g., TD coefficients only exist when proportional hazards are violated; Late window results require sufficient late events)."
  )

  md <- c(
    "# Subtype Aggregate Signature Summary (Quartile High vs Low)", "",
    "Generated from subtype-specific survival analyses (quartiles).", "",
    "## Table", "", md_table(agg_df), legend
  )
  md_file <- file.path(output_root, "subtype_signature_summary.md")
  writeLines(md, md_file)
  html_doc <- md_to_html(md, page_title = "Subtype Aggregate Signature Summary")
  writeLines(html_doc, file.path(output_root, "subtype_signature_summary.html"))

  # Generate individual gene summary tables
  if (length(gene_agg_rows)) {
    gene_df <- dplyr::bind_rows(gene_agg_rows)

    # Format numbers
    gene_numeric_cols <- c(
      "HR_High_vs_Low", "CI_Lower", "CI_Upper", "P_High_vs_Low",
      "Adj_HR_High_vs_Low", "Adj_CI_Lower", "Adj_CI_Upper", "Adj_P_High_vs_Low",
      "PH_Global_p", "PH_Gene_p", "Early_HR", "Early_p", "Late_HR", "Late_p",
      "TD_Early_coef", "TD_Early_p", "TD_Late_coef", "TD_Late_p"
    )
    for (col in gene_numeric_cols) {
      if (col %in% names(gene_df)) {
        gene_df[[col]] <- sapply(gene_df[[col]], format_number)
      }
    }

    # Create summary for each gene
    for (gene in unique(gene_df$Gene)) {
      gene_subset <- gene_df[gene_df$Gene == gene, ]

      # CSV
      gene_csv <- file.path(output_root, paste0(gene, "_summary.csv"))
      write.csv(gene_subset, gene_csv, row.names = FALSE)

      # Markdown with legend
      gene_legend <- c(
        "", "## Column Definitions", "",
        "- **Gene**: Gene symbol",
        "- **Subtype**: BRCA molecular subtype",
        "- **Endpoint**: Overall Survival (OS) or Progression-Free Interval (PFI)",
        "- **N**: Total number of samples analyzed (all three quartile groups)",
        "- **Events**: Number of events across all quartiles",
        "- **HR_High_vs_Low**: Hazard ratio for High vs Low quartile (from model with all 3 groups)",
        "- **CI_Lower/CI_Upper**: 95% confidence interval bounds for HR_High_vs_Low",
        "- **P_High_vs_Low**: P-value for High vs Low comparison",
        "- **Adj_HR_High_vs_Low**: Adjusted HR (controlling for age, stage, race)",
        "- **Adj_CI_Lower/Adj_CI_Upper**: 95% CI for adjusted HR",
        "- **Adj_P_High_vs_Low**: P-value for adjusted comparison",
        "- **PH_Global_p**: Global proportional hazards test p-value",
        "- **PH_Gene_p**: Gene-specific proportional hazards test p-value",
        "- **Early_HR**: Hazard ratio (High vs Low) in early window (≤60 months, all 3 groups)",
        "- **Early_p**: P-value for early window HR",
        "- **Late_HR**: Hazard ratio (High vs Low) in late window (>60 months, all 3 groups)",
        "- **Late_p**: P-value for late window HR",
        "- **TD_Early_coef**: Time-dependent coefficient for early period (from PH violation model)",
        "- **TD_Early_p**: P-value for time-dependent early coefficient",
        "- **TD_Late_coef**: Time-dependent coefficient for late period (from PH violation model)",
        "- **TD_Late_p**: P-value for time-dependent late coefficient",
        "",
        "**Note**: These use quartile-based High vs Low groupings, matching the aggregate signature analysis methodology."
      )

      gene_md <- c(
        paste0("# ", gene, " Summary Across Subtypes (Quartile Analysis)"), "",
        "Individual gene survival analysis results.", "",
        "## Table", "", md_table(gene_subset), gene_legend
      )
      gene_md_file <- file.path(output_root, paste0(gene, "_summary.md"))
      writeLines(gene_md, gene_md_file)

      # HTML
      gene_html <- md_to_html(gene_md, page_title = paste0(gene, " Summary"))
      writeLines(gene_html, file.path(output_root, paste0(gene, "_summary.html")))
    }

    message("Individual gene summaries written for: ", paste(unique(gene_df$Gene), collapse = ", "))
  }
} else {
  message("No subtype analyses produced aggregate results (all skipped?)")
}

message("Subtype-stratified quartile analyses complete. Output in ", output_root)
