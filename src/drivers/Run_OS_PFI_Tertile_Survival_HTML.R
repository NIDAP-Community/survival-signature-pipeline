# Run_OS_PFI_Tertile_Survival_HTML.R
# Overall BRCA survival analyses for OS and PFI using the 4-gene signature (SAT1, SRM, ODC1, AMD1)
# Using TERTILE grouping strategy

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(survival)
  library(survminer)
  library(ggplot2)
  if (!requireNamespace("cowplot", quietly = TRUE) && !requireNamespace("gridExtra", quietly = TRUE)) {
    message("NOTE: 'cowplot' or 'gridExtra' not available; statistics table layout may fallback.")
  }
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
  expr_file = file.path("data", "IODC-0001-TCGA-BRCA.normalized_counts.tsv"),
  meta_file = file.path("data", "Survival_subtype_joined_metadata.csv"),
  output_root = "outputs/tertile_analysis/overall",
  primary_filter = "tcga_primary",
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

write_run_metadata <- function(output_dir, expr_path, meta_path, n_before, n_after, covariates_used) {
  write_run_metadata_guardrail(
    output_dir = output_dir,
    purpose = "Run OS and PFI survival analyses using a 4-gene signature with tertile grouping.",
    driver = "src/drivers/Run_OS_PFI_Tertile_Survival_HTML.R",
    inputs = list(Expression_TSV = expr_path, Metadata_CSV = meta_path),
    parameters = list(
      primary_filter = opt$primary_filter,
      metadata_rows_pre_filter = n_before,
      metadata_rows_post_filter = n_after,
      covariates_used = ifelse(length(covariates_used), paste(covariates_used, collapse = ", "), "<none>"),
      output_root = normalize_path_safe(output_dir),
      allow_legacy_output = allow_legacy_output
    ),
    notes = c(
      "OS and PFI are expected to be encoded as event (0/1) plus time columns in the metadata.",
      if (is_legacy_output_root(output_dir)) "Legacy Output_* root used via explicit override (SURVIVAL_ALLOW_LEGACY_OUTPUT=1 or --allow-legacy-output true)." else NULL
    )
  )
}

signature_genes <- c("SAT1", "SRM", "ODC1", "AMD1")
base_covariates <- c("Diagnosis_Age", "ajcc_pathologic_tumor_stage", "race_collapsed", "Subtype")

# Helper loaders --------------------------------------------------------------
read_expression <- function(path) {
  expr <- readr::read_tsv(path, show_col_types = FALSE)
  if (!"Gene" %in% colnames(expr)) {
    colnames(expr)[1] <- "Gene"
  }
  expr
}

read_metadata <- function(path) {
  meta <- readr::read_csv(path, show_col_types = FALSE)
  colnames(meta) <- gsub("\n", "", colnames(meta))
  meta
}

expr <- read_expression(expr_file)
meta <- read_metadata(meta_file)

if (!"sampnames" %in% colnames(meta)) stop("sampnames column missing in metadata")

meta_primary <- meta
if (identical(opt$primary_filter, "tcga_primary")) {
  meta_primary <- meta_primary %>% filter(grepl("_01", sampnames))
}
if (!nrow(meta_primary)) {
  stop("No metadata rows remain after applying primary_filter=", opt$primary_filter, ".")
}

# Ensure survival columns exist
stopifnot(all(c("OS", "OS_time", "PFI", "PFI_time") %in% colnames(meta_primary)))

covariates_used <- base_covariates[base_covariates %in% colnames(meta_primary)]

message("Running tertile-based survival analysis for overall BRCA cohort")
message("Total samples: ", nrow(meta_primary))

# Function to run and save results for one endpoint
run_endpoint <- function(meta_data, endpoint_event, endpoint_time, endpoint_label) {
  out_dir <- file.path(output_root, endpoint_label)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  lasso_ok <- requireNamespace("glmnet", quietly = TRUE)
  if (!lasso_ok) message("glmnet not available; skipping LASSO for ", endpoint_label)

  res <- analyze_survival_expression(
    expression_data = expr,
    metadata = meta_data,
    sample_col = "sampnames",
    gene_col = "Gene",
    genes = signature_genes,
    event_col = endpoint_event,
    time_col = endpoint_time,
    clinical_covariates = covariates_used,
    split_method = "tertiles",
    include_groups = c("low", "mid", "high"),
    steps_to_run = c("individual", "individual+covariates", "signature", "signature+covariates", "lasso", "km_individual"),
    run_lasso = lasso_ok,
    run_ph_checks = TRUE,
    run_time_dependent = TRUE,
    cut_months = 60,
    save_plots = TRUE,
    export_csv = TRUE,
    output_dir = out_dir,
    plot_params = list(plot_title = paste0(endpoint_label, " Survival Analysis (Tertiles)")),
    verbose = TRUE
  )

  saveRDS(res, file = file.path(out_dir, paste0(endpoint_label, "_results.rds")))

  # Plain text summary
  summary_lines <- c(
    paste0("Endpoint: ", endpoint_label),
    paste0("Samples analyzed: ", ifelse(!is.null(res$sample_groups), length(unique(na.omit(unlist(res$sample_groups)))), NA))
  )
  if (!is.null(res$results_univariate)) {
    summary_lines <- c(summary_lines, "\nUnivariate Cox (first 10 rows):", utils::capture.output(print(head(res$results_univariate, 10))))
  }
  if (!is.null(res$results_adjusted)) {
    summary_lines <- c(summary_lines, "\nAdjusted Cox (first 10 rows):", utils::capture.output(print(head(res$results_adjusted, 10))))
  }
  if (!is.null(res$results_signature)) {
    summary_lines <- c(summary_lines, "\nSignature Cox (unadjusted):", utils::capture.output(print(res$results_signature$unadjusted)))
    if (!is.null(res$results_signature$adjusted)) {
      summary_lines <- c(summary_lines, "\nSignature Cox (adjusted):", utils::capture.output(print(res$results_signature$adjusted)))
    }
  }
  if (!is.null(res$ph_diagnostics)) {
    summary_lines <- c(summary_lines, "\nProportional Hazards Test:", utils::capture.output(print(res$ph_diagnostics$ph_test)))
    if (!is.null(res$ph_diagnostics$time_dependent)) {
      summary_lines <- c(summary_lines, "\nTime-Dependent Model Summary:", utils::capture.output(print(res$ph_diagnostics$time_dependent)))
    }
  }
  writeLines(summary_lines, con = file.path(out_dir, paste0(endpoint_label, "_Summary.txt")))

  # Markdown report
  md <- c(paste0("# ", endpoint_label, " Survival Analysis (Tertiles)"), "")
  md <- c(
    md, "## Samples", "",
    paste0("Total samples (cohort): ", nrow(meta_data)),
    paste0("Analyzed samples (with survival & expression): ", ifelse(!is.null(res$sample_groups), length(unique(na.omit(unlist(res$sample_groups)))), NA)), ""
  )

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
    sep <- paste(rep("---", ncol(df_fmt)), collapse = " | ")
    rows <- apply(df_fmt, 1, function(r) paste(r, collapse = " | "))
    c(paste0("| ", header, " |"), paste0("| ", sep, " |"), paste0("| ", rows, " |"), "")
  }

  if (!is.null(res$results_univariate)) md <- c(md, "## Univariate Cox (top 10)", "", md_table(res$results_univariate, 10))
  if (!is.null(res$results_adjusted)) md <- c(md, "## Adjusted Cox (top 10)", "", md_table(res$results_adjusted, 10))

  if (!is.null(res$results_signature)) {
    md <- c(md, "## Signature Cox", "")
    if (!is.null(res$results_signature$unadjusted) && inherits(res$results_signature$unadjusted, "summary.coxph")) {
      su <- res$results_signature$unadjusted
      co <- as.data.frame(su$coefficients)
      ci <- as.data.frame(su$conf.int)
      tab <- cbind(Term = rownames(co), HR = ci$`exp(coef)`, CI_Lower = ci$`lower .95`, CI_Upper = ci$`upper .95`, P_Value = co$`Pr(>|z|)`)
      md <- c(md, "### Unadjusted", "", md_table(tab))
    }
    if (!is.null(res$results_signature$adjusted) && inherits(res$results_signature$adjusted, "summary.coxph")) {
      sa <- res$results_signature$adjusted
      co <- as.data.frame(sa$coefficients)
      ci <- as.data.frame(sa$conf.int)
      tab2 <- cbind(Term = rownames(co), HR = ci$`exp(coef)`, CI_Lower = ci$`lower .95`, CI_Upper = ci$`upper .95`, P_Value = co$`Pr(>|z|)`)
      md <- c(md, "### Adjusted", "", md_table(tab2))
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

  # Lightweight HTML conversion
  md_lines <- md
  to_html <- function(lines) {
    html <- c()
    i <- 1L
    n <- length(lines)
    esc <- function(x) {
      x <- gsub("&", "&amp;", x, fixed = TRUE)
      x <- gsub("<", "&lt;", x, fixed = TRUE)
      x <- gsub(">", "&gt;", x, fixed = TRUE)
      x
    }
    while (i <= n) {
      ln <- lines[i]
      if (grepl("^\\|", ln)) {
        block <- c()
        while (i <= n && grepl("^\\|", lines[i])) {
          block <- c(block, lines[i])
          i <- i + 1L
        }
        if (length(block) >= 2 && grepl("^\\| *-+", block[2])) {
          header_line <- block[1]
          body <- block[-c(1, 2)]
        } else {
          header_line <- block[1]
          body <- block[-1]
        }
        sf <- function(line) {
          line <- sub("^\\|", "", line)
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
        html <- c(html, paste0("<h3>", esc(sub("^### ", "", ln)), "</h3>"))
      } else if (grepl("^## ", ln)) {
        html <- c(html, paste0("<h2>", esc(sub("^## ", "", ln)), "</h2>"))
      } else if (grepl("^# ", ln)) {
        html <- c(html, paste0("<h1>", esc(sub("^# ", "", ln)), "</h1>"))
      } else if (nzchar(ln)) html <- c(html, paste0("<p>", esc(ln), "</p>"))
      i <- i + 1L
    }
    html
  }
  html_body <- to_html(md_lines)
  esc_html <- function(x) {
    x <- gsub("&", "&amp;", x, fixed = TRUE)
    x <- gsub("<", "&lt;", x, fixed = TRUE)
    x <- gsub(">", "&gt;", x, fixed = TRUE)
    x
  }
  html_doc <- c(
    "<!DOCTYPE html>", '<html lang="en"><head><meta charset="UTF-8"/>',
    sprintf("<title>%s Survival Analysis Summary (Tertiles)</title>", endpoint_label),
    "<style>body{font-family:Arial,Helvetica,sans-serif;margin:1.5rem;line-height:1.35;font-size:14px;} h1,h2,h3{font-weight:600;margin-top:1.8em;} table.report-table{border-collapse:collapse;margin:1em 0;font-size:13px;} table.report-table th,table.report-table td{border:1px solid #999;padding:4px 6px;text-align:left;} table.report-table th{background:#f2f2f2;} p{margin:0.4em 0;} .note{color:#555;font-size:12px;margin-top:3em;border-top:1px solid #ddd;padding-top:0.5em;} </style>",
    "</head><body>", html_body,
    sprintf('<div class="note">Generated on %s</div>', format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "</body></html>"
  )
  writeLines(html_doc, file.path(out_dir, paste0(endpoint_label, "_Summary.html")))

  invisible(res)
}

write_run_metadata(
  output_dir = output_root,
  expr_path = expr_file,
  meta_path = meta_file,
  n_before = nrow(meta),
  n_after = nrow(meta_primary),
  covariates_used = covariates_used
)

# Run both endpoints
message("\n==============================")
message("Processing Overall Survival (OS)")
message("==============================")
run_endpoint(meta_primary, endpoint_event = "OS", endpoint_time = "OS_time", endpoint_label = "OS")

message("\n==============================")
message("Processing Progression-Free Interval (PFI)")
message("==============================")
run_endpoint(meta_primary, endpoint_event = "PFI", endpoint_time = "PFI_time", endpoint_label = "PFI")

message("\nTertile-based overall BRCA survival analysis complete. See output in ", output_root, ".")
