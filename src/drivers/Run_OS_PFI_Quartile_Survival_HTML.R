# Run_OS_PFI_Quartile_Survival_HTML.R
# Regenerate OS and PFI survival analyses using QUARTILES (low/mid/high)
# with Markdown + standalone HTML reports (mirrors tertile driver).

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(survival)
  library(survminer)
  library(ggplot2)
  if (!requireNamespace("cowplot", quietly = TRUE) && !requireNamespace("gridExtra", quietly = TRUE)) {
    message("NOTE: 'cowplot' or 'gridExtra' not available; statistics table may fallback.")
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
# Configuration (Quartiles)
# -----------------------------------------------------------------------------
expr_file <- file.path("data", "IODC-0001-TCGA-BRCA.normalized_counts.tsv")
meta_file <- file.path("data", "Survival_subtype_joined_metadata.csv")
output_root <- "outputs/quartile_analysis/overall"
allow_legacy_output <- identical(Sys.getenv("SURVIVAL_ALLOW_LEGACY_OUTPUT"), "1")
assert_output_guardrails(output_root, allow_legacy_output = allow_legacy_output)
if (!dir.exists(output_root)) dir.create(output_root, recursive = TRUE)

signature_genes <- c("SAT1","SRM","ODC1","AMD1")

covariates <- c("Diagnosis_Age","ajcc_pathologic_tumor_stage","race_collapsed","Subtype")

read_expression <- function(path) {
  expr <- readr::read_tsv(path, show_col_types = FALSE)
  if (!"Gene" %in% colnames(expr)) {
    colnames(expr)[1] <- "Gene"
  }
  expr
}

read_metadata <- function(path) {
  meta <- readr::read_csv(path, show_col_types = FALSE)
  colnames(meta) <- gsub("\n","", colnames(meta))
  meta
}

expr <- read_expression(expr_file)
meta <- read_metadata(meta_file)

if (!"sampnames" %in% colnames(meta)) stop("sampnames column missing in metadata")
meta_primary <- meta %>% dplyr::filter(grepl("_01", sampnames))

message("Samples (all): ", nrow(meta), "; primary tumor subset: ", nrow(meta_primary))
stopifnot(all(c("OS","OS_time","PFI","PFI_time") %in% colnames(meta_primary)))

covariates_used <- covariates[covariates %in% colnames(meta_primary)]

write_run_metadata_guardrail(
  output_dir = output_root,
  purpose = "Run OS and PFI survival analyses using a 4-gene signature with quartile grouping.",
  driver = "src/drivers/Run_OS_PFI_Quartile_Survival_HTML.R",
  inputs = list(Expression_TSV = expr_file, Metadata_CSV = meta_file),
  parameters = list(
    primary_filter = "tcga_primary",
    metadata_rows_pre_filter = nrow(meta),
    metadata_rows_post_filter = nrow(meta_primary),
    covariates_used = ifelse(length(covariates_used), paste(covariates_used, collapse = ", "), "<none>"),
    output_root = normalize_path_safe(output_root),
    allow_legacy_output = allow_legacy_output
  ),
  notes = if (is_legacy_output_root(output_root)) "Legacy Output_* root used via explicit override (SURVIVAL_ALLOW_LEGACY_OUTPUT=1)." else character(0)
)

# -----------------------------------------------------------------------------
# Wrapper (same as tertile script, but split_method = 'quartiles')
# -----------------------------------------------------------------------------
run_endpoint <- function(event_col, time_col, label) {
  out_dir <- file.path(output_root, label)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  lasso_ok <- requireNamespace("glmnet", quietly = TRUE)
  if (!lasso_ok) message("glmnet not available; skipping LASSO for ", label)

  res <- analyze_survival_expression(
    expression_data = expr,
    metadata = meta_primary,
    sample_col = "sampnames",
    gene_col = "Gene",
    genes = signature_genes,
    event_col = event_col,
    time_col = time_col,
    clinical_covariates = covariates_used,
    split_method = "quartiles",
    include_groups = c("low","mid","high"),
    steps_to_run = c("individual","individual+covariates","signature","signature+covariates","lasso","km_individual"),
    run_lasso = lasso_ok,
    run_ph_checks = TRUE,
    run_time_dependent = TRUE,
    cut_months = 60,
    save_plots = TRUE,
    export_csv = TRUE,
    output_dir = out_dir,
    plot_params = list(plot_title = paste0(label, " Survival Analysis")),
    verbose = TRUE
  )

  if (!is.null(res$plots$signature)) {
    message("Signature KM statistics panel mode: bottom_table (if supporting packages present)")
  }

  saveRDS(res, file = file.path(out_dir, paste0(label, "_results.rds")))

  # Plain text summary (legacy)
  summary_lines <- c()
  summary_lines <- c(summary_lines, paste0("Endpoint: ", label))
  summary_lines <- c(summary_lines, paste0("Samples analyzed: ", ifelse(!is.null(res$sample_groups), length(unique(na.omit(unlist(res$sample_groups)))), NA)))
  if (!is.null(res$results_univariate)) {
    summary_lines <- c(summary_lines, "\nUnivariate Cox (first 10 rows):")
    summary_lines <- c(summary_lines, utils::capture.output(print(head(res$results_univariate, 10))))
  }
  if (!is.null(res$results_adjusted)) {
    summary_lines <- c(summary_lines, "\nAdjusted Cox (first 10 rows):")
    summary_lines <- c(summary_lines, utils::capture.output(print(head(res$results_adjusted, 10))))
  }
  if (!is.null(res$results_signature)) {
    summary_lines <- c(summary_lines, "\nSignature Cox (unadjusted):")
    summary_lines <- c(summary_lines, utils::capture.output(print(res$results_signature$unadjusted)))
    if (!is.null(res$results_signature$adjusted)) {
      summary_lines <- c(summary_lines, "\nSignature Cox (adjusted):")
      summary_lines <- c(summary_lines, utils::capture.output(print(res$results_signature$adjusted)))
    }
  }
  if (!is.null(res$lasso_selected)) {
    summary_lines <- c(summary_lines, paste0("LASSO selected genes: ", paste(res$lasso_selected, collapse = ", ")))
  }
  if (!is.null(res$ph_diagnostics)) {
    summary_lines <- c(summary_lines, "\nProportional Hazards Test:")
    summary_lines <- c(summary_lines, utils::capture.output(print(res$ph_diagnostics$ph_test)))
    if (!is.null(res$ph_diagnostics$time_dependent)) {
      summary_lines <- c(summary_lines, "\nTime-Dependent Model Summary:")
      summary_lines <- c(summary_lines, utils::capture.output(print(res$ph_diagnostics$time_dependent)))
    }
  }
  writeLines(summary_lines, con = file.path(out_dir, paste0(label, "_Summary.txt")))

  # ---------------- Markdown & HTML Report -----------------
  md <- c(paste0("# ", label, " Survival Analysis Summary"), "")
  md <- c(md, "## Samples", "",
          paste0("Total samples (metadata): ", nrow(meta_primary)),
          paste0("Analyzed samples (with survival & expression): ", ifelse(!is.null(res$sample_groups), length(unique(na.omit(unlist(res$sample_groups)))), NA)), "")

  md_table <- function(df, max_rows = NULL, digits = 4) {
    if (is.null(df) || length(df) == 0) return(character(0))
    if (!is.data.frame(df)) df <- as.data.frame(df)
    if (nrow(df) == 0) return(character(0))
    if (!is.null(max_rows)) df <- head(df, max_rows)
    df_fmt <- df
    num_cols <- vapply(df_fmt, is.numeric, logical(1))
    df_fmt[num_cols] <- lapply(df_fmt[num_cols], function(x) format(signif(x, digits), trim = TRUE, scientific = FALSE))
    header <- paste(colnames(df_fmt), collapse = " | ")
    sep <- paste(vapply(colnames(df_fmt), function(x) "---", character(1)), collapse = " | ")
    rows <- apply(df_fmt, 1, function(r) paste(r, collapse = " | "))
    c(paste0("| ", header, " |"), paste0("| ", sep, " |"), paste0("| ", rows, " |"), "")
  }

  if (!is.null(res$results_univariate)) md <- c(md, "## Univariate Cox (top 10)", "", md_table(res$results_univariate, 10))
  if (!is.null(res$results_adjusted))  md <- c(md, "## Adjusted Cox (top 10)", "", md_table(res$results_adjusted, 10))

  if (!is.null(res$results_signature)) {
    md <- c(md, "## Signature Cox", "")
    if (!is.null(res$results_signature$unadjusted) && inherits(res$results_signature$unadjusted, "summary.coxph")) {
      su <- res$results_signature$unadjusted
      co <- as.data.frame(su$coefficients); ci <- as.data.frame(su$conf.int)
      tab <- cbind(Term = rownames(co), HR = ci$`exp(coef)`, CI_Lower = ci$`lower .95`, CI_Upper = ci$`upper .95`, P_Value = co$`Pr(>|z|)`)
      md <- c(md, "### Unadjusted", "", md_table(tab))
    }
    if (!is.null(res$results_signature$adjusted) && inherits(res$results_signature$adjusted, "summary.coxph")) {
      sa <- res$results_signature$adjusted
      co <- as.data.frame(sa$coefficients); ci <- as.data.frame(sa$conf.int)
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
      co <- as.data.frame(td$coefficients); ci <- as.data.frame(td$conf.int)
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

  md_path <- file.path(out_dir, paste0(label, "_Summary.md"))
  writeLines(md, md_path)

  # HTML conversion (reuse lightweight converter)
  md_lines <- md
  to_html <- function(lines) {
    html <- c(); i <- 1L; n <- length(lines)
    esc <- function(x){x <- gsub('&','&amp;',x,fixed=TRUE);x <- gsub('<','&lt;',x,fixed=TRUE);x <- gsub('>','&gt;',x,fixed=TRUE);x}
    while (i <= n) {
      ln <- lines[i]
      if (grepl('^\\|', ln)) {
        block <- c(); while(i<=n && grepl('^\\|', lines[i])){block <- c(block, lines[i]); i <- i+1L}
        if (length(block) >=2 && grepl('^\\| *-+', block[2])) {header_line <- block[1]; body <- block[-c(1,2)]} else {header_line <- block[1]; body <- block[-1]}
        sf <- function(line){ line <- sub('^\\|','',line); line <- sub('\\| *$','',line); trimws(strsplit(line,'\\|')[[1]]) }
        heads <- sf(header_line)
        html <- c(html,'<table class="report-table">','<thead><tr>', paste0('<th>',esc(heads),'</th>', collapse=''), '</tr></thead>','<tbody>')
        if (length(body)) for (bl in body){ if(!nzchar(gsub('[\n[:space:]]','',bl))) next; cells <- sf(bl); if(length(cells)==length(heads)) html <- c(html,'<tr>',paste0('<td>',esc(cells),'</td>',collapse=''),'</tr>') }
        html <- c(html,'</tbody></table>'); next
      }
      if (grepl('^### ', ln)) html <- c(html, paste0('<h3>',esc(sub('^### ','',ln)),'</h3>'))
      else if (grepl('^## ', ln)) html <- c(html, paste0('<h2>',esc(sub('^## ','',ln)),'</h2>'))
      else if (grepl('^# ', ln)) html <- c(html, paste0('<h1>',esc(sub('^# ','',ln)),'</h1>'))
      else if (nzchar(ln)) html <- c(html, paste0('<p>',esc(ln),'</p>'))
      i <- i + 1L
    }
    html
  }
  html_body <- to_html(md_lines)
  html_doc <- c(
    '<!DOCTYPE html>','<html lang="en"><head><meta charset="UTF-8"/>',
    sprintf('<title>%s Survival Analysis Summary (Quartiles)</title>', label),
    '<style>body{font-family:Arial,Helvetica,sans-serif;margin:1.5rem;line-height:1.35;font-size:14px;} h1,h2,h3{font-weight:600;margin-top:1.8em;} table.report-table{border-collapse:collapse;margin:1em 0;font-size:13px;} table.report-table th,table.report-table td{border:1px solid #999;padding:4px 6px;text-align:left;} table.report-table th{background:#f2f2f2;} p{margin:0.4em 0;} .note{color:#555;font-size:12px;margin-top:3em;border-top:1px solid #ddd;padding-top:0.5em;} </style>',
    '</head><body>', html_body,
    sprintf('<div class="note">Generated on %s</div>', format(Sys.time(), '%Y-%m-%d %H:%M:%S %Z')),
    '</body></html>'
  )
  writeLines(html_doc, file.path(out_dir, paste0(label, '_Summary.html')))

  invisible(res)
}

# Run endpoints (quartiles)
res_os <- run_endpoint("OS", "OS_time", "OS")
res_pfi <- run_endpoint("PFI", "PFI_time", "PFI")

message("Quartile HTML reports complete. See output in ", output_root)
