# Resolve user-supplied human gene symbols without modifying the expression
# matrix.  Matrix rows remain the source-of-truth feature identifiers; this
# helper only maps requested feature names to rows that already exist.

normalize_gene_resolution_mode <- function(value) {
  cli_parse_enum_flag(
    value = value,
    allowed_values = c("none", "l2p"),
    flag_name = "gene_resolution",
    normalize = "lower"
  )
}

resolve_requested_gene_identifiers <- function(requested_genes, expression_gene_ids,
                                               mode = "l2p") {
  requested_genes <- as.character(requested_genes)
  if (!length(requested_genes)) {
    return(data.frame(
      requested_gene = character(0),
      l2psupp_oldname = character(0),
      l2psupp_newname = character(0),
      l2psupp_change_flag = integer(0),
      l2psupp_status = integer(0),
      l2psupp_is_legit_name = integer(0),
      resolved_gene = character(0),
      resolution = character(0),
      expression_row_found = logical(0),
      resolver = character(0),
      resolver_version = character(0),
      stringsAsFactors = FALSE
    ))
  }

  if (any(!nzchar(trimws(requested_genes)))) {
    stop("Gene identifiers must be non-empty after whitespace trimming.")
  }
  requested_genes <- trimws(requested_genes)
  if (anyDuplicated(expression_gene_ids)) {
    stop("Expression gene identifiers must be unique before gene-input resolution.")
  }

  mode <- normalize_gene_resolution_mode(mode)
  resolver_version <- NA_character_
  resolver <- "none"
  l2p_table <- data.frame(
    oldname = requested_genes,
    newname = requested_genes,
    change_flag = 0L,
    status = NA_integer_,
    is_legit_name = NA_integer_,
    stringsAsFactors = FALSE
  )

  if (mode == "l2p") {
    if (!requireNamespace("l2psupp", quietly = TRUE)) {
      stop(
        "--gene-resolution l2p requires the l2psupp package. ",
        "Use the published container image or install the pinned l2psupp dependency. ",
        "Use --gene-resolution none only when identifiers are already expression-row identifiers."
      )
    }
    resolver <- "l2psupp::updategenes"
    resolver_version <- as.character(utils::packageVersion("l2psupp"))
    updated <- l2psupp::updategenes(requested_genes)
    required_columns <- c("oldname", "newname", "change_flag", "status", "is_legit_name")
    if (!is.data.frame(updated) || nrow(updated) != length(requested_genes) ||
        !all(required_columns %in% names(updated))) {
      stop("Unexpected l2psupp::updategenes() result; refusing to resolve gene inputs silently.")
    }
    l2p_table <- updated[, required_columns, drop = FALSE]
    l2p_table[] <- lapply(l2p_table, as.vector)
  }

  # Prefer an exact user-supplied matrix identifier. This keeps a valid legacy
  # expression row usable while still recording L2P's suggested current name.
  # When no such row exists, use L2P's updated symbol only if it is present.
  original_found <- requested_genes %in% expression_gene_ids
  updated_names <- as.character(l2p_table$newname)
  updated_found <- updated_names %in% expression_gene_ids
  resolved_gene <- ifelse(original_found, requested_genes, updated_names)
  expression_row_found <- original_found | updated_found
  resolution <- ifelse(
    original_found & updated_names != requested_genes,
    "exact_input_preferred_over_l2p_update",
    ifelse(
      original_found,
      "exact_input",
      ifelse(updated_found & updated_names != requested_genes, "l2p_updated_match", "unmatched")
    )
  )

  data.frame(
    requested_gene = requested_genes,
    l2psupp_oldname = as.character(l2p_table$oldname),
    l2psupp_newname = updated_names,
    l2psupp_change_flag = as.integer(l2p_table$change_flag),
    l2psupp_status = as.integer(l2p_table$status),
    l2psupp_is_legit_name = as.integer(l2p_table$is_legit_name),
    resolved_gene = resolved_gene,
    resolution = resolution,
    expression_row_found = expression_row_found,
    resolver = resolver,
    resolver_version = resolver_version,
    stringsAsFactors = FALSE
  )
}

assert_gene_resolution_is_unambiguous <- function(resolution_table) {
  matched <- resolution_table[resolution_table$expression_row_found, , drop = FALSE]
  duplicate_targets <- unique(matched$resolved_gene[duplicated(matched$resolved_gene)])
  if (length(duplicate_targets)) {
    details <- vapply(duplicate_targets, function(target) {
      inputs <- matched$requested_gene[matched$resolved_gene == target]
      paste0(target, " <- ", paste(inputs, collapse = ", "))
    }, character(1))
    stop(
      "Multiple requested genes resolve to the same expression row; refusing to double-count a feature: ",
      paste(details, collapse = "; "), "."
    )
  }
  invisible(resolution_table)
}

replace_requested_gene_identifiers <- function(requested_genes, resolution_table) {
  positions <- match(requested_genes, resolution_table$requested_gene)
  if (anyNA(positions)) stop("Internal error: a requested gene is absent from its resolution table.")
  resolved <- resolution_table$resolved_gene[positions]
  missing <- resolution_table$requested_gene[positions][!resolution_table$expression_row_found[positions]]
  if (length(missing)) {
    stop(
      "Requested gene(s) do not match expression rows after gene resolution: ",
      paste(unique(missing), collapse = ", "), ". See gene_identifier_resolution.csv for L2P results."
    )
  }
  resolved
}

write_gene_resolution_report <- function(resolution_table, output_dir) {
  tables_dir <- file.path(output_dir, "tables")
  dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(resolution_table, file.path(tables_dir, "gene_identifier_resolution.csv"))
}
