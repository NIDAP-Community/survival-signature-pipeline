suppressPackageStartupMessages({
  library(readr)
})

parse_subgroup_comparison_field <- function(x) {
  pieces <- strsplit(as.character(x), "\\s*\\|\\s*")
  dataset <- vapply(pieces, function(part) if (length(part) >= 1) trimws(part[[1]]) else NA_character_, character(1))
  window <- vapply(pieces, function(part) if (length(part) >= 2) trimws(part[[2]]) else NA_character_, character(1))
  list(dataset = dataset, window = window)
}

normalize_subgroup_summary_table <- function(
    subgroup_df,
    subgroup_column = "subgroup",
    subgroup_type_column = "subgroup_type",
    feature_column = "feature",
    endpoint_column = "endpoint",
    comparison_column = "comparison",
    hr_column = "HR",
    ci_low_column = "CI_low",
    ci_high_column = "CI_high",
    p_value_column = "p_value",
    p_value_bh_column = "p_value_bh",
    dataset_column = NULL,
    window_column = NULL,
    allow_legacy_pooled_columns = FALSE) {
  shared_columns <- c(subgroup_column, feature_column, endpoint_column, hr_column, ci_low_column, ci_high_column, p_value_column, p_value_bh_column)
  missing_shared_columns <- setdiff(shared_columns, names(subgroup_df))

  if (!length(missing_shared_columns)) {
    normalized_df <- subgroup_df %>%
      dplyr::mutate(
        subgroup = .data[[subgroup_column]],
        subgroup_type = if (subgroup_type_column %in% names(subgroup_df)) .data[[subgroup_type_column]] else "subgroup",
        feature = .data[[feature_column]],
        endpoint = .data[[endpoint_column]],
        comparison = if (comparison_column %in% names(subgroup_df)) .data[[comparison_column]] else "Overall",
        HR = .data[[hr_column]],
        CI_low = .data[[ci_low_column]],
        CI_high = .data[[ci_high_column]],
        p_value = .data[[p_value_column]],
        p_value_bh = .data[[p_value_bh_column]]
      )

    comparison_parts <- parse_subgroup_comparison_field(normalized_df$comparison)

    normalized_df <- normalized_df %>%
      dplyr::mutate(
        dataset = if (!is.null(dataset_column) && dataset_column %in% names(subgroup_df)) .data[[dataset_column]] else if ("dataset" %in% names(subgroup_df)) .data[["dataset"]] else comparison_parts$dataset,
        window = if (!is.null(window_column) && window_column %in% names(subgroup_df)) .data[[window_column]] else if ("window" %in% names(subgroup_df)) .data[["window"]] else comparison_parts$window,
        subtype = if ("subtype" %in% names(subgroup_df)) .data[["subtype"]] else subgroup
      ) %>%
      dplyr::select(subgroup, subgroup_type, feature, endpoint, comparison, dataset, window, HR, CI_low, CI_high, p_value, p_value_bh, dplyr::everything())

    return(normalized_df)
  }

  legacy_columns <- c("subtype", "feature", "endpoint", "dataset", "window", "HR", "CI_low", "CI_high", "p_value", "p_value_bh")
  if (isTRUE(allow_legacy_pooled_columns) && all(legacy_columns %in% names(subgroup_df))) {
    normalized_df <- subgroup_df %>%
      dplyr::mutate(
        subgroup = .data$subtype,
        subgroup_type = "subtype",
        comparison = paste(.data$dataset, .data$window, sep = " | ")
      ) %>%
      dplyr::select(subgroup, subgroup_type, feature, endpoint, comparison, dataset, window, HR, CI_low, CI_high, p_value, p_value_bh, dplyr::everything())
    return(normalized_df)
  }

  stop(
    "Subgroup summary input is missing required columns for the shared contract: ",
    paste(shared_columns, collapse = ", "),
    if (isTRUE(allow_legacy_pooled_columns)) {
      paste0(
        ". Legacy pooled subtype columns are also accepted: ",
        paste(legacy_columns, collapse = ", "),
        "."
      )
    } else {
      "."
    }
  )
}

validate_subgroup_summary_table <- function(
    subgroup_df,
    context_label = "subgroup_summary_table",
    require_dataset_window = FALSE,
    require_subtype_alias = FALSE) {
  required_columns <- c("subgroup", "subgroup_type", "feature", "endpoint", "comparison", "HR", "CI_low", "CI_high", "p_value", "p_value_bh")
  if (isTRUE(require_dataset_window)) {
    required_columns <- c(required_columns, "dataset", "window")
  }
  if (isTRUE(require_subtype_alias)) {
    required_columns <- c(required_columns, "subtype")
  }

  missing_columns <- setdiff(required_columns, names(subgroup_df))
  if (length(missing_columns)) {
    stop(context_label, " is missing required columns after normalization: ", paste(missing_columns, collapse = ", "))
  }

  invisible(subgroup_df)
}

sanitize_path_component <- function(value) {
  sanitized <- tolower(trimws(as.character(value)))
  sanitized <- gsub("[^a-z0-9]+", "_", sanitized)
  sanitized <- gsub("^_+|_+$", "", sanitized)
  if (!nzchar(sanitized)) "run" else sanitized
}

default_run_output_dir <- function(topic, root = "results", timestamp = Sys.time()) {
  run_id <- format(timestamp, "%Y-%m-%d_run%H%M%S")
  file.path(root, sanitize_path_component(topic), run_id)
}

normalize_path_safe <- function(path) {
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

is_legacy_output_root <- function(path) {
  p <- gsub("\\\\", "/", as.character(path))
  grepl("(^|/)Output_", p)
}

get_git_sha <- function() {
  sha <- tryCatch(system2("git", c("rev-parse", "HEAD"), stdout = TRUE, stderr = FALSE), error = function(e) "")
  if (!length(sha)) return("NA")
  out <- trimws(sha[[1]])
  if (!nzchar(out)) "NA" else out
}

get_command_string <- function() {
  paste(commandArgs(), collapse = " ")
}

file_md5_or_na <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  if (dir.exists(path)) return(NA_character_)
  as.character(tools::md5sum(path)[[1]])
}

assert_output_guardrails <- function(output_root, allow_legacy_output = FALSE) {
  if (is_legacy_output_root(output_root) && !isTRUE(allow_legacy_output)) {
    stop(
      paste0(
        "Guardrail violation: output_root points to a legacy Output_* path (", output_root, "). ",
        "Use a run-scoped results path under results/<topic>/<YYYY-MM-DD_runNN>/. ",
        "If you intentionally need legacy output for backward compatibility, rerun with ",
        "SURVIVAL_ALLOW_LEGACY_OUTPUT=1."
      )
    )
  }
}

write_run_metadata_guardrail <- function(
    output_dir,
    purpose,
    driver,
    inputs = list(),
    parameters = list(),
    notes = character(0),
    scope = "exploratory") {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  logs_dir <- file.path(output_dir, "logs")
  if (!dir.exists(logs_dir)) dir.create(logs_dir, recursive = TRUE)

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
  md <- c(
    "# RUN_METADATA",
    "",
    paste0("Date: ", timestamp),
    paste0("Scope: ", scope),
    paste0("Purpose: ", purpose),
    paste0("Driver: ", driver),
    "",
    "## Execution",
    paste0("- Command: ", get_command_string()),
    paste0("- Git commit: ", get_git_sha()),
    paste0("- R version: ", R.version.string),
    ""
  )

  if (length(inputs)) {
    md <- c(md, "## Inputs")
    for (nm in names(inputs)) {
      p <- as.character(inputs[[nm]])
      path_kind <- if (dir.exists(p)) "directory" else "file"
      md <- c(md, paste0("- ", nm, ": ", normalize_path_safe(p), " (", path_kind, ", md5=", file_md5_or_na(p), ")"))
    }
    md <- c(md, "")
  }

  if (length(parameters)) {
    md <- c(md, "## Parameters")
    for (nm in names(parameters)) {
      md <- c(md, paste0("- ", nm, ": ", as.character(parameters[[nm]])))
    }
    md <- c(md, "")
  }

  if (length(notes)) {
    md <- c(md, "## Notes", paste0("- ", notes), "")
  }

  writeLines(md, con = file.path(output_dir, "RUN_METADATA.md"))
  writeLines(capture.output(sessionInfo()), con = file.path(logs_dir, "sessionInfo.txt"))
}
