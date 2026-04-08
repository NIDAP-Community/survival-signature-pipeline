suppressPackageStartupMessages({
  library(readr)
})

cli_flag_name <- function(name) {
  paste0("--", gsub("_", "-", name))
}

cli_help_requested <- function(args) {
  any(args %in% c("--help", "-h", "help"))
}

cli_default_label <- function(value) {
  value <- as.character(value)
  if (!length(value) || !nzchar(value)) "<empty>" else value
}

print_standard_cli_help <- function(
    title,
    workflow_name,
    script_path,
    defaults,
    summary_lines = character(0),
    required_flags = character(0),
    aliases = character(0)) {
  lines <- c(
    title,
    "",
    "Usage:",
    paste0("  Rscript entrypoint.R ", workflow_name, " [--flag value]..."),
    paste0("  Rscript ", script_path, " [--flag value]..."),
    ""
  )

  if (length(summary_lines)) {
    lines <- c(lines, "Key notes:", paste0("  - ", summary_lines), "")
  }

  if (length(required_flags)) {
    lines <- c(lines, "Required flags:", paste0("  ", cli_flag_name(required_flags)), "")
  }

  optional_flags <- setdiff(names(defaults), required_flags)
  if (length(optional_flags)) {
    lines <- c(
      lines,
      "Optional flags:",
      vapply(optional_flags, function(name) {
        paste0("  ", cli_flag_name(name), " (default: ", cli_default_label(defaults[[name]]), ")")
      }, character(1)),
      ""
    )
  }

  if (length(aliases)) {
    lines <- c(
      lines,
      "Compatibility aliases:",
      vapply(names(aliases), function(alias_name) {
        paste0("  ", cli_flag_name(alias_name), " -> ", cli_flag_name(aliases[[alias_name]]))
      }, character(1)),
      ""
    )
  }

  cat(paste(lines, collapse = "\n"), "\n")
}

stop_cli_invalid_token <- function(token, workflow_name) {
  stop(
    "Unexpected token for ", workflow_name, ": ", token, ". ",
    "Arguments must be supplied as --flag value pairs. ",
    "Run `Rscript entrypoint.R ", workflow_name, " --help` for usage."
  )
}

stop_cli_unknown_argument <- function(argument, workflow_name) {
  stop(
    "Unknown argument for ", workflow_name, ": ", argument, ". ",
    "Run `Rscript entrypoint.R ", workflow_name, " --help` for usage."
  )
}

stop_cli_missing_value <- function(argument, workflow_name) {
  stop(
    "Missing value for ", argument, " in ", workflow_name, ". ",
    "Arguments must be supplied as --flag value pairs."
  )
}

stop_cli_missing_required <- function(required_flags, workflow_name) {
  stop(
    "Missing required argument(s) for ", workflow_name, ": ",
    paste(cli_flag_name(required_flags), collapse = ", "), ". ",
    "Run `Rscript entrypoint.R ", workflow_name, " --help` for usage."
  )
}

cli_parse_csv_list <- function(value, uppercase = FALSE, lowercase = FALSE, allow_empty = FALSE) {
  items <- trimws(unlist(strsplit(as.character(value), ",", fixed = TRUE)))
  items <- items[nzchar(items)]

  if (!length(items) && !isTRUE(allow_empty)) {
    return(character(0))
  }
  if (isTRUE(uppercase) && isTRUE(lowercase)) {
    stop("cli_parse_csv_list cannot apply both uppercase and lowercase normalization.")
  }
  if (isTRUE(uppercase)) {
    items <- toupper(items)
  }
  if (isTRUE(lowercase)) {
    items <- tolower(items)
  }
  items
}

cli_parse_numeric_flag <- function(value, flag_name) {
  numeric_value <- suppressWarnings(as.numeric(value))
  if (!is.finite(numeric_value)) {
    stop(cli_flag_name(flag_name), " must be numeric.")
  }
  numeric_value
}

cli_parse_enum_flag <- function(value, allowed_values, flag_name, normalize = c("none", "lower", "upper")) {
  normalize <- match.arg(normalize)
  parsed_value <- as.character(value)

  if (normalize == "lower") {
    parsed_value <- tolower(parsed_value)
    allowed_values <- tolower(allowed_values)
  } else if (normalize == "upper") {
    parsed_value <- toupper(parsed_value)
    allowed_values <- toupper(allowed_values)
  }

  if (!parsed_value %in% allowed_values) {
    stop(cli_flag_name(flag_name), " must be one of: ", paste(allowed_values, collapse = ", "))
  }

  parsed_value
}

cli_assert_existing_input_files <- function(input_paths, workflow_name) {
  missing_paths <- input_paths[!file.exists(input_paths)]
  if (!length(missing_paths)) {
    return(invisible(input_paths))
  }

  details <- paste0(cli_flag_name(names(missing_paths)), "=", unname(missing_paths))
  stop(
    "Missing required input file(s) for ", workflow_name, ": ",
    paste(details, collapse = ", "), "."
  )
}

read_expression_table <- function(path, gene_column) {
  expression_data <- readr::read_tsv(path, show_col_types = FALSE)
  if (!gene_column %in% colnames(expression_data)) {
    colnames(expression_data)[1] <- gene_column
  }
  expression_data
}

read_metadata_table <- function(path) {
  metadata <- readr::read_csv(path, show_col_types = FALSE)
  colnames(metadata) <- gsub("\\n", "", colnames(metadata))
  metadata
}

resolve_standard_endpoint_columns <- function(endpoint_name, options, workflow_name) {
  endpoint_name <- toupper(endpoint_name)
  if (endpoint_name == "OS") {
    return(list(event = options$os_event_column, time = options$os_time_column))
  }
  if (endpoint_name == "PFI") {
    return(list(event = options$pfi_event_column, time = options$pfi_time_column))
  }

  stop(
    "Unsupported endpoint for ", workflow_name, ": ", endpoint_name, ". ",
    "Supported endpoints in the maintained interface are OS and PFI."
  )
}

parse_standard_cli_args <- function(
    defaults,
    workflow_name,
    help_printer,
    required_flags = character(0),
    aliases = character(0),
    args = commandArgs(trailingOnly = TRUE)) {
  if (!length(args) || cli_help_requested(args)) {
    help_printer()
    quit(save = "no", status = 0)
  }

  index <- 1
  while (index <= length(args)) {
    key <- args[[index]]
    if (!startsWith(key, "--")) {
      stop_cli_invalid_token(key, workflow_name)
    }

    key_clean <- gsub("-", "_", substring(key, 3))
    if (key_clean %in% names(aliases)) {
      key_clean <- aliases[[key_clean]]
    }

    if (!key_clean %in% names(defaults)) {
      stop_cli_unknown_argument(key, workflow_name)
    }
    if (index == length(args)) {
      stop_cli_missing_value(key, workflow_name)
    }

    next_value <- args[[index + 1]]
    if (startsWith(next_value, "--")) {
      stop_cli_missing_value(key, workflow_name)
    }

    defaults[[key_clean]] <- next_value
    index <- index + 2
  }

  if (length(required_flags)) {
    missing_required <- required_flags[vapply(required_flags, function(name) {
      !nzchar(as.character(defaults[[name]]))
    }, logical(1))]

    if (length(missing_required)) {
      stop_cli_missing_required(missing_required, workflow_name)
    }
  }

  defaults
}