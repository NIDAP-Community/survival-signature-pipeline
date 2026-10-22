#!/usr/bin/env Rscript

workflow_specs <- list(
	"piecewise-dataset-interaction" = list(
		script = file.path("drivers", "run_piecewise_dataset_interaction.R"),
		description = "Run pooled piecewise Cox models with dataset interaction across two harmonized datasets."
	),
	"piecewise-subtype-specific" = list(
		script = file.path("drivers", "run_piecewise_subtype_specific.R"),
		description = "Run subtype-specific pooled piecewise Cox models across two harmonized datasets."
	),
	"single-cohort-survival" = list(
		script = file.path("drivers", "run_single_cohort_survival.R"),
		description = "Run single-cohort survival analysis through the reusable core engine."
	),
	"single-cohort-reporting-bundle" = list(
		script = file.path("drivers", "run_single_cohort_reporting_bundle.R"),
		description = "Assemble single-cohort reporting artifacts from a prior single-cohort run directory."
	),
	"subgroup-forest-bundle" = list(
		script = file.path("drivers", "run_subgroup_forest_bundle.R"),
		description = "Build generic subgroup forest plots from a shared subgroup summary table."
	),
	"signature-consistency-summary" = list(
		script = file.path("drivers", "run_signature_consistency_summary.R"),
		description = "Compare overall versus subtype-specific signature results from upstream analysis tables."
	),
	"finalized-pooled-results-bundle" = list(
		script = file.path("drivers", "run_finalized_pooled_results_bundle.R"),
		description = "Assemble pooled forest plots, subtype heatmaps, and summary text from upstream result tables."
	)
)

print_usage <- function() {
	cat(
		paste(
			"Reusable Survival Analysis Pipeline",
			"",
			"Usage:",
			"  Rscript entrypoint.R list",
			"  Rscript entrypoint.R <workflow> [workflow arguments]",
			"  Rscript entrypoint.R <workflow> --help",
			"",
			"Available workflows:",
			paste(vapply(names(workflow_specs), function(name) {
				paste0("  - ", name, ": ", workflow_specs[[name]]$description)
			}, character(1)), collapse = "\n"),
			"",
			"Examples:",
			"  Rscript entrypoint.R piecewise-dataset-interaction --output-dir results/example_run",
			"  Rscript entrypoint.R piecewise-subtype-specific --subtype-list BRCA_Basal --adjust-covariates none",
			"  Rscript entrypoint.R single-cohort-survival --expr-file test_data/tcga_expression.tsv --meta-file test_data/tcga_metadata.csv --endpoints OS",
			"  Rscript entrypoint.R single-cohort-reporting-bundle --input-run-dir results/smoke_test_single_cohort --subgroup-summary-table test_data/single_cohort_subgroup_summary.csv",
			"  Rscript entrypoint.R subgroup-forest-bundle --subgroup-summary-table test_data/single_cohort_subgroup_summary.csv --subgroup-type-filter subtype",
			"  Rscript entrypoint.R signature-consistency-summary --overall-implied-csv <path> --overall-summary-csv <path> --subtype-implied-csv <path> --subtype-summary-csv <path> --subtype-terms-csv <path>",
			"  Rscript entrypoint.R finalized-pooled-results-bundle --overall-signature-table <path> --pooled-term-table <path> --subgroup-summary-table <path>",
			sep = "\n"
		),
		"\n"
	)
}

args <- commandArgs(trailingOnly = TRUE)

if (!length(args) || args[[1]] %in% c("--help", "-h", "help", "list", "--list")) {
	print_usage()
	quit(save = "no", status = 0)
}

workflow_name <- args[[1]]
forwarded_args <- args[-1]

if (!workflow_name %in% names(workflow_specs)) {
	cat("Unknown workflow: ", workflow_name, "\n\n", sep = "")
	print_usage()
	quit(save = "no", status = 1)
}

workflow <- workflow_specs[[workflow_name]]
script_path <- workflow$script

if (!file.exists(script_path)) {
	stop("Workflow script not found: ", script_path)
}

rscript_bin <- file.path(R.home("bin"), "Rscript")
# system2() ultimately constructs a command line on Unix-like platforms. Quote
# each independently parsed argument so values such as comma-separated display
# labels with spaces arrive at the selected driver as one commandArgs() token.
child_args <- shQuote(c(script_path, forwarded_args))
status <- system2(rscript_bin, args = child_args)
quit(save = "no", status = as.integer(status))
