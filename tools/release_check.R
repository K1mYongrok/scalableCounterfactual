#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
strict <- "--strict" %in% args || identical(
  tolower(Sys.getenv("SCALABLECF_RELEASE_STRICT", unset = "false")), "true"
)
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
description_path <- file.path(root, "DESCRIPTION")
if (!file.exists(description_path)) {
  stop("Run tools/release_check.R from the package root", call. = FALSE)
}
description <- read.dcf(description_path)
failures <- character()
warnings <- character()

require_path <- function(path) {
  if (!file.exists(file.path(root, path))) {
    failures <<- c(failures, paste("missing required file:", path))
  }
}
for (path in c(
  ".Rbuildignore", "NEWS.md", "README.md",
  "inst/doc/API_STABILITY.md", "inst/doc/RELEASE_CHECKLIST_1.0.0.md",
  "inst/python/requirements-cuda.txt",
  "inst/python/requirements-cuda-windows-py312.lock"
)) require_path(path)

source_files <- list.files(root, recursive = TRUE, all.files = TRUE,
                           full.names = FALSE, include.dirs = FALSE)
forbidden <- source_files[grepl("(__pycache__|[.]py[co]$)", source_files)]
if (length(forbidden)) {
  warnings <- c(warnings, paste(
    "source tree contains generated Python cache files:",
    paste(forbidden, collapse = ", ")
  ))
}

authors <- description[[1L, "Authors@R"]]
if (grepl("Project.*Team|noreply@example[.]com", authors)) {
  failures <- c(failures, "Authors@R still contains the placeholder maintainer")
}
for (field in c("URL", "BugReports")) {
  if (!field %in% colnames(description) || !nzchar(description[[1L, field]])) {
    warnings <- c(warnings, paste(
      "DESCRIPTION is missing", field,
      "because no public repository has been created"
    ))
  }
}
if (!file.exists(file.path(root, "inst", "CITATION"))) {
  failures <- c(failures, "inst/CITATION is missing")
}
if (!file.exists(file.path(
  root, "inst", "provenance", "stata_parity_verified.csv"
))) {
  failures <- c(failures, "archived Stata parity evidence is missing")
}
if (strict && description[[1L, "Version"]] != "1.0.0") {
  failures <- c(failures, "strict release check requires Version 1.0.0")
}

cat("scalableCounterfactual release preflight\n")
cat("  version:", description[[1L, "Version"]], "\n")
cat("  strict:", strict, "\n")
if (length(warnings)) {
  cat("Warnings:\n", paste0("  - ", warnings, collapse = "\n"), "\n")
}
if (length(failures)) {
  cat("Unresolved release items:\n",
      paste0("  - ", failures, collapse = "\n"), "\n")
  if (strict) quit(save = "no", status = 1L)
} else {
  cat("All release preflight checks passed.\n")
}
