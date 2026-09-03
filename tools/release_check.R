#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
strict <- "--strict" %in% args || identical(
  tolower(Sys.getenv("SCALABLECF_RELEASE_STRICT", unset = "false")), "true"
)
require_stata <- "--require-stata" %in% args || identical(
  tolower(Sys.getenv("SCALABLECF_REQUIRE_STATA", unset = "false")), "true"
)
root <- normalizePath(".", winslash = "/", mustWork = TRUE)
description_path <- file.path(root, "DESCRIPTION")
if (!file.exists(description_path)) {
  stop("Run tools/release_check.R from the package root", call. = FALSE)
}
description <- read.dcf(description_path)
version <- description[[1L, "Version"]]
failures <- character()
warnings <- character()
deferred <- character()

require_path <- function(path) {
  if (!file.exists(file.path(root, path))) {
    failures <<- c(failures, paste("missing required file:", path))
  }
}
for (path in c(
  ".Rbuildignore", "NEWS.md", "README.md",
  "inst/doc/API_STABILITY.md",
  file.path("inst", "doc", paste0("RELEASE_CHECKLIST_", version, ".md")),
  "inst/python/requirements-cuda.txt",
  "inst/python/requirements-cuda-windows-py312.lock"
)) require_path(path)

python_root <- file.path(root, "inst", "python")
python_files <- if (dir.exists(python_root)) {
  list.files(
    python_root, recursive = TRUE, all.files = TRUE,
    full.names = FALSE, include.dirs = FALSE
  )
} else {
  character()
}
forbidden <- python_files[grepl("(^|/)__pycache__/|[.]py[co]$", python_files)]
if (length(forbidden)) {
  message <- paste(
    "source tree contains generated Python cache files:",
    paste(file.path("inst", "python", forbidden), collapse = ", ")
  )
  if (strict) failures <- c(failures, message) else warnings <- c(warnings, message)
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
  message <- "archived Stata parity evidence is missing"
  if (require_stata) {
    failures <- c(failures, message)
  } else {
    deferred <- c(deferred, message)
  }
}
if (!grepl("^[0-9]+[.][0-9]+[.][0-9]+([.-][0-9A-Za-z.-]+)?$", version)) {
  failures <- c(failures, paste("invalid package Version:", version))
}
github_ref <- Sys.getenv("GITHUB_REF", unset = "")
if (strict && startsWith(github_ref, "refs/tags/")) {
  tag <- sub("^refs/tags/", "", github_ref)
  expected_tag <- paste0("v", version)
  if (!identical(tag, expected_tag)) {
    failures <- c(failures, paste(
      "release tag", tag, "does not match DESCRIPTION Version", version
    ))
  }
}

cat("scalableCounterfactual release preflight\n")
cat("  version:", description[[1L, "Version"]], "\n")
cat("  strict:", strict, "\n")
cat("  require Stata parity:", require_stata, "\n")
if (length(warnings)) {
  cat("Warnings:\n", paste0("  - ", warnings, collapse = "\n"), "\n")
}
if (length(deferred)) {
  cat("Deferred release evidence:\n",
      paste0("  - ", deferred, collapse = "\n"), "\n")
}
if (length(failures)) {
  cat("Unresolved release items:\n",
      paste0("  - ", failures, collapse = "\n"), "\n")
  quit(save = "no", status = 1L)
} else {
  cat("All release preflight checks passed.\n")
}
