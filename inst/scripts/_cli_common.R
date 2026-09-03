as_bool <- function(x) {
  normalized <- tolower(trimws(as.character(x)))
  if (length(normalized) != 1L || is.na(normalized) ||
      !normalized %in% c("1", "0", "true", "false", "yes", "no", "y", "n")) {
    stop("expected a true/false value, got: ", x, call. = FALSE)
  }
  normalized %in% c("1", "true", "yes", "y")
}

canonical_cli_path <- function(path) {
  normalized <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (.Platform$OS.type == "windows") tolower(normalized) else normalized
}

assert_distinct_cli_paths <- function(paths) {
  paths <- paths[!vapply(paths, is.null, logical(1L))]
  if (!length(paths)) return(invisible(NULL))
  if (is.null(names(paths)) || any(!nzchar(names(paths)))) {
    stop("paths must be a named list", call. = FALSE)
  }
  canonical <- vapply(paths, canonical_cli_path, character(1L))
  duplicate <- duplicated(canonical) | duplicated(canonical, fromLast = TRUE)
  if (any(duplicate)) {
    conflicts <- split(names(canonical)[duplicate], canonical[duplicate])
    labels <- vapply(conflicts, paste, collapse = " = ", character(1L))
    stop(
      "path collision: ", paste(labels, collapse = "; "),
      call. = FALSE
    )
  }
  invisible(NULL)
}

assert_cli_directories_not_below_files <- function(directories, files) {
  directories <- directories[!vapply(directories, is.null, logical(1L))]
  files <- files[!vapply(files, is.null, logical(1L))]
  if (!length(directories) || !length(files)) return(invisible(NULL))
  directory_paths <- vapply(directories, canonical_cli_path, character(1L))
  file_paths <- vapply(files, canonical_cli_path, character(1L))
  for (directory_name in names(directory_paths)) {
    directory <- directory_paths[[directory_name]]
    conflicts <- names(file_paths)[vapply(file_paths, function(file) {
      identical(directory, file) || startsWith(directory, paste0(file, "/"))
    }, logical(1L))]
    if (length(conflicts)) {
      stop(
        "path collision: ", directory_name,
        " is at or below file path ", conflicts[[1L]],
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}
