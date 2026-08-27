parse_cli_args <- function(args) {
  output <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("unexpected argument: ", token, call. = FALSE)
    key <- gsub("-", "_", substring(token, 3L), fixed = TRUE)
    if (i == length(args) || startsWith(args[[i + 1L]], "--")) {
      output[[key]] <- TRUE
      i <- i + 1L
    } else {
      output[[key]] <- args[[i + 1L]]
      i <- i + 2L
    }
  }
  output
}

cli_value <- function(args, name, default = NULL, cast = identity) {
  value <- args[[name]]
  if (is.null(value)) value <- default
  if (is.null(value)) return(NULL)
  cast(value)
}

validate_cli_args <- function(args, allowed) {
  unknown <- setdiff(names(args), allowed)
  if (length(unknown)) {
    stop(
      "unknown option(s): ",
      paste0("--", gsub("_", "-", unknown), collapse = ", "),
      call. = FALSE
    )
  }
  invisible(args)
}
