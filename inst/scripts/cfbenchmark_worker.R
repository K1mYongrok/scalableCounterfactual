#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
task_index <- match("--task", args)
result_index <- match("--result", args)
if (is.na(task_index) || is.na(result_index) ||
    task_index == length(args) || result_index == length(args)) {
  stop("worker requires --task FILE and --result FILE", call. = FALSE)
}
suppressPackageStartupMessages(library(scalableCounterfactual))
scalableCounterfactual:::isolated_qr_benchmark_task(
  args[[task_index + 1L]], args[[result_index + 1L]]
)
