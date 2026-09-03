# Solver benchmark protocol

`benchmark_qr_solvers()` is intended for computational comparisons, not for
selecting a solver after inspecting preferred empirical results.

For a fair comparison it freezes:

- the analysis rows and group split;
- the outcome and common design matrices;
- normalized sampling weights;
- conditional and reported quantile grids;
- the point-estimation worker count.

The requested point-worker count is recorded separately from the effective
count, which is capped at two group fits. Distribution-regression backend
comparisons also report each group's conditional-grid size and effective
threshold-worker count; the combined field is missing when the two differ.

Each solver receives a stable solver-specific random seed, so PFNB/PFN
preprocessing is reproducible and does not depend on execution order.
Optional stratified subsampling preserves the requested total size while
guaranteeing that each group has more observations than retained design
columns. For highly imbalanced groups, allocation is clamped to the feasible
range instead of requesting more rows than the smaller group contains.

Repeated benchmarks randomize execution order and report both raw timings and
median-based summaries. Warmup repetitions are supported and discarded.

CSV loading and formula parsing are excluded from solver timing. For each
solver, the output records success or the original error, elapsed point time,
wrapper time, peak R heap reported by `gc()`, convergence flags, and maximum
absolute differences from the chosen reference in effects and coefficients.

The in-session functions report R heap. `benchmark_qr_scaling()` implements the
publication path: every solver/repetition runs in a fresh R session and the
parent process RSS is sampled at a configurable interval. Use
`point_workers = 1` when comparing RSS because child-worker memory is not
included in the parent RSS field. Reports should still state CPU, RAM, R and
package versions, poll interval, sample sizes, quantile grid, and worker count.

Scaling checkpoints include the R/platform, BLAS/LAPACK, package, extension,
and applicable GPU runtime identity. Their row keys and stored fit shapes are
validated before reuse. Resumption uses absolute warmup-plus-repetition indices
so the recorded execution order matches an uninterrupted run. An all-missing
RSS condition remains `NA`.

Solvers may differ in numerical robustness. In particular, a solver failure is
reported as a result and is never replaced silently by another backend.
Warnings, diagnostic availability, preconditioning method, and condition
estimates are also recorded. For BR, FN, and PFN, the low-level quantreg return
objects do not expose convergence flags; the benchmark reports this as missing
rather than as zero failures.

CLI raw and summary tables are staged and committed together. Input, raw,
summary, and checkpoint paths must be distinct, which prevents a failed or
misconfigured benchmark from overwriting its source data.
