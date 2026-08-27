# Release validation for 0.15.0

Validation was run on 2026-08-10 with R 4.5.2 on 64-bit Windows. The
machine-level evidence is stored outside the package source under
`output/package_validation/scalable_counterfactual_0.15.0_release`.

The share-ready source tree was rebuilt and checked again on 2026-08-11 after
adding the documentation index and installed examples. `R CMD check
--no-manual` completed with status `OK`, and the self-contained PFNB quick-start
example completed successfully with decomposition, bootstrap, plot, and output
files.

## Package and device checks

- A clean source tarball passed `R CMD check --no-manual` with status `OK`.
- The API contract test passed with output schema version `1.0`.
- The real-device CUDA suite passed using Python 3.12.1, NumPy 2.5.1,
  CuPy 14.1.1 for CUDA 12, and an NVIDIA GeForce RTX 4060 Ti.
- The source tarball excludes Python bytecode, cache directories, CI files,
  development tools, and generated output.

The multi-operating-system CPU workflow and the self-hosted CUDA workflow are
present, but were not executed by a remote CI service for this candidate.

## Analytic simulation

The release suite used three data-generating processes, 100 replications per
process, 1,000 observations per group, weighted estimation, and 49 regression
quantiles. All 300 fits completed. Across the reported structure, composition,
and total effects at Q10, Q50, and Q90, the maximum absolute bias was 0.0163
log point.

## Solver benchmark

The generic benchmark used 5,000 observations, 100 regression quantiles, one
warm-up, and three timed repetitions. BR is the exact reference.

| Solver | Median seconds | Speed relative to BR | Maximum effect difference from BR |
|---|---:|---:|---:|
| BR | 3.40 | 1.00x | 0 |
| FN | 1.51 | 2.25x | 7.61e-09 |
| QFNB | 1.38 | 2.46x | 7.61e-09 |
| PFNB | 1.16 | 2.93x | 2.77e-06 |
| Onestep | 1.25 | 2.72x | 1.92e-02 |

BR, FN, QFNB, and PFNB are treated as exact solvers within their numerical
tolerances. Onestep remains explicitly labelled approximate. BR emitted a
non-unique-solution warning in each timed run, but no run failed.

## Bootstrap coverage smoke checks

Each inference smoke check used three data-generating processes, 30 outer
replications per process, 500 observations per group, and 49 bootstrap draws.
All 90 tasks completed for both bootstrap schemes and the decomposition
identity residual never exceeded 5.56e-17.

| Bootstrap scheme | Pointwise coverage range | Median pointwise coverage | Uniform coverage range | Median uniform coverage |
|---|---:|---:|---:|---:|
| Multiplier | 0.833-0.967 | 0.933 | 0.833-0.967 | 0.900 |
| Counterfactual | 0.967-1.000 | 1.000 | 0.967-1.000 | 1.000 |

These figures validate that pointwise and uniform inference, checkpointing,
and failure recovery execute end to end. They are not publication-grade
coverage evidence: 30 outer replications and 49 bootstrap draws are too small
for a precise 95% coverage assessment. In particular, the multiplier intervals
appear somewhat liberal in some cells and the counterfactual intervals appear
conservative. A methodological release should repeat this exercise with a
larger outer loop and at least 199, preferably 499, bootstrap draws.

## Remaining 1.0.0 blockers

- Replace the placeholder author and maintainer metadata.
- Add the public repository URL, bug-report URL, and package citation.
- Run and archive the licensed Stata cross-language parity test.
- Execute the Windows, Linux, macOS, and self-hosted CUDA workflows in CI.
- Run the larger Monte Carlo coverage study described above.

Version 0.15.0 is therefore a validated release candidate, not a final 1.0.0
release.
