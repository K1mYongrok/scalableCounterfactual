# Release validation for 1.0.0

Version 1.0.0 was built and checked locally on 2026-08-27 with R 4.5.2 on
64-bit Windows. The package author and maintainer are Yongrok Kim, and the
source package includes a machine-readable `inst/CITATION` entry.

## Final source-package check

- `R CMD build` completed successfully from the package source tree.
- The resulting source tarball passed `R CMD check --no-manual` with status
  `OK` under a native Windows UTF-8 locale.
- All examples and the API-contract, benchmark, extension-API, GPU-backend,
  frozen-reference parity, smoke, and Stata-harness tests completed.
- The public API and output schema remain frozen at schema version `1.0`.

The first aborted check was caused by inherited `LC_ALL=C.UTF-8` and
`LC_CTYPE=C.UTF-8` environment variables, which are not valid locale names in
this Windows R installation. Removing those shell variables restored the
native `Korean_Korea.utf8` locale and the same tarball passed without warnings
or errors. This was an execution-environment issue, not a package failure.

## Inherited numerical evidence

No estimator implementation changed between 0.15.0 and 1.0.0. The 0.15.0
release evidence therefore continues to apply: 300 analytic simulation fits
completed, the maximum reported absolute bias was 0.0163 log point, CPU/CUDA
parity passed on an NVIDIA GeForce RTX 4060 Ti, and the bootstrap smoke checks
preserved the decomposition identity to within 5.56e-17. See
`RELEASE_VALIDATION_0.15.0.md` for the complete numerical tables.

## Distribution status and deferred evidence

This is a stable source release and has not been submitted to CRAN. Its public
source repository is `https://github.com/K1mYongrok/scalableCounterfactual`,
and bug reports are accepted through the repository issue tracker.

The following evidence remains appropriate before making broader public or
methodological claims:

- licensed cross-language execution of the optional Stata parity harness;
- remote CPU CI execution on Windows, Linux, and macOS;
- execution of the self-hosted CUDA CI workflow; and
- a publication-grade Monte Carlo coverage study with larger outer and
  bootstrap replication counts.

These deferred checks do not affect the successful local R package check, but
they should be completed before a CRAN submission or a methodological paper
claims comprehensive cross-platform validation.
