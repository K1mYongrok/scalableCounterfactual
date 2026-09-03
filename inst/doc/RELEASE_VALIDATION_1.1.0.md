# scalableCounterfactual 1.1.0 validation record

Validation date: 2026-09-04 (Asia/Seoul)

This record describes the evidence produced from the final 1.1.0 source tree
before tagging. It is not a claim of bitwise equivalence across operating
systems, numerical libraries, or approximate solvers.

## Local build and package check

- Platform: Windows 11 x64, R 4.5.2 (UCRT), quantreg 6.1.
- A clean source tarball was built with `R CMD build` and checked with
  `R CMD check --as-cran --no-manual`.
- Result: 0 errors, 1 warning, and 1 note.
- The warning was environmental: the optional `qpdf` executable was not
  installed, so R could not run its PDF size-reduction check.
- The note was the expected CRAN incoming-check classification for a new
  submission and maintainer.
- Package installation, loading, namespace checks, examples, documentation,
  code checks, and every packaged test completed successfully.

## CPU regression suite

The final source was installed into an isolated library. The following test
programs passed:

- `api_contract.R`
- `reference_parity.R`
- `smoke.R`
- `review_fixes_core.R`
- `review_fixes_inference.R`
- `benchmark.R`
- `extension_api.R`
- `gpu_backend.R` in its CPU/skip configuration
- `review_fixes_gpu_packaging.R` in its CPU/skip configuration
- `stata_parity.R` in its documented optional/skip configuration

These tests cover the versioned API and output schemas, weighted quantile and
frequency-compressed resampling rules, QR/CQR/DR/Cox paths, decomposition
identities, bootstrap recovery and functional inference, custom extensions,
transactional output writing, command-line adapters, and packaging rules.

## Real-device CUDA validation

The CUDA test programs `gpu_backend.R` and
`review_fixes_gpu_packaging.R` also passed with a real CUDA device and an
isolated Python 3.12 environment. They exercised capability detection,
CPU/CUDA marginalization parity, distribution-regression execution,
frequency-aware QR handling, runtime isolation, and packaging checks. CUDA
remains optional; the CPU package neither imports Python nor requires CuPy.

## Simulation and solver validation

The quick release-validation workflow completed 30 of 30 Monte Carlo tasks
across the location-shift, composition-shift, and combined scale-shift
scenarios, with zero failed tasks, zero retries, and a maximum decomposition
identity residual of approximately `5.55e-17`.

All requested QR solver benchmark runs completed. Against BR on the frozen
quick-validation design, the largest absolute decomposition-effect
differences were approximately:

- FN: `9.21e-09`
- QFNB: `9.20e-09`
- PFNB: `1.99e-07`

The `onestep` result was also produced, but its larger difference is expected:
it is an explicitly approximate process estimator and is not included in an
exact-solver equivalence claim. These quick runs validate execution and basic
numerical consistency; they are not publication-quality speed benchmarks.

## Reference manual and release artifacts

- The package reference manual was regenerated from the 1.1.0 source.
- All 20 rendered pages and representative full-resolution pages were visually
  inspected. No clipping, overlap, missing content, or stale version metadata
  was found.
- The final source archive is accompanied by a SHA-256 file in the external
  `release/` directory. That directory is intentionally excluded from the
  source package and Git repository.

## Public CPU CI

GitHub Actions run
[`33808346987`](https://github.com/K1mYongrok/scalableCounterfactual/actions/runs/33808346987)
completed successfully for source commit `4f37ecb` on Windows, Ubuntu Linux,
and macOS with the release R version. This matrix reran the package preflight,
installation, examples, and packaged regression tests independently on all
three platforms.

## Remaining external gates

- The Stata parity harness remains optional and was not run locally because no
  Stata executable was configured. No cross-language numerical-equivalence
  claim is made without that evidence.
- No 1.1.0 tag or GitHub Release was created as part of this validation.
