# scalableCounterfactual 1.0.0 release checklist

## Required before tagging

- [x] Replace the placeholder `Authors@R` and maintainer email.
- [x] Add package `CITATION`.
- [x] Add public repository `URL` and `BugReports`.
- [x] Exclude Python bytecode and cache directories from the source tarball.
- [x] Record the validated CUDA Python requirements and exact Windows lock.
- [x] Freeze the public API and output schema contract.
- [x] Run `R CMD check` on a clean source tarball.
- [x] Run the optional CUDA parity suite on a real NVIDIA device.
- [ ] Run the optional Stata cross-language parity test and archive its output.
- [x] Run CPU CI on Windows, Linux, and macOS.
- [x] Complete a separate end-to-end point-estimation run on macOS.
- [ ] Run the self-hosted CUDA CI workflow.
- [x] Review simulation bias/coverage and generic solver benchmark reports.
- [ ] Run publication-grade coverage validation with a larger outer loop and
      at least 199 bootstrap draws.
- [x] Install the final 1.0.0 tarball into the analysis R library.
- [x] Record the final 1.0.0 tarball hash in the external release bundle.

## Release evidence

The release directory should retain check logs, session information, CUDA
status, dependency versions, simulation summaries, benchmark summaries,
Stata parity output, the source-tarball SHA-256 hash, and the complete output
of the release preflight script.

To create the Stata evidence file on a licensed installation:

```powershell
$env:STATA_EXE = "C:/Program Files/Stata18/StataMP-64.exe"
$env:STATA_PARITY_OUTPUT = "inst/provenance/stata_parity_verified.csv"
Rscript tests/stata_parity.R
```

To generate the generic analytic-simulation and solver benchmark evidence:

```powershell
Rscript tools/run_release_validation.R `
  --workers 4 `
  --output output/release_validation
```

Use `--quick` only for development checks. The release evidence must use the
default 100 simulation replications and three timed benchmark repetitions.

The separate nested-bootstrap coverage check is intentionally smaller because
each outer replication contains a complete bootstrap run:

```powershell
Rscript tools/run_release_coverage.R `
  --workers 4 `
  --replications 30 `
  --bootstrap-reps 49 `
  --bootstrap-scheme counterfactual `
  --output output/release_validation/coverage
```

These 30 outer replications are a release smoke check for inference plumbing,
not a publication-grade Monte Carlo coverage study. A methodological release
should increase both the outer and bootstrap replication counts.
