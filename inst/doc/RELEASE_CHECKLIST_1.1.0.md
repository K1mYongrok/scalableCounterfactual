# scalableCounterfactual 1.1.0 release checklist

## Required before tagging

- [x] Preserve the 1.0 public API as an ordered prefix and version the additive
      `dr_noncrossing` control and `fit_weighted_qr()` frequency argument.
- [x] Preserve the 1.0 numerical defaults for the QR scalar shift and DR
      running-maximum correction.
- [x] Version and test the extended marginalization-diagnostics schema.
- [x] Align matrix, chunked, and CUDA marginalization on the stable weighted
      type-7 rank rule, including extreme weight ranges.
- [x] Harden CQR selection, DR monotonicity correction, bootstrap recovery,
      functional inference, and checkpoint identities.
- [x] Isolate and lock the optional CUDA CI environment and record runtime
      identity and capability probes.
- [x] Synchronize `DESCRIPTION`, `CITATION`, the public repository links, and
      the external `release/BUILD_INFO.txt` metadata with version 1.1.0.
- [x] Run `R CMD check --as-cran --no-manual` on a clean 1.1.0 source tarball.
- [x] Run the complete CPU regression suite from an isolated installed library.
- [x] Run the CUDA parity suite on a real NVIDIA device using the final source.
- [x] Confirm the public CPU CI matrix on Windows, Linux, and macOS.
- [x] Create and review the final source tarball and SHA-256 before tagging.
- [x] Regenerate and visually inspect `release/scalableCounterfactual-reference.pdf`
      from the final 1.1.0 source.
- [x] Record the final evidence and remaining limitations in
      `RELEASE_VALIDATION_1.1.0.md`.

## Optional cross-language evidence

- [ ] Run the Stata parity harness and archive its output before making
      cross-language equivalence claims. Enforce this optional gate with
      `tools/release_check.R --require-stata`.

The checklist records release gates; an unchecked item is not evidence that a
feature is broken. No 1.1.0 tag or GitHub release should be created until the
required gates are complete and the release-validation record is updated.
The optional Stata item is required only for a Stata-equivalence claim.
