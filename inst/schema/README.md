# Versioned public contracts

This directory records the stable 1.0 public contracts in a form that can
be reviewed without reading the implementation.

- `public_api_1.0.csv` fixes the exported function names and ordered formal
  arguments. Methods registered only through S3 are tested separately.
- `s3_methods_1.0.csv` fixes the registered public method signatures.
- `decomposition_1.0.csv` fixes the ordered columns written to
  `decomposition.csv` and their statistical meaning.
- `output_files_1.0.csv` distinguishes files written on every run from files
  written only when the corresponding calculation exists.

The executable contract is `tests/api_contract.R`. A change to these files
must be accompanied by a deliberate schema/API version decision; updating the
test merely to make an unintended interface change pass is not acceptable.
