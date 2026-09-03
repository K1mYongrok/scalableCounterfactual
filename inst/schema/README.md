# Versioned public contracts

This directory records the versioned public contracts in a form that can be
reviewed without reading the implementation.

- `public_api_1.0.csv` preserves the initial 1.0 API. `public_api_1.1.csv`
  appends the distribution-regression noncrossing control to `cf_control()`
  and optional row frequencies to `fit_weighted_qr()` while preserving all
  1.0 functions and ordered arguments. Methods registered only through S3 are
  tested separately.
- `s3_methods_1.0.csv` preserves the initial decomposition-method contract.
  `s3_methods_1.1.csv` also records the existing
  `print.cf_simulation_validation()` signature, which was registered in 1.0
  but omitted from the initial machine-readable inventory.
- `decomposition_1.0.csv` preserves the initial decomposition table;
  `decomposition_1.1.csv` confirms that its ordered columns and meanings are
  unchanged in schema 1.1.
- `output_files_1.0.csv` preserves the initial file set;
  `output_files_1.1.csv` confirms that required and conditional files are
  unchanged in schema 1.1.
- `marginalization_1.1.csv` fixes the ordered marginalization diagnostics,
  including distribution-regression crossing and probability-bound checks.

The executable contract is `tests/api_contract.R`. A change to these files
must be accompanied by a deliberate schema/API version decision; updating the
test merely to make an unintended interface change pass is not acceptable.
