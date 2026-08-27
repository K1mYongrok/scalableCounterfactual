# Optional CUDA backend

## Scope

The CUDA layer changes the computation device, not the decomposition estimand.
The default remains `gpu_backend = "cpu"`. The optional layer contains:

- CuPy float64/float32 dense matrix multiplication for QR, CQR, location, and
  location-scale predictions;
- GPU-resident sorting, cumulative survey weighting, and marginal-quantile
  selection for QR and CQR when `legacy_weighted_quantile = FALSE`;
- GPU-resident probability aggregation for distribution regression;
- batched unpenalized binomial IRLS for logit, probit, and cloglog DR;
- an experimental weighted QR ADMM solver named `cuda_admm`.

QR and CQR marginal quantile selection follows the established
`Hmisc::wtd.quantile(..., normwt = TRUE)` rank and interpolation definition on
the GPU. Explicit `marginal_method = "chunked"` continues to request the
bounded-host-memory implementation; `auto` and `matrix` use GPU-resident
sorting when CUDA is selected. The legacy scale-sensitive definition remains
on the CPU. GPU sorting must fit the conditional draws and sorting workspace in
device memory. Bootstrap replications remain
checkpointed and deterministic and run in one R process; threshold and
quantile columns are batched inside the GPU process.

## Runtime setup

CuPy and `reticulate` are optional dependencies. A project-local Python
installation can be passed through `gpu_python_path`. On Windows, the package
adds the CUDA runtime directories supplied by the Python packages before CuPy
is imported. `gpu_backend_status()` reports the selected Python executable,
device, CuPy version, and initialization error. Imported modules and normalized
runtime paths are cached once per Python executable. Existing `PATH` and
`PYTHONPATH` entries are split on the platform path separator and de-duplicated,
so repeated bootstrap and chunk calls do not grow the process environment.

The validated installation entry point is bundled at
`inst/python/requirements-cuda.txt`:

```powershell
python -m pip install --target tmp/cuda_python `
  -r inst/python/requirements-cuda.txt
```

Release 0.15.0 was validated with Python 3.12.1, NumPy 2.5.1,
CuPy CUDA 12x 14.1.1, and an NVIDIA RTX 4060 Ti. The exact Windows dependency
set is recorded in `requirements-cuda-windows-py312.lock`. Other compatible
CUDA 12 devices and Python versions may work, but must pass
`tests/gpu_backend.R` before substantive analysis.

Model/backend combinations are validated before execution. CUDA distribution
regression is restricted to logit, probit, and cloglog; Cox remains a CPU-only
fit and marginalization path. Run metadata distinguishes the requested backend
from the resolved fit backend and reports fit, prediction, and marginalization
devices independently. CUDA DR additionally stores each threshold's actual
backend and any CPU-fallback reason.

## Numerical policy

`float64` is the default. A candidate backend should be assessed using the
same sample, design matrix, weights, threshold/quantile grid, seed, and
bootstrap law. Record at least:

1. coefficient differences for non-boundary, identified fits;
2. fitted-probability and objective differences at boundary fits;
3. counterfactual distribution and decomposition-effect differences;
4. bootstrap standard-error differences;
5. raw and corrected quantile-crossing diagnostics;
6. elapsed time, R heap, device memory, and precision.

Large coefficient differences under complete or quasi-complete separation do
not by themselves establish a numerical disagreement: the finite coefficient
vector is not identified there. Prediction and objective comparisons are the
primary diagnostics for those thresholds.

## Current limitations

- `cuda_admm` is experimental, approximate, and not selected by `auto`.
- `cuda_admm` reports its internal preconditioner and evaluates convergence at
  the returned final iterate; a column that met tolerance earlier is not frozen
  or automatically declared converged.
- Independent bootstrap replications are not yet fused into one device batch.
- Multiple R GPU workers are rejected to prevent device oversubscription.
- CUDA is beneficial only beyond a hardware- and design-specific break-even
  sample size because Python initialization and host/device transfers are
  material for small samples.
