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

QR and CQR marginal quantile selection uses the same stable normalized
weighted type-7 rule on CPU and GPU. It agrees with
`Hmisc::wtd.quantile(..., normwt = TRUE)` for ordinary weights when cumulative
weighted ranks are numerically distinct, but evaluates ranks directly to avoid
`approx()` tie collapse for extreme weight ranges. Explicit
`marginal_method = "chunked"` continues to request the bounded-host-memory
implementation; `auto` and `matrix` use GPU-resident sorting when CUDA is
selected. The legacy scale-sensitive definition remains on the CPU. GPU
sorting must fit the conditional draws and sorting workspace in device memory.
Bootstrap replications remain
checkpointed and deterministic and run in one R process; threshold and
quantile columns are batched inside the GPU process.

## Runtime setup

CuPy and `reticulate` are optional dependencies. A project-local Python
installation can be passed through `gpu_python_path`. On Windows, the package
adds the CUDA runtime directories supplied by the Python packages before CuPy
is imported. Use an isolated virtual environment; do not mix a `--target`
installation with packages visible in the Python user site, which can load
CuPy from one location and CUDA libraries from another.
`gpu_backend_status()` reports
the selected Python executable, Python/NumPy/CuPy versions and files, device,
CUDA runtime and driver versions, visible CuPy distributions, source-module
SHA-256, and initialization warnings. It marks CUDA available only after small
matrix-multiplication, sorting, and linear-solve probes have completed and the
default stream has synchronized.

CUDA modules are loaded under a name derived from their normalized path and
content hash, and their Python `__file__` is checked against the requested
file. Imported modules and normalized runtime paths are cached once per unique
runtime/module identity. Existing `PATH` and `PYTHONPATH` entries are split on
the platform path separator and de-duplicated, so repeated bootstrap and chunk
calls do not grow the process environment. A custom `gpu_module_path` executes
Python code and must therefore point only to a trusted local file.

The exact Windows/Python 3.12 validation set is bundled at
`inst/python/requirements-cuda-windows-py312.lock`. For example:

```powershell
py -3.12 -m venv tmp/cuda-venv
$env:PYTHONNOUSERSITE = "1"
tmp/cuda-venv/Scripts/python.exe -m pip install `
  -r inst/python/requirements-cuda-windows-py312.lock
```

Pass the virtual environment's Python executable as `gpu_python` and its
`Lib/site-packages` directory as `gpu_python_path`. The shorter
`requirements-cuda.txt` is the normal minimum entry point; the lock file is the
release validation environment. Other compatible CUDA 12 devices and Python
versions may work, but must pass `tests/gpu_backend.R` before substantive
analysis. Setting `PYTHONNOUSERSITE=1` before R initializes Python prevents a
user-site CuPy installation from shadowing the selected environment.

Model/backend combinations are validated before execution. CUDA distribution
regression is restricted to logit, probit, and cloglog; Cox remains a CPU-only
fit and marginalization path. Run metadata distinguishes the requested backend
from the resolved fit backend and reports fit, prediction, and marginalization
devices independently. CUDA DR additionally stores each threshold's actual
backend and any CPU-fallback reason. Point-run metadata retain the selected
Python, NumPy, and CuPy versions and files, device and CUDA versions, module
path and SHA-256, isolation flags, capability-probe results, and runtime
warnings so the executed GPU environment remains auditable.

Checkpoint/runtime identity includes Python, NumPy and CuPy versions, device
and CUDA versions, the loaded module path and SHA-256, capability-probe
results, and visible-runtime warnings. A change in this identity prevents a
checkpoint produced by another GPU implementation from being silently reused.

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

When `float32` is explicitly selected, DR probability clipping uses a
float32-aware epsilon. CPU and CUDA boundary diagnostics use the same
scale-aware signed-margin condition for complete and quasi-complete separation.
This avoids declaring a model separated solely because one legitimate fitted
probability is extreme. Survey-weight cumulative sums for CUDA QR
marginalization remain float64 even when conditional draws use float32.

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
