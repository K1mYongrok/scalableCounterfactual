# Installed examples

- `quick_start.R`: self-contained CPU QR decomposition, bootstrap, summary,
  plot, and CSV output example.
- `cuda_status.R`: checks whether the optional CUDA backend can be resolved and
  runs a small CUDA marginalization example when available.

After installation:

```r
example_dir <- system.file("examples", package = "scalableCounterfactual")
source(file.path(example_dir, "quick_start.R"))
source(file.path(example_dir, "cuda_status.R"))
```

`quick_start.R` writes its figure and CSV/RDS outputs to
`scalableCounterfactual_quick_start_output` under the current working
directory.
