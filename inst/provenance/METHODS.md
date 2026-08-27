# Method and source provenance

This file distinguishes statistical references from implementation ancestry.
“Method” means that the package implements the cited estimand or algorithm;
“delegated” means that numerical fitting is performed by the named R package;
“transcribed” means that source logic was translated into native R. No source
translation is claimed unless it is stated explicitly.

| Component | Package implementation | Statistical or algorithmic source | Software source and relationship |
|---|---|---|---|
| Counterfactual distributions, structure/composition effects, functional inference | `R/decomposition_engine.R`, `R/bootstrap.R`, `R/functional_inference.R` | Chernozhukov, Fernández-Val, and Melly (2013), [doi:10.3982/ECTA10582](https://doi.org/10.3982/ECTA10582) | The estimands and bootstrap organization follow the paper and CRAN `Counterfactual` 1.2; the package is a new modular implementation, not a source copy. |
| Linear quantile regression objective | `R/solver_qr.R`, `R/model_qr.R` | Koenker and Bassett (1978), [doi:10.2307/1913643](https://doi.org/10.2307/1913643) | Exact numerical fits are delegated to public routines in CRAN [`quantreg`](https://cran.r-project.org/package=quantreg). Version 6.1 was used for release validation. |
| Frisch–Newton and preprocessing QR solvers (`fn`, `pfn`, `qfnb`, `pfnb`, `proqreg`, `profn`) | `R/solver_qr.R` | Portnoy and Koenker (1997), [doi:10.1214/ss/1030037960](https://doi.org/10.1214/ss/1030037960); solver-specific details are documented by `quantreg` | Delegated to `quantreg::rq.fit.fnb()`, `rq.fit.pfn()`, `rq.fit.qfnb()`, `rq.fit.pfnb()`, and `rq.fit.ppro()`; no `quantreg` source is copied. |
| Barrodale–Roberts QR (`br`) | `R/solver_qr.R` | Barrodale and Roberts (1974), Algorithm 478, [doi:10.1145/355616.361024](https://doi.org/10.1145/355616.361024) | Delegated to `quantreg::rq.fit.br()`; no source is copied. |
| One-step QR process and bootstrap (`onestep`) | `R/solver_qr.R`, `R/bootstrap.R` | Chernozhukov, Fernández-Val, and Melly (2022), [doi:10.1007/s00181-020-01898-0](https://doi.org/10.1007/s00181-020-01898-0) | Transcribed from MIT-licensed `qrprocess.ado` 1.1.3 at commit `ec56830ef9c84ce54411ad59c5ce94535847d9df`; exact lines and hash are in `qrprocess_onestep.yml`. |
| Censored quantile regression (`cqr`) | `R/model_cqr.R` | Chernozhukov and Hong (2002), [doi:10.1198/016214502388618663](https://doi.org/10.1198/016214502388618663) | Selection stages are implemented in this package; QR stages are delegated to the selected exact `quantreg` solver. |
| Distribution regression (`logit`, `probit`, `cloglog`, `lpm`) | `R/model_dr.R`, `R/model_location.R` | Counterfactual-distribution framework: Chernozhukov, Fernández-Val, and Melly (2013), [doi:10.3982/ECTA10582](https://doi.org/10.3982/ECTA10582) | Base fits use `stats::glm.fit()` or package WLS; optional `fastglm` and `speedglm` backends are delegated when requested. |
| Cox conditional distribution and Kaplan–Meier reference quantiles | `R/model_cox.R` | Cox (1972), [doi:10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x); Kaplan and Meier (1958), [doi:10.1080/01621459.1958.10501452](https://doi.org/10.1080/01621459.1958.10501452) | Cox fitting is delegated to `survival::coxph.fit()`; weighted Breslow and Kaplan–Meier distribution calculations are implemented in this package. |
| Quantile rearrangement | `R/marginal_distribution.R` | Chernozhukov, Fernández-Val, and Galichon (2010), [doi:10.3982/ECTA7880](https://doi.org/10.3982/ECTA7880) | Row-wise rearrangement is implemented in this package. |
| QR residual-sign bootstrap preprocessing (`xy_preprocess`) | `R/bootstrap.R` | Algorithm exposed by CRAN `quantreg` | Preprocessing is delegated to `quantreg::boot.rq.pxy()`; no source is copied. |
| Experimental CUDA ADMM QR (`cuda_admm`) | `R/gpu_backend.R`, `inst/python/scalablecf_cuda.py` | ADMM background: Boyd et al. (2011), [doi:10.1561/2200000016](https://doi.org/10.1561/2200000016) | Package-specific CuPy implementation. It is marked approximate and experimental; it is not a translation of `quantreg` or `qrprocess`. |
| CUDA prediction and marginalization | `R/gpu_backend.R`, `inst/python/scalablecf_cuda.py` | Same counterfactual estimand as the CPU path | Package-specific CuPy implementation; CUDA changes the execution device, not the documented estimand. |

The package `DESCRIPTION`, `inst/NOTICE`, and source comments point here. The
release benchmark records the runtime package versions separately because an
algorithm citation is not a reproducible software environment specification.
