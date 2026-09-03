"""Optional CuPy kernels for scalableCounterfactual."""

from __future__ import annotations

import hashlib
import importlib.metadata
from pathlib import Path
import site
import sys
import time
import numpy as np

try:
    import cupy as cp
except Exception as exc:  # reported by cuda_status
    cp = None
    _IMPORT_ERROR = repr(exc)
    _IMPORT_ERROR_TYPE = type(exc).__name__
else:
    _IMPORT_ERROR = None
    _IMPORT_ERROR_TYPE = None


_STATUS_CACHE = None


def _module_sha256():
    try:
        return hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    except (OSError, TypeError):
        return None


def _cupy_distribution_records():
    """Report all visible CuPy distributions without importing extra packages."""
    records = []
    warnings = []
    try:
        for distribution in importlib.metadata.distributions():
            name = str(distribution.metadata.get("Name", ""))
            canonical = name.lower().replace("_", "-")
            if canonical != "cupy" and not canonical.startswith("cupy-cuda"):
                continue
            try:
                location = str(Path(distribution.locate_file("")).resolve())
            except (OSError, TypeError, ValueError):
                location = str(distribution.locate_file(""))
            records.append(
                f"{name}=={distribution.version} @ {location}"
            )
    except (importlib.metadata.PackageNotFoundError, OSError, ValueError) as exc:
        warnings.append(
            "could not enumerate installed CuPy distributions: "
            f"{type(exc).__name__}: {exc}"
        )
    records = sorted(set(records))
    if len(records) > 1:
        warnings.append(
            "multiple CuPy distributions are visible; use an isolated virtual "
            "environment and set PYTHONNOUSERSITE=1"
        )
    return records, warnings


def _runtime_metadata():
    distributions, warnings = _cupy_distribution_records()
    try:
        user_site = site.getusersitepackages()
    except (AttributeError, TypeError):
        user_site = None
    if isinstance(user_site, (list, tuple)):
        user_sites = [str(Path(value).resolve()) for value in user_site]
    elif user_site:
        user_sites = [str(Path(user_site).resolve())]
    else:
        user_sites = []
    cupy_file = getattr(cp, "__file__", None) if cp is not None else None
    if cupy_file and bool(site.ENABLE_USER_SITE):
        resolved_cupy = str(Path(cupy_file).resolve())
        if any(
            resolved_cupy == value or resolved_cupy.startswith(value + str(Path("/")))
            for value in user_sites
        ):
            warnings.append(
                "CuPy was imported from the Python user site; an explicitly "
                "selected environment may be shadowed"
            )
    return {
        "python_version": sys.version.split()[0],
        "python_executable": sys.executable,
        "numpy_version": np.__version__,
        "numpy_file": getattr(np, "__file__", None),
        "cupy_version": getattr(cp, "__version__", None),
        "cupy_file": cupy_file,
        "cupy_distributions": distributions,
        "python_isolated": bool(sys.flags.isolated),
        "python_no_user_site": bool(sys.flags.no_user_site),
        "user_site_enabled": bool(site.ENABLE_USER_SITE),
        "module_file": str(Path(__file__).resolve()),
        "module_sha256": _module_sha256(),
        "runtime_warnings": warnings,
    }


def cuda_status(refresh=False):
    """Probe the CUDA capabilities used by this module, once per process."""
    global _STATUS_CACHE
    if _STATUS_CACHE is not None and not refresh:
        return dict(_STATUS_CACHE)
    metadata = _runtime_metadata()
    status = {
        **metadata,
        "available": False,
        "error": _IMPORT_ERROR,
        "error_type": _IMPORT_ERROR_TYPE,
        "device_id": None,
        "device": None,
        "cuda_runtime_version": None,
        "cuda_driver_version": None,
        "capability_matmul": False,
        "capability_sort": False,
        "capability_solve": False,
    }
    if cp is None:
        _STATUS_CACHE = status
        return dict(status)
    try:
        if int(cp.cuda.runtime.getDeviceCount()) < 1:
            raise RuntimeError("no CUDA devices")
        device_id = int(cp.cuda.runtime.getDevice())
        props = cp.cuda.runtime.getDeviceProperties(device_id)
        name = props["name"]
        if isinstance(name, bytes):
            name = name.decode("utf-8")
        status.update({
            "device_id": device_id,
            "device": name,
            "cuda_runtime_version": int(cp.cuda.runtime.runtimeGetVersion()),
            "cuda_driver_version": int(cp.cuda.runtime.driverGetVersion()),
        })
        left = cp.asarray([[1.0, 2.0], [3.0, 4.0]], dtype=cp.float64)
        right = cp.asarray([[2.0], [1.0]], dtype=cp.float64)
        product = left @ right
        cp.cuda.Stream.null.synchronize()
        if not bool(cp.allclose(product, cp.asarray([[4.0], [10.0]]))):
            raise RuntimeError("CUDA matrix multiplication returned invalid values")
        status["capability_matmul"] = True
        ordered = cp.sort(cp.asarray([3.0, 1.0, 2.0], dtype=cp.float64))
        cp.cuda.Stream.null.synchronize()
        if not bool(cp.allclose(ordered, cp.asarray([1.0, 2.0, 3.0]))):
            raise RuntimeError("CUDA sort returned invalid values")
        status["capability_sort"] = True
        solution = cp.linalg.solve(
            cp.asarray([[2.0, 0.0], [0.0, 4.0]], dtype=cp.float64),
            cp.asarray([2.0, 8.0], dtype=cp.float64),
        )
        cp.cuda.Stream.null.synchronize()
        if not bool(cp.allclose(solution, cp.asarray([1.0, 2.0]))):
            raise RuntimeError("CUDA linear solve returned invalid values")
        status["capability_solve"] = True
        status.update({
            "available": True,
            "error": None,
            "error_type": None,
        })
    except Exception as exc:
        # This is a capability probe: preserve the exact failure class instead
        # of silently treating a partially usable runtime as available.
        status["error"] = f"{type(exc).__name__}: {exc}"
        status["error_type"] = type(exc).__name__
    _STATUS_CACHE = status
    return dict(status)


def _require_cuda():
    status = cuda_status()
    if not status["available"]:
        raise RuntimeError("CUDA is unavailable: " + str(status["error"]))


def _dtype(precision):
    if precision == "float64":
        return cp.float64
    if precision == "float32":
        return cp.float32
    raise ValueError("precision must be 'float64' or 'float32'")


def _probability_epsilon(dtype):
    return dtype(1e-12 if dtype == cp.float64 else 1e-6)


def _logit_eta_limit(dtype):
    return dtype(709.0 if dtype == cp.float64 else 80.0)


def matmul(x, coefficients, precision="float64"):
    _require_cuda()
    dtype = _dtype(precision)
    result = cp.asarray(x, dtype=dtype, order="C") @ cp.asarray(
        coefficients, dtype=dtype, order="C"
    )
    cp.cuda.Stream.null.synchronize()
    return cp.asnumpy(result).astype(np.float64, copy=False)


def qr_marginal_quantiles(x, coefficients, weights, probs, shift=0.0,
                          precision="float64", normalization_rows=None):
    """Pool conditional QR draws and select weighted quantiles on the GPU.

    The R wrapper supplies finite positive weights normalized to the effective
    row count. This may exceed the stored row count when bootstrap duplicates
    are represented by frequency weights.
    The rank positions and interpolation rule match the package's stable
    weighted type-7 rule (and ordinary-range Hmisc normwt=TRUE results). Only
    the selected quantiles are copied back to the host.
    """
    _require_cuda()
    dtype = _dtype(precision)
    xg = cp.asarray(x, dtype=dtype, order="C")
    bg = cp.asarray(coefficients, dtype=dtype, order="C")
    # Keep cumulative weights in float64 even when predictions use float32.
    wg = cp.asarray(weights, dtype=cp.float64)
    pg = cp.atleast_1d(cp.asarray(probs, dtype=cp.float64))
    draws = xg @ bg
    if float(shift) != 0.0:
        draws += dtype(shift)
    if not bool(cp.all(cp.isfinite(draws))):
        raise ValueError("conditional draws contain non-finite values")

    n, q = draws.shape
    stored_draws = int(n) * int(q)
    effective_rows = float(n) if normalization_rows is None else float(
        normalization_rows
    )
    total_positions = effective_rows * float(q)
    values = draws.ravel(order="F")
    ordering = cp.argsort(values)
    sorted_values = values[ordering]
    sorted_weights = wg[ordering % int(n)]
    cp.cumsum(sorted_weights, out=sorted_weights)
    # The R boundary normalizes row weights to the effective row count. Pin the
    # analytically known last rank so selection is independent of GPU summation
    # order.
    sorted_weights[-1] = cp.float64(total_positions)

    rank = 1.0 + (float(total_positions) - 1.0) * pg
    low = cp.maximum(cp.floor(rank), 1.0)
    high = cp.minimum(low + 1.0, float(total_positions))
    interpolation = cp.mod(rank, 1.0)
    total_weight = cp.float64(total_positions)
    targets = cp.concatenate((low, high))
    # CuPy accumulates these weights in float64, but the sort order can still
    # change the last few bits of a cumulative sum.  Use the same rank-hit
    # convention as the CPU paths so numerically equal boundaries select the
    # same observation after matrix/chunk/GPU reordering.
    rank_tolerance = (
        cp.float64(64.0 * np.finfo(np.float64).eps)
        * cp.maximum(cp.float64(1.0), total_weight)
    )
    selected_index = cp.searchsorted(
        sorted_weights, targets - rank_tolerance, side="left"
    )
    selected_index = cp.minimum(selected_index, stored_draws - 1)
    selected = sorted_values[selected_index]
    k = int(pg.size)
    result = ((1.0 - interpolation) * selected[:k] +
              interpolation * selected[k:])
    cp.cuda.Stream.null.synchronize()
    return cp.asnumpy(result).astype(np.float64, copy=False)


def dr_marginal_cdf(x, coefficients, weights, link, precision="float64",
                    block_columns=16):
    """Compute X beta and its weighted probability average on the GPU."""
    _require_cuda()
    dtype = _dtype(precision)
    xg = cp.asarray(x, dtype=dtype, order="C")
    bg = cp.asarray(coefficients, dtype=dtype, order="C")
    wg = cp.asarray(weights, dtype=dtype)
    output = cp.empty(bg.shape[1], dtype=dtype)
    denominator = cp.sum(wg)
    for first in range(0, bg.shape[1], int(block_columns)):
        last = min(first + int(block_columns), bg.shape[1])
        eta = xg @ bg[:, first:last]
        if link == "logit":
            limit = _logit_eta_limit(dtype)
            probability = 1.0 / (1.0 + cp.exp(-cp.clip(eta, -limit, limit)))
        elif link == "probit":
            probability = 0.5 * cp.erfc(-eta / cp.sqrt(dtype(2.0)))
        elif link == "cloglog":
            probability = 1.0 - cp.exp(-cp.exp(cp.clip(eta, -36.0, 36.0)))
        elif link == "lpm":
            probability = eta
        else:
            raise ValueError("unsupported DR link: " + str(link))
        output[first:last] = cp.sum(probability * wg[:, None], axis=0) / denominator
    cp.cuda.Stream.null.synchronize()
    return cp.asnumpy(output).astype(np.float64, copy=False)


def _link_values(eta, link, dtype):
    if link == "logit":
        limit = _logit_eta_limit(dtype)
        mu = 1.0 / (1.0 + cp.exp(-cp.clip(eta, -limit, limit)))
        derivative = mu * (1.0 - mu)
    elif link == "probit":
        mu = 0.5 * cp.erfc(-eta / cp.sqrt(dtype(2.0)))
        derivative = cp.exp(-0.5 * eta * eta) / cp.sqrt(dtype(2.0 * np.pi))
    elif link == "cloglog":
        exp_eta = cp.exp(cp.clip(eta, -36.0, 36.0))
        mu = 1.0 - cp.exp(-exp_eta)
        derivative = cp.exp(eta - exp_eta)
    else:
        raise ValueError("unsupported DR link: " + str(link))
    eps = _probability_epsilon(dtype)
    return cp.clip(mu, eps, 1.0 - eps), cp.maximum(derivative, eps)


def fit_dr_process(x, y, weights, thresholds, link, precision="float64",
                   maxit=100, tolerance=1e-8, block_thresholds=16):
    """Batched unpenalized binomial IRLS with X resident on the GPU."""
    _require_cuda()
    dtype = _dtype(precision)
    xg = cp.asarray(x, dtype=dtype, order="C")
    yg = cp.asarray(y, dtype=dtype)
    wg = cp.asarray(weights, dtype=dtype)
    # reticulate may simplify a length-one R vector to a Python scalar.  Keep
    # the process dimension explicit so a single regular threshold follows the
    # same batched path as longer grids.
    tg = cp.atleast_1d(cp.asarray(thresholds, dtype=dtype))
    _, p = xg.shape
    total = int(tg.size)
    coefficients = cp.empty((p, total), dtype=dtype)
    iterations = np.zeros(total, dtype=np.int32)
    converged = np.zeros(total, dtype=bool)
    boundary = np.zeros(total, dtype=bool)
    started = time.perf_counter()

    for first in range(0, total, int(block_thresholds)):
        last = min(first + int(block_thresholds), total)
        response = (yg[:, None] <= tg[None, first:last]).astype(dtype)
        width = last - first
        beta = cp.zeros((p, width), dtype=dtype)
        active = cp.ones(width, dtype=cp.bool_)
        previous_deviance = cp.full(width, cp.inf, dtype=dtype)
        for iteration in range(1, int(maxit) + 1):
            eta = xg @ beta
            mu, derivative = _link_values(eta, link, dtype)
            variance = cp.maximum(mu * (1.0 - mu), dtype(1e-12))
            working_weight = wg[:, None] * derivative * derivative / variance
            working_response = eta + (response - mu) / derivative
            gram = cp.einsum("ni,nj,nt->tij", xg, xg, working_weight)
            rhs = cp.einsum("ni,nt,nt->ti", xg, working_weight, working_response)
            try:
                proposed = cp.linalg.solve(gram, rhs[..., None])[..., 0].T
            except np.linalg.LinAlgError:
                proposed = cp.stack(
                    [cp.linalg.pinv(gram[j]) @ rhs[j] for j in range(width)],
                    axis=1,
                )
            proposed_eta = xg @ proposed
            proposed_mu, _ = _link_values(proposed_eta, link, dtype)
            deviance = -dtype(2.0) * cp.sum(
                wg[:, None] * (
                    response * cp.log(proposed_mu) +
                    (dtype(1.0) - response) * cp.log(dtype(1.0) - proposed_mu)
                ), axis=0
            )
            relative_deviance = cp.abs(deviance - previous_deviance) / (
                dtype(0.1) + cp.abs(deviance)
            )
            newly = active & (relative_deviance <= dtype(tolerance))
            if bool(cp.any(newly)):
                idx = cp.asnumpy(cp.where(newly)[0])
                iterations[first + idx] = iteration
            beta = cp.where(active[None, :], proposed, beta)
            previous_deviance = cp.where(active, deviance, previous_deviance)
            active = active & ~newly
            if not bool(cp.any(active)):
                break
        coefficients[:, first:last] = beta
        active_cpu = cp.asnumpy(active)
        converged[first:last] = ~active_cpu
        iterations[first:last][active_cpu] = int(maxit)
        final_eta = xg @ beta
        signed_margin = (dtype(2.0) * response - dtype(1.0)) * final_eta
        margin_tolerance = dtype(
            max(float(tolerance), 64.0 * float(np.finfo(dtype).eps))
        ) * cp.maximum(dtype(1.0), cp.max(cp.abs(final_eta), axis=0))
        boundary[first:last] = cp.asnumpy(
            cp.all(signed_margin >= -margin_tolerance[None, :], axis=0)
            & cp.any(signed_margin > margin_tolerance[None, :], axis=0)
        )
    cp.cuda.Stream.null.synchronize()
    return {
        "coefficients": cp.asnumpy(coefficients).astype(np.float64, copy=False),
        "iterations": iterations,
        "converged": converged,
        "boundary": boundary,
        "elapsed_seconds": time.perf_counter() - started,
    }


def fit_qr_admm(x, y, weights, taus, precision="float64", rho=1.0,
                maxit=5000, tolerance=1e-6, block_quantiles=8):
    """Experimental batched ADMM solver for weighted linear QR."""
    _require_cuda()
    dtype = _dtype(precision)
    xg = cp.asarray(x, dtype=dtype, order="C")
    yg = cp.asarray(y, dtype=dtype)
    wg = cp.asarray(weights, dtype=dtype)
    qg = cp.atleast_1d(cp.asarray(taus, dtype=dtype))
    n, p = xg.shape
    total = int(qg.size)
    coefficients = cp.empty((p, total), dtype=dtype)
    iterations = np.zeros(total, dtype=np.int32)
    converged = np.zeros(total, dtype=bool)
    primal = np.full(total, np.inf, dtype=np.float64)
    dual = np.full(total, np.inf, dtype=np.float64)
    gram = xg.T @ xg
    started = time.perf_counter()

    for first in range(0, total, int(block_quantiles)):
        last = min(first + int(block_quantiles), total)
        tau = qg[first:last]
        width = last - first
        residual = cp.zeros((n, width), dtype=dtype)
        scaled_dual = cp.zeros((n, width), dtype=dtype)
        beta = cp.zeros((p, width), dtype=dtype)
        for iteration in range(1, int(maxit) + 1):
            rhs = xg.T @ (yg[:, None] - residual + scaled_dual)
            beta = cp.linalg.solve(gram, rhs)
            linear_residual = yg[:, None] - xg @ beta
            value = linear_residual + scaled_dual
            previous = residual
            positive = wg[:, None] * tau[None, :] / dtype(rho)
            negative = wg[:, None] * (1.0 - tau[None, :]) / dtype(rho)
            residual = cp.where(
                value > positive,
                value - positive,
                cp.where(value < -negative, value + negative, dtype(0.0)),
            )
            primal_now = linear_residual - residual
            scaled_dual = scaled_dual + primal_now
            primal_norm = cp.max(cp.abs(primal_now), axis=0)
            dual_norm = dtype(rho) * cp.max(
                cp.abs(xg.T @ (residual - previous)), axis=0
            ) / dtype(max(n, 1))
            currently_converged = (primal_norm <= dtype(tolerance)) & (
                dual_norm <= dtype(tolerance)
            )
            if bool(cp.all(currently_converged)):
                break
        coefficients[:, first:last] = beta
        final_converged = (primal_norm <= dtype(tolerance)) & (
            dual_norm <= dtype(tolerance)
        )
        final_converged_cpu = cp.asnumpy(final_converged)
        converged[first:last] = final_converged_cpu
        iterations[first:last] = np.where(
            final_converged_cpu, iteration, int(maxit)
        )
        primal[first:last] = cp.asnumpy(primal_norm)
        dual[first:last] = cp.asnumpy(dual_norm)
    cp.cuda.Stream.null.synchronize()
    return {
        "coefficients": cp.asnumpy(coefficients).astype(np.float64, copy=False),
        "iterations": iterations,
        "converged": converged,
        "primal_residual": primal,
        "dual_residual": dual,
        "elapsed_seconds": time.perf_counter() - started,
    }
