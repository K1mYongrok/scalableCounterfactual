"""Optional CuPy kernels for scalableCounterfactual."""

from __future__ import annotations

import time
import numpy as np

try:
    import cupy as cp
except Exception as exc:  # reported by cuda_status
    cp = None
    _IMPORT_ERROR = repr(exc)
else:
    _IMPORT_ERROR = None


def cuda_status():
    if cp is None:
        return {"available": False, "error": _IMPORT_ERROR,
                "device": None, "cupy_version": None}
    try:
        if int(cp.cuda.runtime.getDeviceCount()) < 1:
            raise RuntimeError("no CUDA devices")
        props = cp.cuda.runtime.getDeviceProperties(0)
        name = props["name"]
        if isinstance(name, bytes):
            name = name.decode("utf-8")
        return {"available": True, "error": None, "device": name,
                "cupy_version": cp.__version__}
    except Exception as exc:
        return {"available": False, "error": repr(exc), "device": None,
                "cupy_version": getattr(cp, "__version__", None)}


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


def matmul(x, coefficients, precision="float64"):
    _require_cuda()
    dtype = _dtype(precision)
    result = cp.asarray(x, dtype=dtype, order="C") @ cp.asarray(
        coefficients, dtype=dtype, order="C"
    )
    cp.cuda.Stream.null.synchronize()
    return cp.asnumpy(result).astype(np.float64, copy=False)


def qr_marginal_quantiles(x, coefficients, weights, probs, shift=0.0,
                          precision="float64"):
    """Pool conditional QR draws and select weighted quantiles on the GPU.

    The rank positions and interpolation rule match Hmisc::wtd.quantile with
    type="quantile" and normwt=TRUE. Only the selected quantiles are copied
    back to the host.
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
    total_positions = int(n) * int(q)
    values = draws.ravel(order="F")
    ordering = cp.argsort(values)
    sorted_values = values[ordering]
    sorted_weights = wg[ordering % int(n)]
    cp.cumsum(sorted_weights, out=sorted_weights)

    rank = 1.0 + (float(total_positions) - 1.0) * pg
    low = cp.maximum(cp.floor(rank), 1.0)
    high = cp.minimum(low + 1.0, float(total_positions))
    interpolation = cp.mod(rank, 1.0)
    total_weight = cp.sum(wg) * float(q)
    targets = cp.concatenate((low, high)) / float(total_positions) * total_weight
    selected_index = cp.searchsorted(sorted_weights, targets, side="left")
    selected_index = cp.minimum(selected_index, total_positions - 1)
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
            probability = 1.0 / (1.0 + cp.exp(-cp.clip(eta, -709.0, 709.0)))
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
        mu = 1.0 / (1.0 + cp.exp(-cp.clip(eta, -709.0, 709.0)))
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
    eps = dtype(1e-12 if dtype == cp.float64 else 1e-6)
    return cp.clip(mu, eps, 1.0 - eps), cp.maximum(derivative, eps)


def fit_dr_process(x, y, weights, thresholds, link, precision="float64",
                   maxit=100, tolerance=1e-8, block_thresholds=16):
    """Batched unpenalized binomial IRLS with X resident on the GPU."""
    _require_cuda()
    dtype = _dtype(precision)
    xg = cp.asarray(x, dtype=dtype, order="C")
    yg = cp.asarray(y, dtype=dtype)
    wg = cp.asarray(weights, dtype=dtype)
    tg = cp.asarray(thresholds, dtype=dtype)
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
            except Exception:
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
        final_mu, _ = _link_values(xg @ beta, link, dtype)
        boundary[first:last] = cp.asnumpy(
            cp.any((final_mu <= dtype(1e-8)) |
                   (final_mu >= dtype(1.0 - 1e-8)), axis=0)
            | cp.any(cp.abs(beta) > dtype(100.0), axis=0)
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
