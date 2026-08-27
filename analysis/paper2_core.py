"""Small, self-contained numerical helpers used by the Paper 2 analyses.

These functions are copied from the frozen analysis implementation so that
the release does not depend on unrelated workstation modules.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.stats import weibull_min


ALPHA_CAP = 1e4


def xy_from_latlon(lat: np.ndarray, lon: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Convert latitude/longitude to local equirectangular coordinates (m)."""
    lat_rad = np.radians(np.asarray(lat, dtype=float))
    lon_rad = np.radians(np.asarray(lon, dtype=float))
    lat0 = np.nanmean(lat_rad)
    lon0 = np.nanmean(lon_rad)
    radius_m = 6_371_000.0
    x = radius_m * np.cos(lat0) * (lon_rad - lon0)
    y = radius_m * (lat_rad - lat0)
    return x, y


@dataclass(frozen=True)
class NDEvent:
    step_id: int
    start_idx: int
    endpoint_idx: int
    trigger_idx: int | None
    is_terminal: bool
    step_length: float
    triggered_by_segment_crossing: bool
    n_points_event_path: int


def point_segment_min_dist_nd(point: np.ndarray, a: np.ndarray, b: np.ndarray) -> float:
    """Euclidean distance from a point to a closed line segment."""
    ab = b - a
    denominator = float(np.dot(ab, ab))
    if denominator <= 0:
        return float(np.linalg.norm(point - a))
    fraction = float(np.dot(point - a, ab) / denominator)
    fraction = min(1.0, max(0.0, fraction))
    return float(np.linalg.norm(point - (a + fraction * ab)))


def extract_continuous_drawdown_events_nd(points: np.ndarray, delta: float) -> list[NDEvent]:
    """Dimension-general continuous-path radial-drawdown detector."""
    pts = np.asarray(points, dtype=float)
    n_points = len(pts)
    events: list[NDEvent] = []
    if n_points < 3:
        return events
    start = 0
    step_id = 0
    while start < n_points - 2:
        center = pts[start]
        rmax = 0.0
        peak_idx = start
        trigger_idx: int | None = None
        crossing = False
        for k in range(start + 1, n_points):
            if rmax >= delta:
                dmin = point_segment_min_dist_nd(center, pts[k - 1], pts[k])
                if dmin <= (rmax - delta):
                    trigger_idx = k
                    crossing = True
                    break
            radius = float(np.linalg.norm(pts[k] - center))
            if radius > rmax:
                rmax = radius
                peak_idx = k
            if rmax >= delta and (rmax - radius) >= delta:
                trigger_idx = k
                crossing = False
                break
        if trigger_idx is None:
            endpoint = n_points - 1
            length = float(np.linalg.norm(pts[endpoint] - center))
            if length > 0:
                events.append(
                    NDEvent(step_id, start, endpoint, None, True, length, False,
                            endpoint - start + 1)
                )
            break
        endpoint = peak_idx
        length = float(np.linalg.norm(pts[endpoint] - center))
        if endpoint <= start or length <= 0:
            start += 1
            continue
        events.append(
            NDEvent(step_id, start, endpoint, trigger_idx, False, length,
                    crossing, trigger_idx - start + 1)
        )
        step_id += 1
        start = endpoint
    return events


def fit_step_distributions(lengths: np.ndarray) -> dict[str, float | int | bool]:
    """Fit the frozen exponential, Lomax, lognormal and Weibull candidates."""
    x = np.asarray(lengths, float)
    x = x[np.isfinite(x) & (x > 0)]
    if len(x) == 0:
        raise ValueError("at least one positive finite step length is required")
    n = len(x)
    lam = 1.0 / float(x.mean())
    ll_exp = float(n * math.log(lam) - lam * x.sum())

    def lomax_nll(theta: np.ndarray) -> float:
        alpha = math.exp(theta[0])
        scale = math.exp(theta[1])
        return float(-(
            n * math.log(alpha) - n * math.log(scale)
            - (alpha + 1.0) * np.log1p(x / scale).sum()
        ))

    lomax_fit = minimize(
        lomax_nll,
        np.array([0.0, math.log(max(float(np.median(x)), 1e-6))]),
        method="L-BFGS-B",
        bounds=[
            (-6, math.log(ALPHA_CAP)),
            (math.log(1e-6), math.log(max(float(x.max()) * 100, 1e-4))),
        ],
        options={"maxiter": 500},
    )
    alpha = float(math.exp(lomax_fit.x[0]))
    scale = float(math.exp(lomax_fit.x[1]))
    ll_lomax = float(-lomax_fit.fun)
    log_x = np.log(x)
    mu = float(log_x.mean())
    sigma = max(float(log_x.std(ddof=0)), 1e-8)
    ll_lognorm = float(np.sum(
        -np.log(x * sigma * math.sqrt(2 * math.pi))
        - 0.5 * ((log_x - mu) / sigma) ** 2
    ))
    try:
        shape_w, _, scale_w = weibull_min.fit(x, floc=0)
        ll_weibull = float(np.sum(weibull_min.logpdf(x, shape_w, loc=0, scale=scale_w)))
    except Exception:
        shape_w, scale_w, ll_weibull = float("nan"), float("nan"), float("nan")
    return {
        "n_steps": n,
        "exp_rate": lam,
        "ll_exp": ll_exp,
        "aic_exp": 2 - 2 * ll_exp,
        "alpha_step": alpha,
        "mu_step": alpha + 1.0,
        "lomax_scale": scale,
        "ll_lomax": ll_lomax,
        "aic_lomax": 4 - 2 * ll_lomax,
        "lomax_success": bool(lomax_fit.success),
        "lognorm_mu": mu,
        "lognorm_sigma": sigma,
        "ll_lognorm": ll_lognorm,
        "aic_lognorm": 4 - 2 * ll_lognorm,
        "weibull_shape": float(shape_w),
        "weibull_scale": float(scale_w),
        "ll_weibull": ll_weibull,
        "aic_weibull": float(4 - 2 * ll_weibull) if np.isfinite(ll_weibull) else float("nan"),
    }
