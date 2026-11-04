"""Feed patterns and the geometric-optics (GO) aperture illumination.

This module turns a feed pattern into the amplitude distribution that lands on
the main reflector aperture. That distribution is the input to the physical
optics (PO) integral in ``po.py``.
"""

from __future__ import annotations

import numpy as np

from go_po.geometry import ParabolicDish


def cos_power_pattern(q: float):
    """Feed power pattern ``cos^q(theta)`` with q = 2 * (E-plane exponent).

    Returns a callable ``g(theta)`` giving the *voltage* pattern.
    """
    if q < 0:
        raise ValueError("q must be non-negative")

    def g(theta):
        theta = np.asarray(theta, dtype=float)
        vals = np.cos(np.clip(theta, 0.0, 0.5 * np.pi)) ** (0.5 * q)
        return np.where(np.abs(theta) <= 0.5 * np.pi, vals, 0.0)

    return g


def theta_of_rho(dish: ParabolicDish, rho):
    """Feed angle subtended by an aperture point at radius rho."""
    return 2.0 * np.arctan(np.asarray(rho, dtype=float) / (2.0 * dish.focal_length))


def slant_range(dish: ParabolicDish, theta):
    """Distance from the focus to the reflector surface (GO ray path)."""
    theta = np.asarray(theta, dtype=float)
    return 2.0 * dish.focal_length / (1.0 + np.cos(theta))


def aperture_illumination(
    dish: ParabolicDish,
    feed_pattern,
    *,
    normalize_edge_db: float | None = None,
):
    """Return a callable ``amp(rho)`` for the aperture voltage distribution.

    GO amplitude on the aperture follows

        A(rho) ~ g(theta) / r(theta)

    because the feed radiates a spherical wave (1/r spreading) and the
    reflector conserves the ray-tube power.

    Parameters
    ----------
    normalize_edge_db:
        If given, the distribution is scaled so that the aperture edge sits
        at this level relative to the centre (e.g. ``-10.0``).
    """

    def raw(rho):
        rho = np.asarray(rho, dtype=float)
        theta = theta_of_rho(dish, rho)
        return feed_pattern(theta) / slant_range(dish, theta)

    if normalize_edge_db is None:
        return raw

    centre = float(np.ravel(raw(np.array([0.0])))[0])
    # raw() 吃的是**口径半径 rho**，不是角度。原来这里传的是 dish.half_angle
    # （1.117 rad），等于把 -10 dB 套在 rho≈1.1λ 处，真正的边缘 rho=R 掉到 -14.5 dB。
    # 2026-11-04 修正为口径半径。
    edge = float(np.ravel(raw(np.array([dish.radius])))[0])
    target = 10.0 ** (normalize_edge_db / 20.0)

    def amp(rho):
        # affine rescale: raw == centre -> 1.0, raw == edge -> target
        a = raw(rho)
        return target + (a - edge) * (1.0 - target) / (centre - edge)

    return amp


def edge_taper_db(dish: ParabolicDish, feed_pattern) -> float:
    """Edge illumination relative to centre, in dB (negative means tapered)."""
    theta0 = dish.half_angle
    centre = feed_pattern(np.array([0.0])) / slant_range(dish, np.array([0.0]))
    edge = feed_pattern(np.array([theta0])) / slant_range(dish, np.array([theta0]))
    ratio = float(np.ravel(edge)[0]) / float(np.ravel(centre)[0])
    return 20.0 * np.log10(abs(ratio))
