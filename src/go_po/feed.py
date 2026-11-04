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


def q_for_edge_taper(theta_edge: float, target_db: float) -> float:
    """给定边缘照射角与期望的边缘锥削，反解馈源的 ``cos^q`` 指数。

    口径照明是 ``A(theta) = g(theta) * cos^2(theta/2) / (M*F)``，双反射面这一路
    ``cos^2(theta/2)`` 是扩散因子。**注意 ``g`` 是电压方向图 ``cos^(q/2)``**
    （``cos_power_pattern(q)`` 返回的就是电压方向图，q 是功率指数），于是边缘锥削

        A(edge)/A(0) = cos^(q/2)(theta_edge) * cos^2(theta_edge/2)

    令它等于 ``10**(target_db/20)``，两边取自然对数：

        (q/2) * ln(cos(theta_edge)) + ln(cos^2(theta_edge/2)) = ln(10**(target_db/20))

    解出

        q = 2 * ln(10**(target_db/20) / cos^2(theta_edge/2)) / ln(cos(theta_edge))

    ⚠️ 分母那个 2 容易漏。漏了会让反解值正好小一半——2026-11-04 第一版就漏了，
    拿 q=8 正算再反解得到 4，一眼看出来。

    这是**设计倒推**用的工具：q 本身仍是喇叭的物理参数，不能反过来把 q 变成
    "由目标决定"的派生量。流水线里 q 依然是输入，这里只是给设计者一个换算。

    Parameters
    ----------
    theta_edge:
        馈源到副面边缘的照射半张角（弧度），必须在 (0, pi/2) 内。
        注意不是 ``Cassegrain.half_angle_feed``（那是到主面边缘的直线张角）。
    target_db:
        期望的边缘锥削（dB），通常为负。为正表示边缘比中心还亮，物理上做不到。

    Raises
    ------
    ValueError
        theta_edge 越界，或目标锥削为正（除以扩散因子后 > 1）。
    """
    th = float(theta_edge)
    if not 0.0 < th < 0.5 * np.pi:
        raise ValueError(f"theta_edge 必须在 (0, pi/2) 内，当前 {th:.6f} rad")

    spread = float(np.cos(0.5 * th) ** 2)
    target = 10.0 ** (float(target_db) / 20.0)
    ratio = target / spread
    if ratio > 1.0:
        raise ValueError(
            f"目标 {float(target_db):.3f} dB 在物理上做不到："
            f"锥削比 {ratio:.6f} > 1（扩散因子 cos²(θ/2) 已经只有 {spread:.6f}），"
            f"边缘不可能比中心亮"
        )
    if ratio <= 0.0:
        raise ValueError(f"锥削比必须为正，当前 {ratio:.3e}")

    # g(theta) = cos^(q/2)，不是 cos^q，所以解出来要乘 2
    return float(2.0 * np.log(ratio) / np.log(np.cos(th)))

