"""物理光学积分器的解析校核。

均匀照射的圆口径有解析解（艾里图案 2J1(u)/u）。
数值积分必须与它吻合，否则方向图的副瓣不可信。
"""

import numpy as np
import pytest

from go_po.feed import aperture_illumination, cos_power_pattern
from go_po.geometry import ParabolicDish
from go_po.po import (
    analytic_uniform_circular_aperture,
    far_field_aperture,
    po_surface_far_field,
)


def _uniform(r):
    return np.ones_like(r)


def test_uniform_aperture_matches_airy(geo):
    k = 2.0 * np.pi
    theta = np.linspace(0.0, np.radians(12.0), 201)

    numeric = far_field_aperture(geo.R, k, _uniform, theta, n_rho=160, n_phi=288)
    exact = analytic_uniform_circular_aperture(geo.R, k, theta)

    assert np.max(np.abs(numeric - exact)) < 5e-3


def test_first_null_position(geo):
    """第一零点应出现在 k*a*sin(theta) = 3.8317。"""
    k = 2.0 * np.pi
    theta = np.linspace(0.0, np.radians(12.0), 1601)
    numeric = np.abs(far_field_aperture(geo.R, k, _uniform, theta, n_rho=200, n_phi=320))

    # 注意：不能用 argmin 找全局最小——艾里图案的零点越往后越深，
    # 全局最小落在 9.3° 附近。第一零点要取「第一个局部极小」。
    first = next(
        i for i in range(1, len(numeric) - 1)
        if numeric[i] < numeric[i - 1] and numeric[i] < numeric[i + 1]
    )
    expected = np.arcsin(3.8317 / (k * geo.R))
    assert theta[first] == pytest.approx(expected, abs=2e-3)


def test_boresight_equals_aperture_integral(geo):
    """视轴场必须等于口径照明的面积分（Kirchhoff 轴向极限）。

    旧版用默认归一化断言 |E(0)| = 1：此时分子分母都是 Σamp dA，恒等于 1，
    任何实现都过。这里关掉归一化（normalize=1.0），直接对解析积分：

        E(0) = ∫ (1 - rho^2/(2R^2)) dA = 3πR²/4

    中点法误差 O(drho²)，n_rho=220 时实测相对误差 3.4e-6（容差 1e-5）。
    漏乘 dA、dA 写成 rho*drho（少 dphi）或网格半径错位都会破坏这个等式。
    """
    k = 2.0 * np.pi
    theta = np.array([0.0])

    def tapered(r):
        return 1.0 - 0.5 * (np.asarray(r, dtype=float) / geo.R) ** 2

    val = far_field_aperture(geo.R, k, tapered, theta, n_rho=220, n_phi=360,
                             normalize=1.0)
    exact = 3.0 * np.pi * geo.R**2 / 4.0
    assert abs(val[0] - exact) / exact < 1e-5


def test_explicit_normalisation_scales_result(geo):
    """normalize 是绝对除数，默认值是口径面元之和；
    显式传 2 倍参考时，结果应恰好减半。"""
    from go_po.geometry import build_polar_grid

    k = 2.0 * np.pi
    theta = np.array([0.0])

    # 与 far_field_aperture 的默认采样保持一致，才能复现它的内部参考值
    _, _, dA = build_polar_grid(geo.R, 220, 360)
    ref = float(dA.sum())

    a = far_field_aperture(geo.R, k, _uniform, theta, normalize=ref)
    b = far_field_aperture(geo.R, k, _uniform, theta, normalize=2.0 * ref)
    assert abs(b[0] - a[0] / 2.0) < 1e-9


def test_complex_illumination_preserves_power(geo):
    """纯相位照明不改变总辐射功率（帕塞瓦尔定理）。

    2026-11-03 审查第 2 条：旧版把 illumination 的返回值截成实部，
    ``cos(phi)`` 相当于给口径加了幅度权重，总辐射功率会掉一大截。
    这条判据比"峰值是否离开视轴"有力：取实部之后波束照样会偏，
    但口径功率不可能守恒。
    """
    k = 2.0 * np.pi
    A = 5.0
    theta = np.radians(np.linspace(0.0, 80.0, 801))
    weight = np.sin(theta)

    def flat(r):
        return np.ones_like(r, dtype=complex)

    def phase_only(r):
        return np.exp(1j * A * (r / geo.R) ** 2)

    norm0 = far_field_aperture(geo.R, k, flat, np.array([0.0]), normalize=1.0)[0]
    u = np.abs(
        far_field_aperture(geo.R, k, flat, theta, normalize=norm0, n_rho=120, n_phi=240)
    )
    p = np.abs(
        far_field_aperture(geo.R, k, phase_only, theta, normalize=norm0, n_rho=120, n_phi=240)
    )

    ratio = float(np.sum(p ** 2 * weight) / np.sum(u ** 2 * weight))
    assert 0.9 < ratio < 1.1, f"纯相位照明改变了总功率（比值 {ratio:.3f}），说明复场被截成了实部"


# ---------------------------------------------------------------------------
# 显式面电流 PO（po_surface_far_field）——之前零覆盖
# ---------------------------------------------------------------------------

SURFACE_DISH = ParabolicDish(diameter=20.0, f_over_d=0.4)


def test_surface_po_matches_aperture_po_in_main_lobe():
    """面电流 PO 与口径 Kirchhoff PO 在主瓣区必须一致（< 1 dB）。

    两条通路是对同一个物理光学远场的不同写法：一个在反射面上积 (n×H)e^{jkr·r'}dS，
    一个在投影口径上积照明 e^{jkr·rho}dA。数值上它们只差离散误差和口径通路
    缺的斜度因子 (1+cos theta)/2（8° 处约 0.04 dB）。第一零点在 10.1°，
    取到 8°（方向图已降到 -42 dB）仍属主瓣；实测最大偏差 0.52 dB。
    """
    k = 2.0 * np.pi
    feed = cos_power_pattern(8.0)
    illumination = aperture_illumination(SURFACE_DISH, feed)
    theta = np.radians(np.linspace(0.0, 8.0, 81))

    aperture = far_field_aperture(
        SURFACE_DISH.radius, k, illumination, theta, n_rho=140, n_phi=240
    )
    surface = po_surface_far_field(
        SURFACE_DISH, k, feed, theta, n_rho=140, n_phi=240
    )
    surface = np.linalg.norm(surface, axis=1)

    ap_db = 20.0 * np.log10(np.abs(aperture) / np.abs(aperture[0]))
    surf_db = 20.0 * np.log10(surface / surface[0])
    assert float(np.max(np.abs(surf_db - ap_db))) < 1.0


def test_surface_po_normalisation_is_call_independent():
    """归一化参考必须是轴向幅度，不能是「传入数组的第一个角度」。

    theta 数组不同（一个含 0、一个不含）时，同一个角度的返回值必须相同。
    旧实现用 e_far[0] 归一化：theta0=0 与 theta0=6° 在同一个 5° 点相差 24 倍。
    这里两次调用在相同 theta 上的求和完全相同，只允许浮点等价（1e-12）。
    """
    k = 2.0 * np.pi
    feed = cos_power_pattern(8.0)

    with_zero = po_surface_far_field(
        SURFACE_DISH, k, feed, np.array([0.0, 3.0, 5.0]), n_rho=100, n_phi=160
    )
    without_zero = po_surface_far_field(
        SURFACE_DISH, k, feed, np.array([3.0, 5.0]), n_rho=100, n_phi=160
    )

    assert abs(float(np.linalg.norm(with_zero[0])) - 1.0) < 1e-12
    assert np.allclose(with_zero[1:], without_zero, rtol=0.0, atol=1e-12)


