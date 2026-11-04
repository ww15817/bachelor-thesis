"""馈源方向图与 GO 口径照明的解析校核。

这些函数是 po.py 口径积分的输入，之前零覆盖。判据全部取自解析几何/解析积分，
不是「函数返回它自己定义的值」：

* 抛物面的焦点张角满足 tan(theta) = rho / (F - z)；
* 焦点到反射面的距离就是 GO 光程 2F/(1+cos theta)；
* 射线管功率守恒给出 ∫|A|^2 dA = ∫|g|^2 sin(theta) dtheta；
* 1/r 扩散 = cos^2(theta/2)/F，由此得到边缘锥削的解析式。
"""

import numpy as np
import pytest

from go_po.feed import (
    aperture_illumination,
    cos_power_pattern,
    edge_taper_db,
    slant_range,
    theta_of_rho,
)
from go_po.geometry import ParabolicDish, build_polar_grid

DISH = ParabolicDish(diameter=20.0, f_over_d=0.4)


def test_cos_power_pattern_boresight_is_one():
    """theta=0 处馈源电压方向图必须为 1，且与 q 无关。"""
    for q in (0.0, 2.0, 4.0, 8.0):
        g = cos_power_pattern(q)
        assert abs(float(g(0.0)) - 1.0) < 1e-15


def test_cos_power_pattern_half_power_angle():
    """半功率角 theta_hp = arccos(2^(-1/q)) 处，电压方向图 = 1/sqrt(2)。

    函数返回的是**电压**方向图 cos^(q/2)；如果误用 cos^q，这里会得到 1/2。
    """
    for q in (1.0, 2.0, 4.0, 8.0):
        theta_hp = np.arccos(2.0 ** (-1.0 / q))
        g = cos_power_pattern(q)
        assert float(g(theta_hp)) == pytest.approx(2.0 ** -0.5, rel=1e-12)


def test_cos_power_pattern_decreases_monotonically():
    """q>0 时电压方向图从 1 严格单调降到 0；q=0 时在 [0, pi/2] 上恒为 1。"""
    theta = np.linspace(0.0, 0.5 * np.pi - 1e-6, 401)
    for q in (2.0, 8.0):
        vals = cos_power_pattern(q)(theta)
        assert np.all(np.diff(vals) < 0.0)
        assert float(vals[0]) == pytest.approx(1.0, abs=1e-15)
        assert float(cos_power_pattern(q)(0.5 * np.pi)) < 1e-15

    uniform = cos_power_pattern(0.0)(theta)
    assert np.allclose(uniform, 1.0, atol=1e-15)


def test_slant_range_is_the_focus_distance():
    """slant_range(theta(rho)) 必须等于焦点到抛物面上该点的欧氏距离。"""
    focus = np.array([0.0, 0.0, DISH.focal_length])
    for rho in (0.0, 2.0, 5.0, 9.0, DISH.radius):
        theta = theta_of_rho(DISH, rho)
        point = DISH.surface_point(np.array([rho]), np.array([0.0]))[0]
        geometric = float(np.linalg.norm(point - focus))
        assert abs(float(slant_range(DISH, theta)) - geometric) < 1e-12


def test_theta_of_rho_matches_focus_geometry():
    """馈源张角由焦点几何决定： tan(theta) = rho / (F - z)，且 theta(R) = 半张角。"""
    for rho in (0.0, 2.0, 5.0, 9.0, DISH.radius):
        theta = float(theta_of_rho(DISH, rho))
        z = float(DISH.z_of_rho(rho))
        expected = float(np.arctan(rho / (DISH.focal_length - z)))
        assert theta == pytest.approx(expected, rel=1e-12, abs=1e-15)

    assert float(theta_of_rho(DISH, DISH.radius)) == pytest.approx(
        DISH.half_angle, rel=1e-12
    )


def test_aperture_illumination_matches_go_amplitude():
    """GO 口径幅度 A = g(theta)/r(theta) = g(theta) cos^2(theta/2) / F。"""
    q = 8.0
    feed = cos_power_pattern(q)
    amp = aperture_illumination(DISH, feed)

    rho = np.array([0.0, 2.0, 5.0, 8.0, DISH.radius])
    theta = theta_of_rho(DISH, rho)
    expected = feed(theta) * np.cos(theta / 2.0) ** 2 / DISH.focal_length
    assert np.allclose(amp(rho), expected, rtol=1e-12, atol=1e-15)


def test_aperture_illumination_conserves_feed_power():
    """射线管功率守恒：∫|A|^2 dA 必须等于馈源总功率 ∫|g|^2 sin(theta) dtheta dphi。

    抛物面的射线管 Jacobian 是 rho d rho = r^2 sin(theta) dtheta，所以 A=g/r
    恰好让口径功率等于馈源功率。1/r 因子、theta(rho) 映射或面元 dA 任一写错，
    这个积分都会偏离解析值。中点法的离散误差是 O(drho^2)，实测 3.7e-6。
    """
    q = 8.0
    feed = cos_power_pattern(q)
    amp = aperture_illumination(DISH, feed)

    rho, _, dA = build_polar_grid(DISH.radius, 400, 360)
    numeric = float(np.sum(np.abs(amp(rho)) ** 2 * dA))

    theta0 = DISH.half_angle
    analytic = 2.0 * np.pi * (1.0 - np.cos(theta0) ** (q + 1.0)) / (q + 1.0)
    assert abs(numeric - analytic) / analytic < 1e-5


def test_edge_taper_matches_go_spreading():
    """边缘锥削的解析值 = 10q log10(cos theta0) + 40 log10(cos(theta0/2))。

    前一项是馈源方向图，后一项是 1/r 扩散（1/r = cos^2(theta/2)/F）。
    若把电压当功率（少一个 2 倍）或漏掉扩散因子，这里立刻对不上。
    """
    theta0 = DISH.half_angle
    for q in (0.0, 2.0, 8.0):
        feed = cos_power_pattern(q)
        expected = 10.0 * q * np.log10(np.cos(theta0)) + 40.0 * np.log10(
            np.cos(theta0 / 2.0)
        )
        assert edge_taper_db(DISH, feed) == pytest.approx(expected, abs=1e-9)


def test_aperture_illumination_edge_normalisation():
    """normalize_edge_db 必须把口径边缘 rho=R 精确放到目标电平。

    旧实现把 ``dish.half_angle``（1.117 rad）当成 rho 传给 raw()，
    结果目标电平落在 rho≈1.1λ，真正的边缘掉到 +23 dB。这里用恒等式
    amp(R)/amp(0) = 10^(target/20) 来卡；中间点还要跟解析 GO 幅度的
    仿射插值一致（排除「在 dB 域线性插值」这类只在端点对的实现）。
    """
    q = 8.0
    feed = cos_power_pattern(q)
    f = DISH.focal_length

    def analytic_amp(rho):
        theta = theta_of_rho(DISH, rho)
        return feed(theta) * np.cos(theta / 2.0) ** 2 / f

    for target_db in (-3.0, -10.0, -20.0):
        amp = aperture_illumination(DISH, feed, normalize_edge_db=target_db)
        target = 10.0 ** (target_db / 20.0)

        centre = float(amp(0.0))
        edge = float(amp(DISH.radius))
        assert centre == pytest.approx(1.0, abs=1e-12)
        assert edge / centre == pytest.approx(target, rel=1e-12)

        rho_mid = 0.5 * DISH.radius
        expected_mid = target + (analytic_amp(rho_mid) - analytic_amp(DISH.radius)) * (
            1.0 - target
        ) / (analytic_amp(0.0) - analytic_amp(DISH.radius))
        assert float(amp(rho_mid)) == pytest.approx(expected_mid, rel=1e-12)

