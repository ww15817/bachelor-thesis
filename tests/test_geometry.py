"""卡塞格伦几何的物理自检。

几何量之间有很多「互相推导出来」的恒等式（z_focus_near = z_center + c、
M = (c+a)/(c-a) 之类），把定义代回去必然成立，测不到建模是否正确。
这里改成验卡塞格伦的物理定义：射线反向过近焦点、边缘三点共线、
M 作为系统焦距的物理含义、双曲线的焦点—准线定义。
"""

import numpy as np
import pytest

from common.cassegrain import Cassegrain


def _edge_ray(geo):
    """返回恰好打到主面边缘的那一条射线。

    ``trace_rays`` 取 (i+0.5)/n_radial * theta_max 的中点采样；
    n_radial=1 且 theta_max=2*theta_edge 时，唯一射线的角度正好是 theta_edge。
    """
    from go_po.rays import feed_half_angle_for_edge, trace_rays

    theta_edge = feed_half_angle_for_edge(geo)
    bundle = trace_rays(geo, n_radial=1, n_azimuth=1, theta_max=2.0 * theta_edge)
    assert bundle.n_valid == 1
    return bundle


def test_near_focus_coincides_with_parabola_focus(geo, bundle):
    """副面反射线的反向延长必须穿过 z = F（抛物面焦点 = 双曲面近焦点）。

    旧版断言 ``z_focus_near == geo.F``，但 ``z_focus_near`` 定义为
    ``z_center + c``，而 ``c - a = F - L`` 也是定义，代回去恒等于 F——
    a、c 取任何值都能过。这里用射线追踪的物理结果：每一条
    「副面 -> 主面」射线所在的直线都必须经过 (0, 0, F)。
    """
    assert abs(geo.z_focus_near - geo.F) < 1e-12

    focus = np.array([0.0, 0.0, geo.F])
    p_sub = bundle.p_sub[bundle.valid]
    p_main = bundle.p_main[bundle.valid]

    direction = p_main - p_sub
    direction /= np.linalg.norm(direction, axis=-1, keepdims=True)
    offset = focus[None, :] - p_main
    distance = np.linalg.norm(np.cross(offset, direction), axis=-1)

    assert float(distance.max()) < 1e-12, (
        f"反射线未过近焦点 (0,0,{geo.F})，最大垂距 {distance.max():.3e} λ"
    )


def test_sub_rim_main_rim_and_near_focus_are_collinear(geo):
    """副面边缘、主面边缘与 F 必须三点共线。

    这是卡塞格伦的定义（副面反射线反向过近焦点，近焦点 = 主面焦点）。
    同时用追踪到的原射线去校核 ``sub_rim()`` 的独立二分结果。
    """
    ray = _edge_ray(geo)
    p_sub = ray.p_sub[0]
    p_main = ray.p_main[0]

    # 追踪的落点确实在主面边缘
    assert np.hypot(p_main[0], p_main[1]) == pytest.approx(geo.R, abs=1e-6)
    assert p_main[2] == pytest.approx(geo.z_rim_primary, abs=1e-6)

    # 追踪到的副面落点 = 几何模块独立二分出的副面边缘
    rim_r, rim_z = geo.sub_rim()
    assert np.hypot(p_sub[0], p_sub[1]) == pytest.approx(rim_r, abs=1e-6)
    assert p_sub[2] == pytest.approx(rim_z, abs=1e-6)

    # F、副面边缘、主面边缘共线
    focus = np.array([0.0, 0.0, geo.F])
    direction = p_main - p_sub
    direction /= np.linalg.norm(direction)
    distance = float(np.linalg.norm(np.cross(focus - p_main, direction)))
    assert distance < 1e-12


def test_magnification_is_system_focal_length_ratio(geo):
    """M 必须满足 rho = 2 (M F) tan(theta/2)。

    卡塞格伦等效于一个焦距 ``f_sys = M F`` 的抛物面（20λ 算例为 32λ）：
    馈源张角 theta 的射线落到口径半径 rho。这是 M 的物理含义，也是流水线
    构造口径幅度时用的 dA/dOmega。旧版断言 ``(c+a)/(c-a) == M``，
    只是把 a、c 的定义代回去。
    """
    from go_po.rays import _trace_one_direction

    theta = np.radians(np.array([3.0, 7.0, 11.0, 15.0]))
    _, p_main, _, _, _, ok = _trace_one_direction(geo, theta, np.zeros_like(theta))
    assert bool(ok.all())

    rho = np.hypot(p_main[:, 0], p_main[:, 1])
    expected = 2.0 * geo.M * geo.F * np.tan(theta / 2.0)
    assert np.allclose(rho, expected, rtol=1e-12, atol=1e-12)


def test_eccentricity_from_focus_directrix_definition(geo):
    """双曲线右支上任意点到近焦点与到准线的距离比必须等于 e。

    准线在 z = z_center + a/e。点由 ``rim_radius_at()``（依赖 b²）给出，
    所以 b²、a、c 必须同时正确；这不是把 ``c/a`` 的定义再代一遍。
    """
    focus = np.array([0.0, 0.0, geo.F])
    directrix_z = geo.z_center + geo.a / geo.e

    z = np.linspace(geo.z_center + geo.a + 0.05, geo.sub_rim()[1], 9)
    r = np.array([geo.rim_radius_at(float(zz)) for zz in z])
    points = np.stack([r, np.zeros_like(r), z], axis=1)

    near = np.linalg.norm(points - focus[None, :], axis=1)
    to_directrix = np.abs(z - directrix_z)
    assert np.allclose(near / to_directrix, geo.e, rtol=1e-12)


def test_consistency_report_all_pass(geo):
    report = geo.consistency_report()
    assert report, "自检报告不应为空"
    failed = [name for name, _, ok in report if not ok]
    assert not failed, f"以下自检项未通过: {failed}"


def test_blockage_is_reasonable(geo):
    """副面遮挡比必须为正且小于 15%（工程上可接受的上限）。"""
    assert 0.0 < geo.blocking_ratio < 0.15


def test_sub_rim_is_computed_not_zero(geo):
    """副面边缘不能退化为 0——早期版本曾在顶点处求半径，结果恒为 0。"""
    rim = geo.sub_rim()
    assert rim is not None
    r, z = rim
    assert r > 0.5
    assert z > geo.L


def test_invalid_L_rejected():
    with pytest.raises(ValueError):
        Cassegrain(D=20.0, F=8.0, M=4.0, L=9.0)   # L 必须 < F


def test_invalid_M_rejected():
    with pytest.raises(ValueError):
        Cassegrain(D=20.0, F=8.0, M=0.8)
