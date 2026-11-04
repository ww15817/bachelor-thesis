"""几何光学射线追踪的两个物理不变量。

这是本项目最重要的两个判据：
  1. 准直误差 —— 出射方向应平行于视轴
  2. 等光程   —— r1 + r2 - z_main 应为常数
只要双曲面取支、焦点位置或法向有任何错误，这两项会立刻失效。
"""

import numpy as np
import pytest


def test_all_rays_complete_the_chain(bundle):
    assert bundle.n_valid == bundle.n_rays, "存在未走完 馈源->副面->主面 的射线"


def test_collimation_error_is_at_machine_precision(bundle):
    err = bundle.collimation_error()
    assert err < 1e-12, f"准直误差过大: {err:.3e}"


def test_optical_path_is_constant(bundle):
    ripple = bundle.path_length_ripple()
    assert ripple < 1e-10, f"等光程起伏过大: {ripple:.3e}"


def test_outgoing_direction_is_along_axis(bundle):
    d = bundle.d_out[bundle.valid]
    assert np.allclose(d[:, 2], 1.0, atol=1e-12)


def test_edge_ray_lands_on_aperture_rim(geo):
    """边缘射线必须落到 rho = R（容差 1e-6 λ），而不是「差不多在口径里」。

    旧版断言 ``r.max() <= R*(1+1e-9)``：valid 掩码本身就按
    ``r <= R*(1+1e-9)`` 过滤，而且中点采样最外圈只到 9.577λ（R=10λ），
    根本没碰到边缘，断言必然成立。这里把射线角真取到 theta_max：
    ``trace_rays`` 在 theta_max/n_radial 的中点采样，n_radial=1 时
    唯一射线的角度是 theta_max/2，所以令 theta_max=2*theta_edge 即命中边缘。
    """
    from go_po.rays import feed_half_angle_for_edge, trace_rays

    theta_edge = feed_half_angle_for_edge(geo)
    bundle = trace_rays(geo, n_radial=1, n_azimuth=1, theta_max=2.0 * theta_edge)

    assert bundle.n_valid == 1
    assert abs(float(bundle.aperture_radius[0]) - geo.R) < 1e-6
    assert abs(float(bundle.p_main[0, 2]) - geo.z_rim_primary) < 1e-6


def test_theta_max_matches_edge_solver(geo):
    """默认 theta_max 必须就是 feed_half_angle_for_edge。

    采样角是 (i+0.5)/n_radial * theta_max，从包里的最大角反解 theta_max，
    与边缘求解器一致到 1e-12。
    """
    from go_po.rays import feed_half_angle_for_edge, trace_rays

    n_radial = 4
    bundle = trace_rays(geo, n_radial=n_radial, n_azimuth=1)
    recovered = float(bundle.theta.max()) * n_radial / (n_radial - 0.5)

    assert recovered == pytest.approx(feed_half_angle_for_edge(geo), rel=1e-12)


def test_ray_count_scales_with_sampling(geo):
    from go_po.rays import trace_rays

    b = trace_rays(geo, n_radial=6, n_azimuth=8)
    assert b.n_rays == 48
