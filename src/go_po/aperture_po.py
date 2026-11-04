"""物理光学 (PO) 口径积分：由口径场分布计算远场方向图。

原理
----
对电大尺寸反射面，最实用的做法是**等效口径法**：
把两次反射后的场等效成口径面上的一个面分布，再用物理光学辐射积分求远场。

设口径面位于 z = 0，口径场为 E_a(x, y)（标量近似，单一极化），
沿 z 轴观察方向为 (sinθ, 0, cosθ)，则归一化远场为

    F(θ) ∝ ∫∫ E_a(x, y) exp( j k x sinθ ) dx dy

即口径场的二维傅里叶变换在 k_x = k sinθ 处的取值。

**为什么要先验证它**：均匀照射的圆形口径有解析解
（Airy 图案，含 (1+cosθ)/2 斜度因子），可以用来校验积分器是否正确。
这是整个射频链路的第一层验证，也是论文里必须交代的一步。

验证目标（理论值）
------------------
* 第一零点： k R sinθ = 3.8317  →  sinθ = 0.6098 * λ / R
* 第一副瓣： 相对峰值 -17.6 dB（含斜度因子后约 -17.6 ~ -18 dB）
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["ApertureField", "uniform_aperture", "tapered_aperture", "far_field_pattern"]


@dataclass
class ApertureField:
    """口径面上的等效场分布。"""

    x: np.ndarray
    """采样点 x 坐标，形状 (N,)，单位 λ。"""

    y: np.ndarray
    """采样点 y 坐标，形状 (N,)，单位 λ。"""

    e: np.ndarray
    """复口径场，形状 (N,)。幅度 + 相位（相位含结构变形引入的误差）。"""

    radius: float
    """口径半径，单位 λ。"""

    weights: np.ndarray | None = None
    """积分权重；若为 None 则视为等面积均匀采样。"""

    @property
    def n_points(self) -> int:
        return self.x.size

    def power(self) -> float:
        """口径总功率（用于归一化）。"""
        w = self.weights if self.weights is not None else np.ones_like(self.e)
        return float(np.sum(w * np.abs(self.e) ** 2))


def uniform_aperture(radius: float, n_per_lambda: float = 12.0) -> ApertureField:
    """均匀照射的圆形口径（用于解析校验）。

    注意：为了积分精度必须用**面积元**而不是等间距网格，
    这里直接在圆内做极坐标均匀采样。
    """
    n_r = max(int(np.ceil(radius * n_per_lambda)), 8)
    n_phi = max(int(np.ceil(4.0 * n_r)), 24)

    r_edges = np.linspace(0.0, radius, n_r + 1)
    r_mid = 0.5 * (r_edges[1:] + r_edges[:-1])
    dr = r_edges[1] - r_edges[0]
    dphi = 2.0 * np.pi / n_phi

    RR, PP = np.meshgrid(r_mid, np.arange(n_phi) * dphi, indexing="ij")
    RR = RR.ravel()
    PP = PP.ravel()

    x = RR * np.cos(PP)
    y = RR * np.sin(PP)
    e = np.ones_like(x, dtype=complex)
    # 极坐标面积元 r dr dphi
    w = RR * dr * dphi

    return ApertureField(x=x, y=y, e=e, radius=radius, weights=w)


def tapered_aperture(
    radius: float,
    edge_taper_db: float = -10.0,
    phase_error: np.ndarray | callable | None = None,
    n_per_lambda: float = 12.0,
) -> ApertureField:
    """带轴对称幅度锥削（以及可选相位误差）的圆形口径。

    幅度模型：A(r) = 1 + (t - 1) * (r / R)^2，t = 10^(edge_taper_db / 20)。
    这是反射面馈源照射最常见的二次锥削近似。

    参数
    ----
    edge_taper_db : float
        边缘照射电平（dB，负值）。
    phase_error : ndarray | callable | None
        相位误差。可传半径 r 的函数（弧度），或与采样点等长的数组。
        结构变形引起的口径相位扰动通常写成 dphi = 2k w(r)。
    """
    af = uniform_aperture(radius, n_per_lambda=n_per_lambda)
    r = np.hypot(af.x, af.y)

    t = 10.0 ** (edge_taper_db / 20.0)
    amp = 1.0 + (t - 1.0) * (r / radius) ** 2
    af.e = amp.astype(complex)

    if phase_error is not None:
        if callable(phase_error):
            dphi = phase_error(r)
        else:
            dphi = np.asarray(phase_error, dtype=float)
            if dphi.shape != r.shape:
                raise ValueError(
                    f"phase_error 形状 {dphi.shape} 与采样点 {r.shape} 不符"
                )
        af.e = af.e * np.exp(1j * dphi)

    return af


def far_field_pattern(
    af: ApertureField,
    theta_deg: np.ndarray,
    obliquity: bool = True,
    *,
    normalize: float | None = None,
) -> np.ndarray:
    """由口径场计算 E 面（phi = 0）远场方向图。

    参数
    ----
    theta_deg : ndarray
        观察角（度）。
    obliquity : bool
        是否乘上口径天线的斜度因子 (1 + cosθ) / 2。
    normalize : float | None
        显式归一化除数（绝对值）。

        - ``None``（默认）：按本次结果的峰值归一化，**增益信息会被抹掉**——
          拿这个比较「有变形 / 无变形」永远得到 0 dB。
        - 传一个数：用它来除，例如传未变形口径的轴向幅度，
          这样 A/B 两组结果共用同一参考，增益损失才出得来。

        2026-11-03 审查第 6 条：原来这里只有峰值归一化一条路，
        ``tapered_aperture(phase_error=...)`` 又是唯一带复相位的口径工厂，
        于是「变形 -> 增益下降」这条路走不通。
    """
    k = 2.0 * np.pi  # 长度以 λ 归一化，故 k = 2π
    th = np.deg2rad(theta_deg)
    kx = k * np.sin(th)

    w = af.weights if af.weights is not None else np.ones(af.n_points)

    # 逐 theta 循环，不一次性构造 (n_theta, n_points) 的复数组。
    # 2026-11-04 审查 C-01：旧写法在 validate_po.py 默认参数下要 4.4 GB，
    # 整个脚本峰值 10 GB，8 GB 笔记本直接挂。数值逐位不变（同一个求和，分批做）。
    we = w * af.e
    field = np.empty(th.size, dtype=complex)
    for i in range(th.size):
        field[i] = np.sum(np.exp(1j * kx[i] * af.x) * we)

    if obliquity:
        field = field * (1.0 + np.cos(th)) / 2.0

    if normalize is None:
        peak = np.max(np.abs(field))
        if peak <= 0.0:
            raise ValueError("远场峰值为零，口径场可能全为零")
        return np.abs(field) / peak

    if abs(normalize) <= 0.0:
        raise ValueError("normalize 不能为 0")
    return np.abs(field) / abs(normalize)


def axial_amplitude(af: ApertureField, obliquity: bool = True) -> float:
    """口径场的轴向（theta=0）远场幅度，用作 ``normalize`` 的参考值。

    口径场带相位扰动时，这个值会随扰动下降——这正是要量的增益损失。
    """
    return float(np.abs(far_field_pattern(af, np.array([0.0]), obliquity=obliquity, normalize=1.0))[0])



def analytic_uniform_circular(radius: float, theta_deg: np.ndarray) -> np.ndarray:
    """均匀圆形口径的解析方向图（Airy 图案 + 斜度因子），用于校验。"""
    from scipy.special import j1

    k = 2.0 * np.pi
    th = np.deg2rad(theta_deg)
    u = k * radius * np.sin(th)

    with np.errstate(divide="ignore", invalid="ignore"):
        core = np.where(np.abs(u) < 1e-12, 1.0, 2.0 * j1(u) / np.where(u == 0, 1.0, u))
    pat = np.abs(core * (1.0 + np.cos(th)) / 2.0)
    return pat / np.max(pat)


if __name__ == "__main__":
    R = 10.0  # 20λ 口径
    theta = np.linspace(-12.0, 12.0, 2401)

    print("=== 校验：均匀圆口径 PO 积分 vs 解析 Airy ===")
    af = uniform_aperture(R, n_per_lambda=14.0)
    print(f"采样点数 = {af.n_points}")

    num = far_field_pattern(af, theta)
    ana = analytic_uniform_circular(R, theta)

    err = np.max(np.abs(num - ana))
    print(f"最大绝对偏差 = {err:.3e}")

    # 第一零点位置
    def first_null(pat, th):
        i0 = np.argmax(pat)
        for i in range(i0 + 1, len(pat) - 1):
            if pat[i + 1] > pat[i]:
                return th[i]
        return np.nan

    null_num = first_null(num, theta)
    null_ana = first_null(ana, theta)
    theory = np.rad2deg(np.arcsin(3.8317 / (2 * np.pi * R)))
    print(f"第一零点: 数值={null_num:.4f}°  解析={null_ana:.4f}°  理论={theory:.4f}°")

    # 第一副瓣
    def first_sidelobe_db(pat, th):
        i0 = np.argmax(pat)
        # 找第一个极小
        j = i0
        while j < len(pat) - 1 and pat[j + 1] <= pat[j]:
            j += 1
        # 再找下一个极大
        k2 = j
        while k2 < len(pat) - 1 and pat[k2 + 1] >= pat[k2]:
            k2 += 1
        return 20 * np.log10(pat[k2]), th[k2]

    sll_num, th_num = first_sidelobe_db(num, theta)
    sll_ana, th_ana = first_sidelobe_db(ana, theta)
    print(f"第一副瓣: 数值={sll_num:.2f} dB @ {th_num:.3f}°   解析={sll_ana:.2f} dB @ {th_ana:.3f}°")
    print()
    print("=== 锥削影响（-10 dB 边缘照射）===")
    for taper in (0.0, -10.0, -15.0):
        a = tapered_aperture(R, edge_taper_db=taper, n_per_lambda=12.0)
        p = far_field_pattern(a, theta)
        s, t = first_sidelobe_db(p, theta)
        n = first_null(p, theta)
        print(f"  边缘 {taper:6.1f} dB -> 第一零点 {n:7.4f}°  第一副瓣 {s:7.2f} dB")
