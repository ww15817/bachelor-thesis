"""Milestone M3: 20λ 卡塞格伦天线的 GO/PO 双反射面链路。

链路
----
    馈源(双曲面远焦点) --副面--> --主面--> 平行出射
        -> 口径场(含副面遮挡) -> 物理光学远场

物理自校验
----------
* 准直误差：出射方向应精确平行于 z 轴（数值上 ~1e-15）
* 光程起伏：r1 + r2 - z_main 应为常数（数值上 ~1e-14）
* 几何自检：近焦点必须精确落在抛物面焦点 z = F

运行：
    .venv\\Scripts\\python.exe examples\\cassegrain_20lambda\\m3_cassegrain.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 中文字体，否则图里中文会变方框
matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

from common.cassegrain import Cassegrain, default_20lambda
from go_po.feed import cos_power_pattern
from go_po.geometry import build_polar_grid
from go_po.po import analytic_uniform_circular_aperture, far_field_aperture
from go_po.rays import feed_half_angle_for_edge, trace_rays


def db(x):
    return 20.0 * np.log10(np.maximum(np.abs(x), 1e-12))


def first_sidelobe_db(pat, theta):
    """先找到主瓣之后的第一个局部极小，再取其后的最大值。"""
    m = db(pat)
    for i in range(2, len(m) - 1):
        if m[i] <= m[i - 1] and m[i] < m[i + 1]:
            return float(m[i + 1:].max())
    return float("nan")


def aperture_from_rays(geo: Cassegrain, bundle, feed_pattern):
    """把射线结果转成口径幅度分布 A(rho)。

    口径功率守恒要求 |E_a| = |g| * sqrt(dOmega/dA)。对卡塞格伦，
    rho = 2*M*F*tan(theta/2)，故 dA/dOmega = (M*F)^2/cos^4(theta/2)，

        A = g(theta) * cos^2(theta/2) / (M*F)

    旧版写 g/(r1+r2) 只在单抛物面上成立（那边 r = F/cos^2(theta/2)）。
    """
    ok = bundle.valid
    rho = bundle.aperture_radius[ok]
    th = bundle.theta[ok]
    amp = feed_pattern(th) * np.cos(th / 2.0) ** 2 / (geo.M * geo.F)

    order = np.argsort(rho)
    rho_s, amp_s = rho[order], amp[order]
    # 把重复半径做平均，得到单调的插值表
    uniq, idx = np.unique(np.round(rho_s, 6), return_index=True)
    avg = np.array([amp_s[idx[i]:idx[i + 1]].mean() if i + 1 < len(idx) else amp_s[idx[i]:].mean()
                    for i in range(len(idx))])
    # 末端补一个点，保证插值覆盖到口径边缘。
    # 2026-11-05 审查 C-18：这里原来是把最外环常数外推（avg[-1]），
    # 但 θ 中点采样到不了 θe，外推值其实是 0.98R 处的；补解析值才对。
    theta_edge = float(feed_half_angle_for_edge(geo))
    edge_amp = float(
        feed_pattern(np.array([theta_edge]))[0]
        * np.cos(theta_edge / 2.0) ** 2 / (geo.M * geo.F)
    )
    uniq = np.append(uniq, geo.R)
    avg = np.append(avg, edge_amp)
    return uniq, avg


def main() -> int:
    geo = default_20lambda(M=4.0)
    print(geo.summary())
    print()

    # ---- GO 射线追踪 -------------------------------------------------------
    bundle = trace_rays(geo, n_radial=26, n_azimuth=72)
    print(f"[GO] 射线 {bundle.n_valid}/{bundle.n_rays} 完成全链路")
    print(f"[GO] 准直误差 (max |d_out.xy|) = {bundle.collimation_error():.3e}")
    print(f"[GO] 准直 RMS                  = {bundle.collimation_rms():.3e}")
    print(f"[GO] 光程起伏 r1+r2-z_main     = {bundle.path_length_ripple():.3e}")
    print()

    ok = bundle.valid
    rho_main = bundle.aperture_radius[ok]
    print(f"[GO] 主面落点半径 {rho_main.min():.3f} .. {rho_main.max():.3f} λ (口径半径 {geo.R})")
    r_sub = np.hypot(bundle.p_sub[:, 0], bundle.p_sub[:, 1])[ok]
    print(f"[GO] 副面落点半径 {r_sub.min():.3f} .. {r_sub.max():.3f} λ")
    print()

    # ---- 口径场 -> PO 远场 -------------------------------------------------
    k = 2.0 * np.pi / geo.lam
    q = 8.0
    feed = cos_power_pattern(q)
    rho_tab, amp_tab = aperture_from_rays(geo, bundle, feed)

    block_radius = float(geo.sub_rim()[0])
    print(f"[PO] 副面遮挡半径(近似) = {block_radius:.4f} λ  ->  面积比 {geo.blocking_ratio*100:.2f}%")

    def illum_blocked(rho):
        a = np.interp(np.asarray(rho, dtype=float), rho_tab, amp_tab)
        return np.where(np.asarray(rho, dtype=float) >= block_radius, a, 0.0)

    def illum_clear(rho):
        return np.interp(np.asarray(rho, dtype=float), rho_tab, amp_tab)

    theta = np.linspace(0.0, np.radians(10.0), 401)

    # 关键：有/无遮挡两条曲线必须用同一个归一化常数，
    # 否则遮挡会被归一化掉，增益损失恒为 0 dB。
    grid_rho, _, grid_dA = build_polar_grid(geo.R, 200, 320)
    ref_norm = float(np.sum(np.interp(grid_rho, rho_tab, amp_tab) * grid_dA))

    pat_blocked = far_field_aperture(geo.R, k, illum_blocked, theta,
                                     n_rho=200, n_phi=320, normalize=ref_norm)
    pat_clear = far_field_aperture(geo.R, k, illum_clear, theta,
                                   n_rho=200, n_phi=320, normalize=ref_norm)

    # 均匀口径做参考，确认积分器
    pat_ref = far_field_aperture(geo.R, k, lambda r: np.ones_like(r), theta, n_rho=200, n_phi=320)
    exact = analytic_uniform_circular_aperture(geo.R, k, theta)
    err_ref = np.max(np.abs(pat_ref - exact))
    print(f"[PO] 均匀口径自检: max|数值-解析| = {err_ref:.3e}")

    g_blocked = db(pat_blocked)[0]
    g_clear = db(pat_clear)[0]
    print(f"[PO] 遮挡造成的口径场幅损失 = {g_blocked - g_clear:.3f} dB"
      f"（场幅比；功率等效约 {(g_blocked - g_clear) / 2:.3f} dB）")
    print(f"[PO] 无遮挡第一副瓣 = {first_sidelobe_db(pat_clear, theta):.2f} dB")
    print(f"[PO] 有遮挡第一副瓣 = {first_sidelobe_db(pat_blocked, theta):.2f} dB")
    print()

    # ---- 出图 ---------------------------------------------------------------
    out_dir = ROOT / "data" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    np.savez(out_dir / "m3_cassegrain.npz", theta=theta, blocked=pat_blocked,
             clear=pat_clear, rho_tab=rho_tab, amp_tab=amp_tab,
             p_main=bundle.p_main[ok], p_sub=bundle.p_sub[ok],
             d_out=bundle.d_out[ok])

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    # (a) 子午面几何 + 射线
    a0 = ax[0, 0]
    zz = np.linspace(0.0, geo.z_rim_primary, 200)
    a0.plot(np.sqrt(4 * geo.F * zz), zz, "k-", lw=2, label="主面")
    a0.plot(-np.sqrt(4 * geo.F * zz), zz, "k-", lw=2)
    zs_sub = np.linspace(geo.L, float(geo.sub_rim()[1]), 120)
    rr_sub = np.array([geo.rim_radius_at(z) for z in zs_sub])
    a0.plot(rr_sub, zs_sub, "b-", lw=2, label="副面")
    a0.plot(-rr_sub, zs_sub, "b-", lw=2)
    for i in range(0, bundle.n_rays, 37):
        if not bundle.valid[i]:
            continue
        pts = np.vstack([[0, 0, geo.z_focus_far], bundle.p_sub[i],
                         bundle.p_main[i], bundle.p_main[i] + [0, 0, 3.0]])
        a0.plot(pts[:, 0], pts[:, 2], "-", color="tab:orange", lw=0.4, alpha=0.6)
    a0.axhline(geo.F, color="g", ls=":", lw=1)
    a0.plot([0], [geo.z_focus_far], "r*", ms=12, label="馈源")
    a0.plot([0], [geo.F], "go", ms=6, fillstyle="none", label="近焦点 (=F)")
    a0.set_xlim(-11, 11)
    a0.set_ylim(-1, 11)
    a0.set_xlabel("x / r  [λ]")
    a0.set_ylabel("z  [λ]")
    a0.set_title("(a) 子午面几何与射线路径")
    a0.grid(alpha=0.3)
    a0.legend(fontsize=8, loc="upper right")

    # (b) 口径幅度
    a1 = ax[0, 1]
    a1.plot(rho_tab, 20 * np.log10(np.maximum(amp_tab / amp_tab.max(), 1e-6)), "b-")
    a1.axvline(block_radius, color="r", ls="--", lw=1, label="副面遮挡边界")
    a1.set_xlabel("口径半径 rho  [λ]")
    a1.set_ylabel("相对幅度 [dB]")
    a1.set_title("(b) GO 口径幅度分布")
    a1.grid(alpha=0.3)
    a1.legend(fontsize=8)

    # (c) 方向图
    a2 = ax[1, 0]
    a2.plot(np.degrees(theta), db(pat_clear), "--", label="无遮挡")
    a2.plot(np.degrees(theta), db(pat_blocked), "-", label="含副面遮挡")
    a2.set_xlabel("theta [deg]")
    a2.set_ylabel("归一化 [dB]")
    a2.set_ylim(-45, 2)
    a2.set_title("(c) 卡塞格伦远场方向图")
    a2.grid(alpha=0.3)
    a2.legend(fontsize=8)

    # (d) 遮挡比对增益损失
    a3 = ax[1, 1]
    ratios = np.linspace(0.02, 0.16, 12)
    losses = []
    base = db(pat_clear)[0]
    for br in ratios:
        rad = np.sqrt(br) * geo.R
        p = far_field_aperture(geo.R, k,
                               lambda r: np.interp(np.asarray(r, float), rho_tab, amp_tab)
                               * (np.asarray(r, float) >= rad),
                               np.array([1e-6]), n_rho=180, n_phi=280, normalize=ref_norm)
        losses.append(db(p)[0])
    a3.plot(ratios * 100, losses, "o-")
    a3.axvline(geo.blocking_ratio * 100, color="r", ls="--", lw=1, label="本设计点")
    a3.set_xlabel("遮挡面积比 [%]")
    a3.set_ylabel("口径场幅损失 [dB]（场幅比，非增益）")
    a3.set_title("(d) 遮挡 -> 增益损失")
    a3.grid(alpha=0.3)
    a3.legend(fontsize=8)

    fig.tight_layout()
    fig_path = out_dir / "m3_cassegrain.png"
    fig.savefig(fig_path, dpi=140)
    print(f"saved: {fig_path}")
    return 0 if (bundle.collimation_error() < 1e-12
                 and bundle.path_length_ripple() < 1e-10
                 and err_ref < 5e-3) else 1


if __name__ == "__main__":
    raise SystemExit(main())
