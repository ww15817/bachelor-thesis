"""流水线：按配置把各环节串起来，一条命令跑完。

    几何建模 -> 射线追踪 -> 口径场 -> 远场积分 -> 保存结果

每一步都会记录校验指标，跑完统一汇报。
"""

from __future__ import annotations

import json
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from antennasim.config import CaseConfig
from common.cassegrain import Cassegrain
from go_po.feed import cos_power_pattern
from go_po.geometry import build_polar_grid
from go_po.po import analytic_uniform_circular_aperture, far_field_aperture
from go_po.rays import feed_half_angle_for_edge, trace_rays

__all__ = ["PipelineResult", "run_pipeline"]


@dataclass
class PipelineResult:
    case: str
    ok: bool
    seconds: float
    metrics: dict = field(default_factory=dict)
    outputs: dict = field(default_factory=dict)
    messages: list = field(default_factory=list)

    def summary(self) -> str:
        lines = [f"算例 {self.case}  用时 {self.seconds:.1f}s  状态 {'通过' if self.ok else '未通过'}"]
        for k, v in self.metrics.items():
            lines.append(f"    {k:<26} {v}")
        for k, v in self.outputs.items():
            lines.append(f"    输出 {k:<21} {v}")
        return "\n".join(lines)


def _build_geometry(cfg: CaseConfig) -> Cassegrain:
    g = cfg.geometry
    return Cassegrain(
        lam=1.0,
        D=g.diameter,
        F=cfg.focal_length,
        M=g.magnification,
        L=cfg.sub_vertex,
    )


def _aperture_from_rays(geo: Cassegrain, bundle, feed_pattern):
    """把射线结果转成口径幅度分布。

    口径平面的功率守恒要求 |E_a| = |g| * sqrt(dOmega/dA)。对卡塞格伦双反射面，
    口径半径与本几何的关系是 rho = 2*M*F*tan(theta/2)（数值验到 5e-15），
    于是 dA/dOmega = (M*F)^2 / cos^4(theta/2)，正确幅度是

        A = g(theta) * cos^2(theta/2) / (M*F)

    注意旧版写的是 g/(r1+r2)，那只在**单抛物面**上成立（因为 r = F/cos^2(theta/2)）。
    双反射面的等效扩散距离是 M*F，不是 r1+r2，2026-11-03 审查确认这里不守恒。
    """
    ok = bundle.valid
    rho = bundle.aperture_radius[ok]
    th = bundle.theta[ok]
    amp = feed_pattern(th) * np.cos(th / 2.0) ** 2 / (geo.M * geo.F)

    order = np.argsort(rho)
    rho_s, amp_s = rho[order], amp[order]
    uniq, idx = np.unique(np.round(rho_s, 6), return_index=True)
    avg = np.array([
        amp_s[idx[i]:idx[i + 1]].mean() if i + 1 < len(idx) else amp_s[idx[i]:].mean()
        for i in range(len(idx))
    ])
    return np.append(uniq, geo.R), np.append(avg, avg[-1])


def run_pipeline(
    cfg: CaseConfig,
    *,
    project_root: str | Path | None = None,
    make_figure: bool | None = None,
    verbose: bool = True,
) -> PipelineResult:
    """按配置跑完整条链路。"""
    t0 = time.time()
    root = Path(project_root) if project_root else Path.cwd()
    s = cfg.sampling
    messages: list[str] = []

    def say(msg: str) -> None:
        messages.append(msg)
        if verbose:
            print(msg)

    # ---------------------------------------------------------- 1 几何
    say("[1/5] 建模几何")
    geo = _build_geometry(cfg)
    geo.assert_consistency()
    say(f"      主面 D={geo.D:g}λ  F={geo.F:g}λ  M={geo.M:g}  L={geo.L:.4f}λ")
    say(f"      副面口径 {geo.d_rim:.4f}λ   遮挡比 {geo.blocking_ratio*100:.3f}%")

    # ---------------------------------------------------------- 2 射线
    say("[2/5] 射线追踪")
    bundle = trace_rays(geo, n_radial=s.n_radial, n_azimuth=s.n_azimuth)
    collim = bundle.collimation_error()
    ripple = bundle.path_length_ripple()
    say(f"      命中 {bundle.n_valid}/{bundle.n_rays}   准直误差 {collim:.2e}   光程起伏 {ripple:.2e}")

    # ---------------------------------------------------------- 3 口径
    say("[3/5] 构造口径场")
    feed = cos_power_pattern(cfg.feed.q)
    rho_tab, amp_tab = _aperture_from_rays(geo, bundle, feed)
    block_radius = float(geo.sub_rim()[0])

    # 边缘锥削是**结果**，不是输入：由 q、几何和双反射面的扩散因子共同决定。
    # 2026-11-03 审查第 3 条：配置里的 edge_taper_db 以前从没被读过，
    # YAML 写 -10 dB 而实际只有 -1.8 dB，谁也不提。现在把它当诊断量打印出来，
    # 配置里若写了期望值，差了就告警。
    taper_db = 20.0 * np.log10(
        max(abs(float(amp_tab[-1])) / max(abs(float(amp_tab[0])), 1e-12), 1e-12)
    )
    say(f"      口径边缘锥削 {taper_db:.2f} dB  (q={cfg.feed.q:g})")

    # 溢出发射：馈源波束比副面张角宽时，有一部分功率根本没打到副面。
    # 2026-11-04 WorkBuddy 审查指出：流水线用 ∫A dA 归一化，只认「打到口径的那部分功率」，
    # 溢出被约掉了，所以任何数字里都看不到它。q=8 时这一项是 -4.5 dB，比遮挡和锥削都大。
    #
    # 口径方向图的功率按 cos^q 计，打到副面的比例：
    #   ∫_0^{θe} cos^q sinθ dθ / ∫_0^{π/2} cos^q sinθ dθ = 1 - cos^(q+1)(θe)
    theta_edge = float(feed_half_angle_for_edge(geo))
    spill_eta = 1.0 - np.cos(theta_edge) ** (float(cfg.feed.q) + 1.0)
    say(
        f"      溢出发射效率 {spill_eta:.4f}  ({10.0 * np.log10(spill_eta):+.2f} dB)"
        f"   副面张角 {np.degrees(theta_edge):.2f}°"
    )
    if cfg.feed.edge_taper_db is not None:
        gap = abs(taper_db - float(cfg.feed.edge_taper_db))
        if gap > 1.0:
            say(
                f"      [警告] 配置期望 {float(cfg.feed.edge_taper_db):.1f} dB，"
                f"实际 {taper_db:.2f} dB，差 {gap:.1f} dB"
                f"（想要那个值得反解 q，不是写在这儿）"
            )

    # ---------------------------------------------------------- 4 远场
    say("[4/5] 远场积分")
    k = 2.0 * np.pi / geo.lam
    theta = np.linspace(0.0, np.radians(10.0), 401)

    grid_rho, _, grid_dA = build_polar_grid(geo.R, s.n_rho, s.n_phi)
    ref_norm = float(np.sum(np.interp(grid_rho, rho_tab, amp_tab) * grid_dA))

    def illum(rho):
        a = np.interp(np.asarray(rho, dtype=float), rho_tab, amp_tab)
        return np.where(np.asarray(rho, dtype=float) >= block_radius, a, 0.0)

    pat = far_field_aperture(geo.R, k, illum, theta,
                             n_rho=s.n_rho, n_phi=s.n_phi, normalize=ref_norm)

    # 均匀口径自检：数值 vs 解析艾里图案
    pat_ref = far_field_aperture(geo.R, k, lambda r: np.ones_like(r), theta,
                                 n_rho=s.n_rho, n_phi=s.n_phi)
    exact = analytic_uniform_circular_aperture(geo.R, k, theta)
    airy_err = float(np.max(np.abs(pat_ref - exact)))
    say(f"      均匀口径自检  max|数值-解析| = {airy_err:.2e}")

    peak_db = 20.0 * np.log10(max(abs(pat[0]), 1e-12))
    say(f"      峰值增益损失 {peak_db:.3f} dB")

    metrics = {
        "近焦点误差": f"{abs(geo.z_focus_near - geo.F):.2e}",
        "遮挡面积比": f"{geo.blocking_ratio*100:.3f} %",
        "射线命中率": f"{bundle.n_valid}/{bundle.n_rays}",
        "准直误差": f"{collim:.2e}",
        "光程起伏": f"{ripple:.2e}",
        "艾里自检偏差": f"{airy_err:.2e}",
        "峰值增益损失": f"{peak_db:.3f} dB",
    }

    checks = {
        "几何自检": True,
        "射线全部命中": bundle.n_valid == bundle.n_rays,
        "准直误差达标": collim < 1e-12,
        "等光程达标": ripple < 1e-10,
        "艾里自检达标": airy_err < 5e-3,
    }
    ok = all(checks.values())

    # ---------------------------------------------------------- 5 输出
    say("[5/5] 保存结果")
    out_dir = root / cfg.output.directory
    out_dir.mkdir(parents=True, exist_ok=True)
    tag = cfg.output.tag

    npz_path = out_dir / f"{tag}.npz"
    np.savez(npz_path, theta=theta, pattern=pat, rho=rho_tab, amplitude=amp_tab,
             z_focus_near=geo.z_focus_near, blocking_ratio=geo.blocking_ratio)

    json_path = out_dir / f"{tag}.json"
    json_path.write_text(
        json.dumps(
            {"case": cfg.name, "metrics": metrics, "checks": checks,
             "config": asdict(cfg)},
            ensure_ascii=False, indent=2,
        ),
        encoding="utf-8",
    )

    outputs = {"数据": str(npz_path), "指标": str(json_path)}

    do_fig = cfg.output.figure if make_figure is None else make_figure
    if do_fig:
        try:
            import matplotlib
            matplotlib.use("Agg")
            matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
            matplotlib.rcParams["axes.unicode_minus"] = False
            import matplotlib.pyplot as plt

            fig, ax = plt.subplots(1, 2, figsize=(11, 4))
            ax[0].plot(rho_tab, 20 * np.log10(np.maximum(amp_tab / amp_tab.max(), 1e-6)))
            ax[0].axvline(block_radius, color="r", ls="--", lw=1, label="副面遮挡")
            ax[0].set_xlabel("口径半径 [λ]"); ax[0].set_ylabel("相对幅度 [dB]")
            ax[0].set_title(f"{cfg.name} 口径幅度"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)

            ax[1].plot(np.degrees(theta), 20 * np.log10(np.maximum(abs(pat), 1e-9)))
            ax[1].set_xlabel("θ [deg]"); ax[1].set_ylabel("归一化 [dB]"); ax[1].set_ylim(-45, 2)
            ax[1].set_title("远场方向图（含遮挡）"); ax[1].grid(alpha=0.3)
            fig.tight_layout()
            fig_path = out_dir / f"{tag}.png"
            fig.savefig(fig_path, dpi=130)
            plt.close(fig)
            outputs["图"] = str(fig_path)
        except Exception as exc:  # 出图失败不影响主流程
            messages.append(f"      出图跳过：{exc}")

    for v in outputs.values():
        say(f"      {v}")

    return PipelineResult(
        case=cfg.name, ok=ok, seconds=time.time() - t0,
        metrics=metrics, outputs=outputs, messages=messages,
    )
