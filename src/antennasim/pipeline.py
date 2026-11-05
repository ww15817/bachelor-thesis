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
    # 表尾补 R 处的**解析**幅度，而不是把最外环常数外推过去。
    # 2026-11-05 审查 C-18：θ 用中点采样（(i+0.5)/n·θmax）永远到不了 θe，
    # 表尾那个值其实是 0.98R 处的，把「口径边缘锥削」印成 −1.83 dB，
    # 而解析真值（cos⁴(θe)·cos²(θe/2)）是 −1.91 dB，系统性偏亮 0.07 dB。
    # 幅度模型本来就是 θ 的解析函数，直接用 θe 取一次即可。
    theta_edge = float(feed_half_angle_for_edge(geo))
    edge_amp = float(
        feed_pattern(np.array([theta_edge]))[0]
        * np.cos(theta_edge / 2.0) ** 2 / (geo.M * geo.F)
    )
    return np.append(uniq, geo.R), np.append(avg, edge_amp)


def _load_feed_pattern(cfg: CaseConfig, root: Path):
    """按 ``cfg.feed.source`` 返回馈源电压方向图 ``g(theta)``。

    - ``analytic``：解析 ``cos^q``
    - ``npz``：openEMS 全波算出的方向图（P2）。读进来后按 theta=0 归一到 1，
      再做成插值函数，这样下游 `_aperture_from_rays` 一行都不用改。

    ⚠️ 矩形喇叭的 E 面和 H 面锥削差很多（ISSUE-020），而流水线假设轴对称。
    选哪个面由 ``cfg.feed.plane`` 决定，两个面会给出不同的增益，报告里要写清。
    """
    if cfg.feed.source == "analytic":
        return cos_power_pattern(cfg.feed.q)

    if cfg.feed.source != "npz":
        raise ValueError(f"未知的 feed.source：{cfg.feed.source}（只认 analytic / npz）")
    if not cfg.feed.file:
        raise ValueError("feed.source=npz 时必须给 feed.file")

    path = Path(root) / cfg.feed.file
    if not path.exists():
        raise FileNotFoundError(f"找不到馈源方向图文件：{path}")

    data = np.load(path, allow_pickle=True)
    theta_deg = np.asarray(data["theta"], dtype=float)
    e_norm = np.asarray(data["E_norm"], dtype=float)
    phi = np.asarray(data["phi"], dtype=float)

    # phi=0 是 x-z 面（H 面，受口径 x 向宽度控制）；phi=90 是 y-z 面（E 面）
    want = 90.0 if str(cfg.feed.plane).upper() == "E" else 0.0
    iphi = int(np.argmin(np.abs(np.mod(phi, 360.0) - want)))

    e = np.abs(e_norm[:, iphi])
    e = e / max(float(e[0]), 1e-12)        # 归一到 theta=0
    th_rad = np.deg2rad(theta_deg)

    return lambda th: np.interp(th, th_rad, e)


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
    feed = _load_feed_pattern(cfg, root)
    if cfg.feed.source == "npz":
        say(f"      馈源：{cfg.feed.file}（{cfg.feed.plane} 面）")
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

    # ---------------------------------------------------------- 结构侧
    # 把 FEniCSx 算出来的面板变形接进口径相位（dφ = 2kw = 4πw/λ）。
    # 变形场是散点（板网格节点），要插值到口径的极坐标网格上。
    phase_func = None
    if cfg.structure.enabled:
        if not cfg.structure.file:
            raise ValueError("structure.enabled=true 时必须给 structure.file")
        st_path = Path(root) / cfg.structure.file
        if not st_path.exists():
            raise FileNotFoundError(f"找不到变形场文件：{st_path}")
        dfl = np.load(st_path)
        fx, fy, fw = dfl["x"], dfl["y"], dfl["w"]
        lam_st = float(cfg.structure.lambda_m)
        from scipy.interpolate import griddata

        def phase_func(rho, phi_a):  # noqa: F811 - 有意覆盖
            # 口径网格是 λ 单位，变形场是米，两边都换成 λ 再插值
            x_lam = np.asarray(rho) * np.cos(phi_a)
            y_lam = np.asarray(rho) * np.sin(phi_a)
            w_m = griddata((fx / lam_st, fy / lam_st), fw,
                           (x_lam, y_lam), method="linear", fill_value=0.0)
            return 4.0 * np.pi * (w_m / lam_st)

        say(f"      结构：{cfg.structure.file}（{fw.size} 节点，λ={lam_st*1000:.1f} mm）")

    def illum(rho):
        a = np.interp(np.asarray(rho, dtype=float), rho_tab, amp_tab)
        return np.where(np.asarray(rho, dtype=float) >= block_radius, a, 0.0)

    pat = far_field_aperture(geo.R, k, illum, theta,
                             n_rho=s.n_rho, n_phi=s.n_phi, normalize=ref_norm,
                             phase_func=phase_func)

    # 均匀口径自检：数值 vs 解析艾里图案
    pat_ref = far_field_aperture(geo.R, k, lambda r: np.ones_like(r), theta,
                                 n_rho=s.n_rho, n_phi=s.n_phi)
    exact = analytic_uniform_circular_aperture(geo.R, k, theta)
    airy_err = float(np.max(np.abs(pat_ref - exact)))
    say(f"      均匀口径自检  max|数值-解析| = {airy_err:.2e}")

    peak_db = 20.0 * np.log10(max(abs(pat[0]), 1e-12))
    # 2026-11-04 审查 C-09：这个数不随网格加密收敛（n_rho 200/400/800/1600 给出
    # -0.4786 / -0.4898 / -0.4842 / -0.4870，抖动 ±0.006）。原因是遮挡台阶
    # rho=Rb 落在网格环之间，整个环被整体判给或判走。所以只给两位小数，
    # 论文里写「约 -0.49 dB」；要三位小数得对 Rb 所在环做亚网格面积加权。
    say(f"      口径场幅损失 {peak_db:.2f} dB（遮挡台阶在网格之间，只到 ±0.01 dB）")

    # 效率预算：把溢出、照明、遮挡三项并起来，给绝对 dBi。
    # 2026-11-04 审查指出「-0.479 dB」只是遮挡 + 锥削的相对损失，
    # 溢出发射那 -4.49 dB 一直不在账上。这里补全。
    #
    # 口径照明效率（标准定义，均匀照明时 = 1）：
    #   η_illum = |∫A dA|² / (面积 · ∫|A|² dA)
    amp_on_grid = np.interp(grid_rho, rho_tab, amp_tab)
    aperture_area = np.pi * geo.R**2
    sum_a = float(np.sum(amp_on_grid * grid_dA))
    sum_a2 = float(np.sum(amp_on_grid**2 * grid_dA))
    eta_illum = sum_a**2 / (aperture_area * sum_a2)
    eta_block = 1.0 - float(geo.blocking_ratio)

    # 归一约定（2026-11-05 审查 C-11 更正）：
    #   方向性 D = 4π|∫A dA|² / (λ² ∫|A|² dA) —— 分母是 ∫|A|²，不是 (∫A)²/A。
    #   所以效率必须回到二重积分里算，**不能从场幅比 peak_db 反推**：
    #   peak_db 已经是「未遮挡/全口径」的场幅比，再平方等于把遮挡算了两遍
    #   （均匀口径下给出 (1-b)²，而教科书遮挡口径是 1-b，正好差一倍 dB）。
    #   正确做法：分子和分母都只取未遮挡区。
    amp_un = np.where(grid_rho >= block_radius, amp_on_grid, 0.0)
    eta_taper_block = float(np.sum(amp_un * grid_dA)) ** 2 / (
        aperture_area * float(np.sum(amp_un**2 * grid_dA))
    )

    # 结构侧：相位扰动对轴向场的影响。
    # 上面那个 eta_taper_block 只用了**幅度**，不含相位；而结构变形是纯相位扰动，
    # 光靠幅度预算看不出它。这里用同一次归一化，比一次「有相位 / 无相位」的轴向幅度。
    eta_structure = 1.0
    if phase_func is not None:
        pat_no_phase = far_field_aperture(
            geo.R, k, illum, np.array([0.0]),
            n_rho=s.n_rho, n_phi=s.n_phi, normalize=ref_norm,
        )
        ratio = float(abs(pat[0]) / max(abs(pat_no_phase[0]), 1e-12))
        eta_structure = ratio ** 2
        say(f"      结构相位损失 {10.0 * np.log10(max(eta_structure, 1e-12)):.3f} dB"
            f"（轴向场比 {ratio:.5f}）")

    eta_total = spill_eta * eta_taper_block * eta_structure
    # 理想均匀口径（η=1）的增益：(πD/λ)²
    gain_ideal_db = 10.0 * np.log10((np.pi * geo.D / geo.lam) ** 2)
    gain_db = gain_ideal_db + 10.0 * np.log10(max(eta_total, 1e-12))
    say(
        f"      口径效率 {eta_taper_block:.4f}"
        f"（锥削 {eta_illum:.4f} × 遮挡 {eta_block:.4f}，诊断用）"
        f"  → 总效率 {eta_total:.4f}  绝对增益 {gain_db:.2f} dBi"
        f"  (理想 {gain_ideal_db:.2f} dBi)"
    )



    metrics = {
        "近焦点误差": f"{abs(geo.z_focus_near - geo.F):.2e}",
        "遮挡面积比": f"{geo.blocking_ratio*100:.3f} %",
        "射线命中率": f"{bundle.n_valid}/{bundle.n_rays}",
        "准直误差": f"{collim:.2e}",
        "光程起伏": f"{ripple:.2e}",
        "艾里自检偏差": f"{airy_err:.2e}",
        "口径场幅损失": f"{peak_db:.2f} dB",
        "溢出发射效率": f"{spill_eta:.4f}",
        "照明效率": f"{eta_illum:.4f}",
        "总效率": f"{eta_total:.4f}",
        "绝对增益": f"{gain_db:.2f} dBi",
        "理想增益": f"{gain_ideal_db:.2f} dBi",
    }

    checks = {
        "几何自检": True,
        "射线全部命中": bundle.n_valid == bundle.n_rays,
        "准直误差达标": collim < 1e-12,
        "等光程达标": ripple < 1e-10,
        "艾里自检达标": airy_err < 5e-3,
        # 2026-11-04 审查 C-03：前五项验的全是「几何和积分器对不对」，
        # 75% 功率溢出的算例照样判「通过」。这条管能量分配。
        # 阈值 0.2 只挡明显病态——fd05 对照算例（溢出 0.2452）要放行，
        # 它是故意做差的对照，不是错误。
        "溢出可接受": bool(spill_eta > 0.2),
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
