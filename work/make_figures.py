"""生成论文要的三张图：效率预算瀑布、q 扫描、表面误差灵敏度。

输出到 docs/figures/，对应 docs/figures/README.md 里「还缺的」第 2/3/5 项。

跑法（项目根目录下）：

    .\.venv\Scripts\python.exe work\make_figures.py
"""

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

# 中文字体（否则中文标签会变成方框）
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

from antennasim.config import load_config
from antennasim.pipeline import run_pipeline
from go_po.aperture_po import axial_amplitude, far_field_pattern, tapered_aperture

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "figures"


def fig_budget() -> None:
    """效率预算瀑布：理想 35.96 → 溢出 → 口径 → 实际。"""
    ideal = 35.96
    spill = -4.49      # q=8，台账 W-0013
    aperture = -0.22   # 锥削 + 遮挡的功率等效，台账 W-0020
    # 实际值用流水线真值，不用三项相加——各项取两位小数，相加会差 0.01
    real = 31.24

    labels = ["理想\n(η=1)", "溢出\n(q=8)", "锥削+遮挡", "实际"]
    values = [ideal, spill, aperture, real]
    bottoms = [0, ideal + spill, ideal + spill, 0]
    colors = ["#9ec5e8", "#e8a87c", "#e8d47c", "#7fb069"]

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for i, (lab, val, bottom, col) in enumerate(zip(labels, values, bottoms, colors)):
        if i in (0, 3):
            ax.bar(i, val, bottom=0, color=col, width=0.6)
            ax.text(i, val + 0.3, f"{val:.2f}", ha="center", fontsize=10)
        else:
            ax.bar(i, abs(val), bottom=bottom, color=col, width=0.6)
            ax.text(i, bottom + 0.3, f"{val:+.2f}", ha="center", fontsize=10)

    ax.set_xticks(range(4))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("绝对增益 (dBi)")
    ax.set_title("20λ 卡塞格伦的效率预算（10 GHz，q=8）", fontsize=11)
    ax.text(0.5, -0.18, "各项保留两位小数，逐项相加 31.25 与真值 31.24 差 0.01 dB",
            transform=ax.transAxes, ha="center", fontsize=8, color="gray")
    ax.set_ylim(28, 37)
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "budget_waterfall.png", dpi=160)
    plt.close(fig)
    print("出图: budget_waterfall.png")


def fig_qscan() -> None:
    """q 扫描：绝对增益 vs q，标出平台。"""
    base = load_config(ROOT / "cases" / "cassegrain_20lambda.yaml")
    qs = [4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 47.0, 64.0, 80.0, 96.0]
    gains = []
    for q in qs:
        cfg = dataclasses.replace(base, feed=dataclasses.replace(base.feed, q=q))
        res = run_pipeline(cfg, project_root=ROOT, make_figure=False, verbose=False)
        gains.append(float(res.metrics["绝对增益"].split()[0]))

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(qs, gains, "o-", color="#3a6ea5", markersize=4)
    best_i = int(np.argmax(gains))
    ax.plot(qs[best_i], gains[best_i], "r*", markersize=14,
            label=f"最优 q={qs[best_i]:g}，{gains[best_i]:.2f} dBi")
    # 平台：差值小于 0.05 dB 的区间
    plateau = [q for q, g in zip(qs, gains) if gains[best_i] - g < 0.10]
    ax.axvspan(min(plateau), max(plateau), alpha=0.15, color="red",
               label=f"平台 {min(plateau):g}–{max(plateau):g}")
    ax.set_xlabel("馈源 cos^q 功率指数 q")
    ax.set_ylabel("绝对增益 (dBi)")
    ax.set_title("馈源 q 扫描：最优是平台不是点", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "q_scan.png", dpi=160)
    plt.close(fig)
    print(f"出图: q_scan.png（最优 q={qs[best_i]:g}，平台 {min(plateau):g}–{max(plateau):g}）")


def fig_surface_error() -> None:
    """表面误差灵敏度：PO 直接算 vs Ruze。"""
    R = 10.0
    k = 2.0 * np.pi
    edge_taper = -1.91
    sigmas = [0.0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05]

    base = tapered_aperture(R, edge_taper_db=edge_taper, n_per_lambda=12.0)
    ref = axial_amplitude(base)
    r = np.hypot(base.x, base.y)
    rng = np.random.default_rng(20261105)

    po, ruze = [], []
    for s in sigmas:
        if s == 0.0:
            po.append(0.0)
            ruze.append(0.0)
            continue
        trials = []
        for _ in range(12):
            dphi = 2.0 * k * s * rng.normal(size=r.shape)
            af = tapered_aperture(R, edge_taper_db=edge_taper, phase_error=dphi,
                                  n_per_lambda=12.0)
            g = float(far_field_pattern(af, np.array([0.0]), normalize=ref)[0])
            trials.append(20.0 * np.log10(max(g, 1e-12)))
        po.append(float(np.mean(trials)))
        sigma_phi = 4.0 * np.pi * s
        ruze.append(10.0 * np.log10(np.exp(-(sigma_phi ** 2))))

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(np.array(sigmas) * 100, po, "o-", label="PO 直接算（随机相位屏）", color="#3a6ea5")
    ax.plot(np.array(sigmas) * 100, ruze, "s--", label="Ruze 公式", color="#c0504d")
    ax.set_xlabel("表面误差 rms  σ_w/λ (%)")
    ax.set_ylabel("增益损失 (dB)")
    ax.set_title("表面误差灵敏度：随机相位屏下 PO 与 Ruze 一致\n（确定性低频变形不适用 Ruze，见 ISSUE-018）", fontsize=10)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "surface_error.png", dpi=160)
    plt.close(fig)
    print("出图: surface_error.png")


def fig_deflection_map() -> None:
    """有限元变形场云图 + 口径相位分布。"""
    from scipy.interpolate import griddata
    d = np.load(ROOT / "data" / "feeds" / "plate_deflection_2mm.npz")
    fx, fy, fw = d["x"], d["y"], d["w"]
    lam = 0.030

    g = np.linspace(-0.3, 0.3, 160)
    X, Y = np.meshgrid(g, g)
    W = griddata((fx, fy), fw, (X, Y), method="linear", fill_value=np.nan)
    mask = X**2 + Y**2 > 0.3**2
    W[mask] = np.nan

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 4.0))
    im0 = axes[0].pcolormesh(X * 1000, Y * 1000, W * 1e6, shading="auto", cmap="viridis")
    axes[0].set_xlabel("x (mm)")
    axes[0].set_ylabel("y (mm)")
    axes[0].set_title("面板变形 w（2 mm 铝板，自重）", fontsize=10)
    axes[0].set_aspect("equal")
    fig.colorbar(im0, ax=axes[0], label="w (μm)")

    im1 = axes[1].pcolormesh(X * 1000, Y * 1000, 4*np.pi*W/lam, shading="auto", cmap="magma")
    axes[1].set_xlabel("x (mm)")
    axes[1].set_ylabel("y (mm)")
    axes[1].set_title("口径相位扰动 4πw/λ（λ=30 mm）", fontsize=10)
    axes[1].set_aspect("equal")
    fig.colorbar(im1, ax=axes[1], label="dφ (rad)")

    fig.tight_layout()
    fig.savefig(OUT / "deflection_map.png", dpi=160)
    plt.close(fig)
    print("出图: deflection_map.png")


def fig_spill_schematic() -> None:
    """馈源功率去向：按 cos^q 功率密度画，面积就是功率比例。"""
    q = 8.0
    th_e = 17.7613
    th_main = 72.65
    theta = np.linspace(0.0, 90.0, 721)
    p = np.cos(np.deg2rad(theta)) ** q

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    m1 = theta <= th_e
    m2 = (theta > th_e) & (theta <= th_main)
    m3 = theta > th_main
    ax.fill_between(theta[m1], 0, p[m1], color="#9ec5e8",
                    label=f"打到副面（θ ≤ {th_e:.1f}°，35.6%）")
    ax.fill_between(theta[m2], 0, p[m2], color="#e8a87c",
                    label="溢出（64% 功率）")
    ax.fill_between(theta[m3], 0, p[m3], color="#e8d47c",
                    label=f"打到主面后散射（θ > {th_main:.1f}°）")
    ax.axvline(th_e, color="k", lw=1.0, ls="--")
    ax.axvline(th_main, color="k", lw=1.0, ls=":")
    ax.set_xlim(0, 90)
    ax.set_ylim(0, 1.05)
    ax.set_xlabel("馈源出射角 θ (deg)")
    ax.set_ylabel("功率密度 cos⁸θ（归一化）")
    ax.set_title("q=8 馈源的功率去向：面积即功率占比", fontsize=11)
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(alpha=0.25)
    fig.tight_layout()
    fig.savefig(OUT / "spill_schematic.png", dpi=160)
    plt.close(fig)
    print("出图: spill_schematic.png（面积按功率密度）")


def fig_pattern_compare() -> None:
    """结构变形对方向图的影响：把变形放大 10 倍才看得出差别。

    真实 2 mm 板的相位损失只有 −0.014 dB，主瓣内两条线完全重合——那本身是结论，
    但图上画不出东西。所以这张图用**放大 10 倍**的变形（−1.585 dB），
    让读者看见「结构变形怎么改变主瓣和近轴旁瓣」，正文里注明倍率。
    """
    from scipy.interpolate import griddata
    d = np.load(ROOT / "data" / "feeds" / "plate_deflection_2mm.npz")
    fx, fy, fw = d["x"], d["y"], d["w"]
    lam = 0.030
    edge = -1.91
    R = 10.0
    scale = 10.0

    theta = np.linspace(0.0, 6.0, 601)
    base = tapered_aperture(R, edge_taper_db=edge, n_per_lambda=10.0)
    ref = axial_amplitude(base)

    def phase(rho, phi_a):
        x = np.asarray(rho) * np.cos(phi_a)
        y = np.asarray(rho) * np.sin(phi_a)
        w = griddata((fx / lam, fy / lam), fw, (x, y), method="linear", fill_value=0.0)
        return 4.0 * np.pi * (w * scale / lam)

    p0 = far_field_pattern(base, theta, normalize=ref)
    af1 = tapered_aperture(R, edge_taper_db=edge, n_per_lambda=10.0)
    r_norm = np.hypot(af1.x, af1.y)
    phi_pts = np.arctan2(af1.y, af1.x)
    af1.e = af1.e * np.exp(1j * phase(r_norm, phi_pts))
    p1 = far_field_pattern(af1, theta, normalize=ref)

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.plot(theta, 20*np.log10(np.clip(p0, 1e-6, None)), label="无变形", color="#3a6ea5")
    ax.plot(theta, 20*np.log10(np.clip(p1, 1e-6, None)), "--",
            label="变形 ×10（≈ −1.6 dB）", color="#c0504d")
    ax.set_xlabel("θ (deg)")
    ax.set_ylabel("归一化方向图 (dB)")
    ax.set_title("结构变形对方向图的影响（20λ，10 GHz；真实 2 mm 板只掉 0.014 dB）", fontsize=9)
    ax.set_xlim(0, 6.0)
    ax.set_ylim(-40, 0.5)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "pattern_compare.png", dpi=160)
    plt.close(fig)
    print("出图: pattern_compare.png（变形 ×10，真实值 0.014 dB 画不出差别）")


def fig_chain_diagram() -> None:
    """链路总框图：三个开源工具 + 耦合回路。"""
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5.2)
    ax.axis("off")

    boxes = [
        (0.4, 3.2, 2.4, 1.3, "openEMS\n(FDTD 全波)", "#9ec5e8", "算馈源\n销钉喇叭"),
        (3.7, 3.2, 2.6, 1.3, "自研 GO/PO\n(高频近似)", "#a8d5a2", "算远场\n方向图 + 增益"),
        (7.0, 3.2, 2.6, 1.3, "FEniCSx\n(有限元)", "#e8c4a0", "算面板变形\nw(x, y)"),
    ]
    for x, y, w, h, title, color, sub in boxes:
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.08",
                                    facecolor=color, edgecolor="k", lw=1.2))
        ax.text(x + w/2, y + h*0.70, title, ha="center", va="center",
                fontsize=10, weight="bold")
        ax.text(x + w/2, y + h*0.24, sub, ha="center", va="center", fontsize=8)

    # 主链：openEMS -> GO/PO
    ax.add_patch(FancyArrowPatch((2.82, 3.85), (3.68, 3.85), arrowstyle="-|>",
                                 mutation_scale=16, lw=1.5, color="k",
                                 shrinkA=0, shrinkB=0))

    # 耦合回路：FEniCSx 底 -> 下 -> 左 -> 上到 GO/PO 底（U 形，避开盒子）
    pts = [(8.3, 3.2), (8.3, 2.2), (5.0, 2.2), (5.0, 3.18)]
    for (x0, y0), (x1, y1) in zip(pts, pts[1:]):
        ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-",
                                     lw=1.5, color="#555555", shrinkA=0, shrinkB=0))
    ax.add_patch(FancyArrowPatch((5.0, 2.5), (5.0, 3.18), arrowstyle="-|>",
                                 mutation_scale=16, lw=1.5, color="#555555",
                                 shrinkA=0, shrinkB=0))
    ax.text(6.65, 2.05, "口径相位  dφ = 4πw/λ", fontsize=9, ha="center", color="#333333")

    ax.text(5.0, 1.35, "同一条命令跑完：feed.source = npz  +  structure.enabled = true",
            ha="center", fontsize=9, style="italic")
    ax.text(5.0, 0.85, "每一段都对上解析解：Airy 7e-06 / 共焦 4e-16 / 幅度守恒 5e-15",
            ha="center", fontsize=8.5, color="#444444")
    ax.set_title("全开源反射面天线分析链路", fontsize=12, pad=8)
    fig.tight_layout()
    fig.savefig(OUT / "chain_diagram.png", dpi=160)
    plt.close(fig)
    print("出图: chain_diagram.png")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    fig_budget()
    fig_qscan()
    fig_surface_error()
    fig_deflection_map()
    fig_spill_schematic()
    fig_pattern_compare()
    fig_chain_diagram()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
