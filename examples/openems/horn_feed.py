"""用 openEMS 算一个角锥喇叭馈源，导出方向图，跟解析 cos^q 对比。

这是 P2（openEMS 馈源接入）的第一步：先把「全波算馈源」这条链路跑通，
再谈把它接进 GO/PO 流水线。

模型：WR-90 矩形波导（22.86×10.16 mm，TE10 激励）+ 阶梯喇叭（三段递扩，
最终口径 80×60 mm），中心频率 10 GHz。阶梯是为了避开斜面的建模麻烦，
方向图会有一点台阶引起的畸变，但作为「集成演示」够用；要精确设计请换
官方 Horn_Antenna 教程的斜面模型。

跑法（WSL 里）：

    source ~/openems-install/venv/bin/activate
    python /mnt/c/.../examples/openems/horn_feed.py

输出：
    data/feeds/horn_10ghz.npz   （theta/phi/E_norm/Prad/Dmax）

注意：
- 网格分辨率 λ/20，域边界离喇叭口 ≥ λ/2
- 端口必须在网格内部，且两端只沿 z 变化
"""

import os
import tempfile

import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS as OpenEMS
from openEMS.physical_constants import C0

# ---------------------------------------------------------------- 参数
unit = 1e-3
f0 = 10e9
lam = C0 / f0 / unit          # 30 mm
RES = lam / 20.0

a = 22.86                     # WR-90 宽边
b = 10.16                     # WR-90 窄边
t = 0.5                       # 壁厚

z_wg0, z_wg1 = -45.0, 0.0     # 波导段
steps = [                      # (z_start, z_end, 宽, 高)
    (0.0, 20.0, 36.0, 22.0),
    (20.0, 40.0, 58.0, 40.0),
    (40.0, 60.0, 80.0, 60.0),
]

print(f"f0 = {f0/1e9:.2f} GHz, lambda = {lam:.2f} mm, RES = {RES:.3f} mm")

# ---------------------------------------------------------------- FDTD
FDTD = OpenEMS(NrTS=30000, EndCriteria=1e-4)
FDTD.SetGaussExcite(f0, f0 / 2)
FDTD.SetBoundaryCond(["PML_8"] * 6)

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

metal = CSX.AddMetal("wall")

# ---------------------------------------------------------------- 几何
def add_tube(z0, z1, w, h, thick=t):
    """四壁围成的矩形管（PEC 薄板）。

    壁厚向外长，**内腔保持 w×h**——否则端口尺寸 a×b 跟腔体对不上，
    TE10 模式激励不出来（Prad 会一直是 0）。
    """
    hw, hh = w / 2.0 + thick, h / 2.0 + thick
    metal.AddBox([-hw, -hh, z0], [hw, -hh + thick, z1])      # 下壁
    metal.AddBox([-hw, hh - thick, z0], [hw, hh, z1])        # 上壁
    metal.AddBox([-hw, -hh, z0], [-hw + thick, hh, z1])      # 左壁
    metal.AddBox([hw - thick, -hh, z0], [hw, hh, z1])        # 右壁


add_tube(z_wg0, z_wg1, a, b)
# 阶梯之间用「封口环」补上台阶（喇叭段本身也建管，台阶面靠相邻管的端面遮挡）
for z0, z1, w, h in steps:
    add_tube(z0, z1, w, h, thick=t)

# 喇叭外壁到远端边界之间留空
z_end = steps[-1][1]

# ---------------------------------------------------------------- 网格
margin = 1.0 * lam
mesh.AddLine("x", [-margin, -a / 2, a / 2, margin])
mesh.AddLine("y", [-margin, -b / 2, b / 2, margin])
z_lines = [z_wg0 - margin, z_wg0, z_wg1,
           z_wg0 + 2.0, z_wg0 + 6.0] + [s[1] for s in steps] + [z_end + margin]
mesh.AddLine("z", sorted(set(z_lines)))
# 喇叭里的台阶尺寸也要进网格线，否则壁面会被网格跨越
for z0, z1, w, h in steps:
    mesh.AddLine("x", [-w / 2, w / 2])
    mesh.AddLine("y", [-h / 2, h / 2])
mesh.SmoothMeshLines("all", RES, 1.4)
print(f"网格: x={mesh.GetQtyLines(0)} y={mesh.GetQtyLines(1)} z={mesh.GetQtyLines(2)}")

# ---------------------------------------------------------------- 端口
# 端口盒：start/stop 给对角两点，p_dir 方向必须**有长度**（零长度会报
# "Port length in excitation direction may not be zero"），这里给 4 mm。
# 端口贴在波导起始面上（z_wg0），不要往腔里挪——挪进去激励不出来（Prad≈0）。
port = FDTD.AddRectWaveGuidePort(
    1,
    np.array([-a / 2, -b / 2, z_wg0]),
    np.array([a / 2, b / 2, z_wg0 + 4.0]),
    "z",
    a,
    b,
    "TE10",
    1.0,
)

# ---------------------------------------------------------------- 近远场
nf2ff = FDTD.CreateNF2FFBox("nf2ff", opt_resolution=[RES, RES, RES], frequency=[f0])

theta = np.arange(0, 91, 1.0)
phi = np.arange(0, 181, 2.0)

# ---------------------------------------------------------------- 运行
sim_path = os.path.join(tempfile.gettempdir(), "openems_horn_feed")
os.makedirs(sim_path, exist_ok=True)
print(f"工作目录: {sim_path}")
FDTD.Run(sim_path, cleanup=True, verbose=0)

# ---------------------------------------------------------------- 后处理
nf2ff_res = nf2ff.CalcNF2FF(
    sim_path, f0, np.deg2rad(theta), np.deg2rad(phi),
    outfile=os.path.join(sim_path, "nf2ff.h5"),
)

E_norm = np.squeeze(nf2ff_res.E_norm[0])
Prad = float(np.squeeze(nf2ff_res.Prad[0]))
Dmax = float(np.squeeze(nf2ff_res.Dmax[0]))
print()
print(f"Prad = {Prad:.4e} W")
print(f"Dmax = {Dmax:.4f} ({10*np.log10(Dmax):.2f} dBi)")

# E 面（phi=0）归一化方向图，拟合 cos^q
i_phi0 = 0
e_face = np.abs(E_norm[:, i_phi0])
e_face = e_face / np.max(e_face)
mask = (theta >= 0) & (theta <= 60) & (e_face > 0)
fit_theta = np.deg2rad(theta[mask])
fit_e = e_face[mask]
# 电压方向图 cos^(q/2) → 线性拟合 log(cos) vs log(E)
q_half, _ = np.polyfit(np.log(np.cos(fit_theta)), np.log(fit_e), 1)
print(f"E 面拟合 q ≈ {2.0 * q_half:.2f}   (cos^q 功率指数)")
for t_deg in (10, 20, 30, 40):
    i = int(np.argmin(np.abs(theta - t_deg)))
    print(f"  theta={t_deg:2d}°  E/E0 = {e_face[i]:.4f}  ({20*np.log10(max(e_face[i],1e-6)):+.2f} dB)")

out_dir = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "data", "feeds")
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, "horn_10ghz.npz")
np.savez(out, theta=theta, phi=phi, E_norm=E_norm, Prad=Prad, Dmax=Dmax,
         q_est=2.0 * q_half, f0=f0)
print(f"\n结果已保存: {out}")
