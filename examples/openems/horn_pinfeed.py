"""openEMS 全波计算角锥喇叭馈源（销钉馈电，10 GHz，WR-90）。

模型沿用官方教程 `python/Tutorials/Horn_Antenna.py`（openEMS v0.37）的拓扑：
WR-90 矩形波导 + 后短路板 + 同轴销钉（集总端口）激励 TE10，喇叭四壁用
`AddLinPoly` 斜面 + 旋转得到。改动只有几何尺寸（10 GHz、口径可调）和
后处理（拟合等效 cos^q、导出供 GO/PO 流水线用的方向图）。

为什么用销钉而不用矩形波导口：ISSUE-019 里波导口激励不出来，根因是
`AddRectWaveGuidePort` 的 a/b 必须以**米**传入（官方教程写 a*unit），
原脚本传了 mm 数值。波导口版本也修好了，见 `horn_feed_wgport.py`，但
它的馈电端开口朝自由空间，会向后泄漏约一半功率，Dmax 偏低；作为馈源
模型应当有后短路板，所以主模型用本文件。

用法：
    set OPENEMS_ROOT=D:\\openEMS
    python horn_pinfeed.py                       # 默认 80x60 mm 口径
    python horn_pinfeed.py --ax 29 --tag a58     # 换口径做扫描

输出：
    data/feeds/horn_10ghz[_tag].npz  theta/phi/E_norm/Prad/Dmax/s11/f0
"""
import argparse
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _env  # noqa: F401  DLL 引导，必须在 import CSXCAD 之前

import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

# ---------------------------------------------------------------- 命令行
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
ap.add_argument("--ax", type=float, default=40.0, help="口径 x 向半宽 [mm]")
ap.add_argument("--ay", type=float, default=30.0, help="口径 y 向半宽 [mm]")
ap.add_argument("--length", type=float, default=85.0, help="喇叭段长度 [mm]")
ap.add_argument("--edge-deg", type=float, default=17.761,
                help="关注的下反射面边缘角 [deg]，用于报边缘锥削")
ap.add_argument("--tag", default="", help="输出文件名后缀")
ap.add_argument("--post-only", action="store_true", help="跳过 FDTD，只做后处理")
args = ap.parse_args()

# ---------------------------------------------------------------- 参数
unit = 1e-3                     # 全部长度以 mm 计
a_wg, b_wg, wg_t = 22.86, 10.16, 1.0
feed_length = 40.0
horn_length = args.length
ax_ap, ay_ap = args.ax, args.ay

f_start, f_stop, f0 = 8.5e9, 11.5e9, 10e9

tag = args.tag or f"a{2*ax_ap:.0f}b{2*ay_ap:.0f}"
Sim_Path = os.path.join(HERE, "out_" + tag)
os.makedirs(Sim_Path, exist_ok=True)

# 由口径反算两个面的张角。两个面张角不同，斜面在旋转前要比 horn_length 长
# 一点（除以 cos），否则旋转投影后两组斜面收不到同一个口径面。
ang_x = np.degrees(np.arctan((ax_ap - a_wg / 2) / horn_length))
ang_y = np.degrees(np.arctan((ay_ap - b_wg / 2) / horn_length))
len_x = horn_length / np.cos(np.radians(ang_x))
len_y = horn_length / np.cos(np.radians(ang_y))

fc_wg = C0 / (2 * a_wg * unit)
lambda_g = C0 / f0 / unit / np.sqrt(1 - (fc_wg / f0) ** 2)
back_short = lambda_g / 4          # 后短路板到销钉 = λg/4
pin_length = 0.55 * b_wg           # 销钉长度，教程给的经验起点
pin_r, port_h, feed_R = 0.5, 1.0, 50.0

z_back = -feed_length
z_pin = z_back + back_short

A_aperture = (2 * ax_ap * unit) * (2 * ay_ap * unit)

print(f"口径 {2*ax_ap:.1f} x {2*ay_ap:.1f} mm "
      f"({2*ax_ap*unit*f0/C0:.2f} x {2*ay_ap*unit*f0/C0:.2f} lambda), "
      f"喇叭长 {horn_length} mm")
print(f"张角 x {ang_x:.2f} deg / y {ang_y:.2f} deg | TE10 截止 {fc_wg/1e9:.2f} GHz "
      f"| λg {lambda_g:.2f} mm | 后短路 {back_short:.2f} mm")

# ---------------------------------------------------------------- FDTD
FDTD = openEMS(NrTS=40000, EndCriteria=1e-4)
FDTD.SetGaussExcite(0.5 * (f_start + f_stop), 0.5 * (f_stop - f_start))
FDTD.SetBoundaryCond(["PML_8"] * 6)

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

max_res = C0 / f_stop / unit / 20          # 1.25 mm ≈ λ/20 @12 GHz
lam0 = C0 / f0 / unit
margin = lam0 / 2 + 9 * max_res

x_max, y_max = ax_ap + margin, ay_ap + margin
z_sim_back = lam0 / 4 + 9 * max_res
z_sim_front = margin

mesh.AddLine("x", [-x_max, -ax_ap, -(a_wg / 2 + wg_t), -a_wg / 2, -pin_r, pin_r,
                   a_wg / 2, a_wg / 2 + wg_t, ax_ap, x_max])
mesh.AddLine("y", [-y_max, -ay_ap, -(b_wg / 2 + wg_t), -b_wg / 2,
                   -b_wg / 2 + port_h, -b_wg / 2 + pin_length,
                   b_wg / 2, b_wg / 2 + wg_t, ay_ap, y_max])
mesh.AddLine("z", [z_back - z_sim_back, z_back - lam0 / 4, z_back - wg_t, z_back,
                   z_pin - pin_r, z_pin + pin_r,
                   0, horn_length, horn_length + z_sim_front])

# ---------------------------------------------------------------- 金属
horn = CSX.AddMetal("horn")

# 馈电波导四壁 + 后短路板
horn.AddBox(priority=10, start=[-a_wg/2 - wg_t, -b_wg/2, z_back], stop=[-a_wg/2, b_wg/2, 0])
horn.AddBox(priority=10, start=[ a_wg/2, -b_wg/2, z_back], stop=[ a_wg/2 + wg_t, b_wg/2, 0])
horn.AddBox(priority=10, start=[-a_wg/2 - wg_t,  b_wg/2, z_back], stop=[ a_wg/2 + wg_t,  b_wg/2 + wg_t, 0])
horn.AddBox(priority=10, start=[-a_wg/2 - wg_t, -b_wg/2 - wg_t, z_back], stop=[ a_wg/2 + wg_t, -b_wg/2, 0])
horn.AddBox(priority=10,
            start=[-a_wg/2 - wg_t, -b_wg/2 - wg_t, z_back - wg_t],
            stop= [ a_wg/2 + wg_t,  b_wg/2 + wg_t, z_back])

# 四片斜面。旋转符号不能省：绕 x 轴时 -y 壁取 +ang、+y 壁取 -ang；
# 绕 y 轴时 -x 壁取 -ang、+x 壁取 +ang。符号写反 = 喇叭朝内收成封闭
# 角锥，表现为 S11≈0 dB、Prad≈0、能量困在腔内缓慢衰减。
z_tb = np.array([0, len_y, len_y, 0])
x_tb = np.array([a_wg/2, ax_ap, -ax_ap, -a_wg/2])

p = horn.AddLinPoly(points=[z_tb, x_tb], norm_dir="y",
                    elevation=-wg_t/2, length=wg_t, priority=10)
p.AddTransform("RotateAxis", "x", ang_y)
p.AddTransform("Translate", [0, -b_wg/2 - wg_t/2, 0])

p = horn.AddLinPoly(points=[z_tb, x_tb], norm_dir="y",
                    elevation=-wg_t/2, length=wg_t, priority=10)
p.AddTransform("RotateAxis", "x", -ang_y)
p.AddTransform("Translate", [0, b_wg/2 + wg_t/2, 0])

y_lr = np.array([b_wg/2 + wg_t, ay_ap + wg_t, -ay_ap - wg_t, -b_wg/2 - wg_t])
z_lr = np.array([0, len_x, len_x, 0])

p = horn.AddLinPoly(points=[y_lr, z_lr], norm_dir="x",
                    elevation=-wg_t/2, length=wg_t, priority=10)
p.AddTransform("RotateAxis", "y", -ang_x)
p.AddTransform("Translate", [-a_wg/2 - wg_t/2, 0, 0])

p = horn.AddLinPoly(points=[y_lr, z_lr], norm_dir="x",
                    elevation=-wg_t/2, length=wg_t, priority=10)
p.AddTransform("RotateAxis", "y", ang_x)
p.AddTransform("Translate", [a_wg/2 + wg_t/2, 0, 0])

# ---------------------------------------------------------------- 销钉馈电
port = FDTD.AddLumpedPort(
    1, feed_R,
    [-pin_r, -b_wg/2, z_pin - pin_r],
    [pin_r, -b_wg/2 + port_h, z_pin + pin_r],
    "y", 1.0, priority=5)

pin = CSX.AddMetal("pin")
pin.AddBox(priority=11,
           start=[-pin_r, -b_wg/2 + port_h, z_pin - pin_r],
           stop=[pin_r, -b_wg/2 + pin_length, z_pin + pin_r])

# ---------------------------------------------------------------- 收尾
mesh.SmoothMeshLines("all", max_res, 1.4)
print(f"网格 x={mesh.GetQtyLines(0)} y={mesh.GetQtyLines(1)} z={mesh.GetQtyLines(2)}")

nf2ff = FDTD.CreateNF2FFBox(frequency=[f0])

if not args.post_only:
    t0 = time.time()
    FDTD.Run(Sim_Path, cleanup=True)
    print(f"FDTD 用时 {time.time() - t0:.1f} s")

# ---------------------------------------------------------------- 后处理
theta = np.arange(0.0, 91.0, 1.0)
phi = np.array([0.0, 90.0])
res = nf2ff.CalcNF2FF(Sim_Path, f0, theta, phi,
                      outfile=os.path.join(Sim_Path, "nf2ff.h5"))

Prad = float(np.squeeze(res.Prad[0]))
Dmax = float(np.squeeze(res.Dmax[0]))

freq = np.linspace(f_start, f_stop, 201)
port.CalcPort(Sim_Path, freq)
s11 = np.abs(port.uf_ref / port.uf_inc)
i0 = int(np.argmin(np.abs(freq - f0)))
s11_dB = 20 * np.log10(max(s11[i0], 1e-12))

E = np.squeeze(res.E_norm[0])
if E.ndim == 1:
    E = E[:, None]
E = np.abs(E) / np.max(np.abs(E))

# 坐标系：TE10 的 E 沿 y，所以 phi=90（y-z 面）是 E 面，phi=0（x-z 面）
# 是 H 面。别按名字猜，反过来会得到"E 面比 H 面窄"的物理上说不通的结论。
planes = (("H 面 (phi=0, x-z)", 0), ("E 面 (phi=90, y-z)", 1))
i_edge = int(np.argmin(np.abs(theta - args.edge_deg)))

print(f"\nPrad = {Prad:.4e} W   （未按馈入功率归一，绝对值无意义，只看方向图/Dmax）")
print(f"Dmax = {Dmax:.4f} = {10*np.log10(Dmax):.2f} dBi")
print(f"S11  @ {f0/1e9:.2f} GHz = {s11_dB:.2f} dB")

q_fit = {}
for name, col in planes:
    e = E[:, col] if E.shape[1] > col else E[:, 0]
    m = (theta > 0) & (theta <= 25) & (e > 0.05)
    if m.sum() > 3:
        k, _ = np.polyfit(np.log10(np.cos(np.radians(theta[m]))), np.log10(e[m]), 1)
        q_fit[name] = 2.0 * k
    print(f"  {name}: 边缘({args.edge_deg:.2f}deg) {20*np.log10(max(e[i_edge],1e-9)):+.2f} dB"
          + (f" | 等效 cos^q 的 q = {q_fit[name]:.1f}" if name in q_fit else ""))

G_ap = 4 * np.pi * A_aperture / (C0 / f0) ** 2
print(f"理想均匀口径增益 {10*np.log10(G_ap):.2f} dBi -> 口径效率 {100*Dmax/G_ap:.1f}%")

# 物理自检：这三个数不对就说明模型或激励有问题，别急着用方向图
if s11_dB > -8:
    print("! 警告：S11 偏高，销钉长度/后短路距离没调好，反射会污染方向图主瓣附近")
if 10 * np.log10(max(Dmax, 1e-12)) < 10:
    print("! 警告：Dmax 过低，检查喇叭口是否被封死（斜面旋转符号）")

out_dir = os.path.join(ROOT, "data", "feeds")
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, f"horn_10ghz{'_' + args.tag if args.tag else ''}.npz")
np.savez(out, theta=theta, phi=phi, E_norm=E, Prad=Prad, Dmax=Dmax,
         s11_freq=freq, s11=s11, f0=f0,
         ax_ap=ax_ap, ay_ap=ay_ap, horn_length=horn_length,
         q_fit=np.array([q_fit.get(k, np.nan) for k, _ in planes]))
print(f"结果已保存: {out}")
