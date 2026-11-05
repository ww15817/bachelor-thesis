"""openEMS 全波计算角锥喇叭馈源（矩形波导口激励，备用方案）。

与 `horn_pinfeed.py` 同一个喇叭，改用 WR-90 矩形波导口直接激励 TE10，
不建销钉和后短路板。保留这个文件的理由是：它的 S11 明显更好
（-21 dB vs -6.9 dB），如果后面要用它做阻抗匹配研究，从这条路走更省事。

已知局限：馈电波导的尾端直接开进 PML，等于波导两端都开口，约有相当
一部分功率向后漏进自由空间（nf2ff 盒把它一起收进去），所以 Dmax 比销钉
版本低约 1.3 dB、口径效率只有 23%。要做成真正的馈源，尾端应该加短路板。

用法同 horn_pinfeed.py：
    set OPENEMS_ROOT=D:\\openEMS
    python horn_feed_wgport.py
"""
import os
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _env  # noqa: F401

import numpy as np

from CSXCAD import ContinuousStructure
from openEMS import openEMS
from openEMS.physical_constants import C0

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

unit = 1e-3
a_wg, b_wg, wg_t = 22.86, 10.16, 1.0
feed_length = 40.0
horn_length = 85.0
ax_ap, ay_ap = 40.0, 30.0

f_start, f_stop, f0 = 8.5e9, 11.5e9, 10e9
edge_deg = 17.761

Sim_Path = os.path.join(HERE, "out_wgport")
os.makedirs(Sim_Path, exist_ok=True)

ang_x = np.degrees(np.arctan((ax_ap - a_wg / 2) / horn_length))
ang_y = np.degrees(np.arctan((ay_ap - b_wg / 2) / horn_length))
len_x = horn_length / np.cos(np.radians(ang_x))
len_y = horn_length / np.cos(np.radians(ang_y))

fc_wg = C0 / (2 * a_wg * unit)
lambda_g = C0 / f0 / unit / np.sqrt(1 - (fc_wg / f0) ** 2)
z_back = -feed_length
z_port0, z_port1 = -30.0, -25.0

print(f"端口 TE10，截面 {a_wg} x {b_wg} mm，z = [{z_port0}, {z_port1}] mm")
print(f"注意：a/b 参数必须以米传入，{a_wg*unit:.5f} / {b_wg*unit:.5f}")

FDTD = openEMS(NrTS=40000, EndCriteria=1e-4)
FDTD.SetGaussExcite(0.5 * (f_start + f_stop), 0.5 * (f_stop - f_start))
FDTD.SetBoundaryCond(["PML_8"] * 6)

CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(unit)

max_res = C0 / f_stop / unit / 20
lam0 = C0 / f0 / unit
margin = lam0 / 2 + 9 * max_res

mesh.AddLine("x", [-ax_ap - margin, -ax_ap, -(a_wg/2 + wg_t), -a_wg/2, a_wg/2,
                   a_wg/2 + wg_t, ax_ap, ax_ap + margin])
mesh.AddLine("y", [-ay_ap - margin, -ay_ap, -(b_wg/2 + wg_t), -b_wg/2, b_wg/2,
                   b_wg/2 + wg_t, ay_ap, ay_ap + margin])
mesh.AddLine("z", [z_back - lam0/4 - 9*max_res, z_back - lam0/4, z_back - wg_t, z_back,
                   z_port0, z_port1, 0, horn_length, horn_length + margin])

horn = CSX.AddMetal("horn")
horn.AddBox(priority=10, start=[-a_wg/2 - wg_t, -b_wg/2, z_back], stop=[-a_wg/2, b_wg/2, 0])
horn.AddBox(priority=10, start=[ a_wg/2, -b_wg/2, z_back], stop=[ a_wg/2 + wg_t, b_wg/2, 0])
horn.AddBox(priority=10, start=[-a_wg/2 - wg_t,  b_wg/2, z_back], stop=[ a_wg/2 + wg_t,  b_wg/2 + wg_t, 0])
horn.AddBox(priority=10, start=[-a_wg/2 - wg_t, -b_wg/2 - wg_t, z_back], stop=[ a_wg/2 + wg_t, -b_wg/2, 0])

# 斜面旋转符号见 horn_pinfeed.py 的注释（两个面的手性相反）
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

port = FDTD.AddRectWaveGuidePort(
    1,
    [-a_wg/2, -b_wg/2, z_port0],
    [ a_wg/2,  b_wg/2, z_port1],
    "z",
    a_wg*unit, b_wg*unit,     # 米！写成 22.86 / 10.16 会被当成 22.86 m 宽
    "TE10", 1.0)

mesh.SmoothMeshLines("all", max_res, 1.4)
print(f"网格 x={mesh.GetQtyLines(0)} y={mesh.GetQtyLines(1)} z={mesh.GetQtyLines(2)}")

nf2ff = FDTD.CreateNF2FFBox(frequency=[f0])

t0 = time.time()
FDTD.Run(Sim_Path, cleanup=True)
print(f"FDTD 用时 {time.time() - t0:.1f} s")

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

E = np.squeeze(res.E_norm[0])
if E.ndim == 1:
    E = E[:, None]
E = np.abs(E) / np.max(np.abs(E))
i_edge = int(np.argmin(np.abs(theta - edge_deg)))

print(f"\nPrad = {Prad:.4e} W (未归一)")
print(f"Dmax = {10*np.log10(Dmax):.2f} dBi")
print(f"S11  @ {f0/1e9:.2f} GHz = {20*np.log10(max(s11[i0],1e-12)):.2f} dB")
for name, col in (("H 面 (phi=0, x-z)", 0), ("E 面 (phi=90, y-z)", 1)):
    e = E[:, col] if E.shape[1] > col else E[:, 0]
    m = (theta > 0) & (theta <= 25) & (e > 0.05)
    k, _ = np.polyfit(np.log10(np.cos(np.radians(theta[m]))), np.log10(e[m]), 1)
    print(f"  {name}: 边缘({edge_deg:.2f}deg) {20*np.log10(max(e[i_edge],1e-9)):+.2f} dB"
          f" | 等效 cos^q 的 q = {2*k:.1f}")

G_ap = 4 * np.pi * (2*ax_ap*unit) * (2*ay_ap*unit) / (C0/f0) ** 2
print(f"理想均匀口径增益 {10*np.log10(G_ap):.2f} dBi -> 口径效率 {100*Dmax/G_ap:.1f}%")

out_dir = os.path.join(ROOT, "data", "feeds")
os.makedirs(out_dir, exist_ok=True)
out = os.path.join(out_dir, "horn_10ghz_wgport.npz")
np.savez(out, theta=theta, phi=phi, E_norm=E, Prad=Prad, Dmax=Dmax,
         s11_freq=freq, s11=s11, f0=f0,
         ax_ap=ax_ap, ay_ap=ay_ap, horn_length=horn_length)
print(f"结果已保存: {out}")
