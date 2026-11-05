# 天线仿真开源软件集成与应用

全开源电磁-结构双域仿真链路（毕业设计）。统一算例：20λ 卡塞格伦天线。

## 链路总览

| 域 | 工具 | 职责 |
|---|---|---|
| 射频-全波 | openEMS (FDTD) | 馈源建模、近场求解 |
| 射频-高频近似 | 自研 GO/PO 积分 | 主/副反射面电大尺寸散射与远场 |
| 结构 | FEniCSx (线弹性有限元) | 反射面面板变形与表面误差 |

核心矛盾：电大尺寸天线直接做全波仿真，网格量爆炸、内存扛不住。
对策：馈源仍用全波（网格量小），反射面改用 GO/PO 高频近似（不需要体网格），
结构域用有限元单独求解，两侧通过「表面误差 -> 相位扰动」耦合。

## 目录结构

```
docs/                          方案、技术路线、环境搭建文档
scripts/                       环境安装与一键脚本
src/
  common/                      公共常量、坐标与单位约定
  go_po/                       自研几何光学 / 物理光学积分
  fenics/                      FEniCSx 结构有限元与面板变形（当前只有环境自检）
examples/cassegrain_20lambda/  统一算例
data/results/                  运行时生成的结果（不入库）
work/                          进度台账、审查报告、扫描脚本
```

> `src/openems/`（openEMS 馈源建模）与 `notebooks/` 目前**尚未提交**，属计划中部分。


## 快速开始

```powershell
# 1. Windows 侧 Python 环境（无需管理员）
.\scripts\setup-python-env.ps1

# 2. WSL + openEMS + FEniCSx（需管理员 PowerShell）
.\scripts\install-wsl-openems.ps1
```

详细步骤见 `docs/02-环境搭建.md`。

## 当前进度

> **先看 [`STATUS.md`](STATUS.md)** —— 那里有目标、已完成情况、是否理想的判断、以及下一步做什么。
> 下面是里程碑层面的快照，细节以 `STATUS.md` 和 [`work/进度台账.md`](work/进度台账.md) 为准。

| 里程碑 | 状态 |
|---|---|
| M1 环境搭建 | **已完成**（Windows + WSL 双层，openEMS v0.37.0 + FEniCSx 0.11.0） |
| M2 抛物面 GO/PO 验证 | 已完成，见下 |
| M3 卡塞格伦双反射面 | **已完成**，见下 |
| M4 结构耦合 | 耦合链路通了（解析挠度版），FEniCSx 真实变形场未接 |
| M5 论文 | 未开始 |

> 回来要做的步骤、以及我改过的问题，都写在 [`HANDOFF.md`](HANDOFF.md)。
>
> **怎么用这些工具，看 [`docs/03-工具使用手册.md`](docs/03-工具使用手册.md)。**
> 环境怎么搭的、踩过哪些坑，见 [`docs/02-环境搭建.md`](docs/02-环境搭建.md)。

### M2 验证结果

运行：

```powershell
.\.venv\Scripts\python.exe examples\cassegrain_20lambda\validate_dish.py
```

结果（20λ 抛物面，F/D = 0.4）：

- 均匀口径的数值积分与解析 Airy 图案最大偏差 **5.8e-06**，第一零点 3.5°；
- 馈源加权口径的方向图正常（-10 dB 边缘照射）；
- 面元 PO 积分可跑，零点位置与口径法一致。

已知待办：深抛物面（张角 64°）下，口径场法省略了斜入射投影因子，
导致两条曲线的副瓣电平有差异。这个差异需要量化和解释，是论文的一个点。

### M3 验证结果

运行：

```powershell
.\.venv\Scripts\python.exe examples\cassegrain_20lambda\m3_cassegrain.py
```

结果（D = 20λ，F/D = 0.4，M = 4）：

- 射线 **1872/1872** 走完「馈源 -> 副面 -> 主面」全链路；
- 准直误差 **4.4e-16**、光程起伏 **7.1e-15**（都是机器精度）；
- 副面口径 4.433λ，遮挡 **4.913%**，馈源位于主面顶点 z = 0；
- 遮挡造成的口径场幅损失 **-0.48 dB**（场幅比，**不是**增益；功率等效约 −0.24 dB）（2026-11-03 修正口径幅度模型后；旧值 -0.525 dB 基于不守恒的 g/(r1+r2)）；

产出：`data\results\m3_cassegrain.png`（四联图：几何+射线 / 口径幅度 / 方向图 / 遮挡-增益）。

## 本机环境现状

自检：`powershell -NoProfile -ExecutionPolicy Bypass -File scripts\verify-env.ps1`
→ **20 项通过，1 项失败**（失败项是 MSVC，已由 WSL 替代，不影响）。

**Windows 侧**

- Python 3.11.9，项目虚拟环境 `.venv`（numpy/scipy/matplotlib/gmsh/pyvista/meshio/h5py/pandas）
- Qt 6.8.3 + GCC 13.1.0 + CMake 4.x（`C:\Qt`），`work\qt_probe\qt_probe.exe` 已编译

**WSL 侧（Ubuntu 22.04.5）**

- **openEMS v0.37.0**（`~/openems-install`），链接 CSXCAD v0.7.0
- **FEniCSx dolfinx 0.11.0 + petsc4py 3.25.6**（`~/miniforge3/envs/fenicsx`）

```powershell
wsl -d Ubuntu-22.04
```
