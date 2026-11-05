# openEMS 馈源（P2）

用全波仿真算喇叭馈源的方向图，替代 GO/PO 流水线里假设的解析 `cos^q`。

## 安装（Windows 原生，不需要 WSL）

以前这条路是在 WSL 里装 openEMS 源码跑的，机器上 WSL 用不了就没法跑。
现在用官方 MSVC 预编译包，解压即可：

1. 下载 `openEMS_x64_v0.37.0_msvc.zip`
   <https://github.com/thliebig/openEMS-Project/releases/download/v0.37.0/openEMS_x64_v0.37.0_msvc.zip>
2. 解压，记下里面 `openEMS` 目录的路径（含 `openEMS.exe`、`CSXCAD.dll`）
3. 装 Python 绑定（zip 内 `openEMS/python/` 下，按 Python 版本选 cp313 / cp314）：

   ```bat
   pip install h5py
   pip install --no-index csxcad-0.7.0-cp313-cp313-win_amd64.whl openems-0.37.0-cp313-cp313-win_amd64.whl
   ```

4. 设环境变量：

   ```bat
   set OPENEMS_ROOT=D:\openEMS
   ```

`_env.py` 会做 `os.add_dll_directory`。这一步不能省：Python 3.8 之后
PATH 不再参与扩展模块的 DLL 搜索，只加 PATH 一定会得到
`ImportError: DLL load failed`。

## 跑

```bat
cd examples\openems
python horn_pinfeed.py                        :: 主模型，80x60 mm 口径
python horn_pinfeed.py --ax 29 --tag a58      :: 换口径做扫描
python horn_feed_wgport.py                    :: 备用方案（波导口激励）
```

单次 10 GHz 算例约 107x91x142 网格、4000 时间步，本机 4~5 分钟。

输出 `data/feeds/horn_10ghz*.npz`：`theta/phi/E_norm/Prad/Dmax/s11/f0`
以及带下划线的几何参数和 `q_fit`（两个面的等效 cos^q 指数）。

## 模型说明

`horn_pinfeed.py` 沿用官方教程 `python/Tutorials/Horn_Antenna.py` 的拓扑：
WR-90 + 后短路板（距销钉 λg/4）+ 同轴销钉（集总端口，50 Ω）激励 TE10，
喇叭四壁用 `AddLinPoly` 斜面旋转而成。改动只有：频率换到 10 GHz、口径
可参数化、后处理加 cos^q 拟合与物理自检。

`horn_feed_wgport.py` 用矩形波导口直接激励。它的 S11 更好（−21 dB），
但波导尾端开口朝自由空间会向后漏功率，Dmax 低约 1.3 dB，所以只作备用。

## 坐标系

`phi=0` 是 x-z 面（H 面，受口径 x 向宽度控制），`phi=90` 是 y-z 面
（E 面，TE10 的 E 沿 y，受 y 向宽度控制）。别按名字直觉判断——反过来
会得出"E 面比 H 面窄"这种物理上说不通的结论。

## 结果（10 GHz，f0 处）

| 模型 | 口径 [mm] | Dmax | 口径效率 | S11 | H 面锥削@17.76° | E 面锥削@17.76° |
| --- | --- | --- | --- | --- | --- | --- |
| 销钉（主） | 80×60 | 13.21 dBi | 31.2% | −6.85 dB | −5.50 dB（q=22.9） | −1.58 dB（q=8.1） |
| 销钉 | 58×60 | 12.02 dBi | 32.7% | −6.58 dB | −2.99 dB（q=11.7） | −1.37 dB（q=6.8） |
| 波导口（备） | 80×60 | 11.92 dBi | 23.2% | −21.22 dB | −3.92 dB（q=17.3） | −7.39 dB（q=33.0） |

流水线的解析馈源假设 `cos^8`，在 17.76° 处锥削 −1.70 dB——正好对上 80×60
销钉模型的 **E 面**（−1.58 dB / q=8.1），但 H 面过锥削 3.8 dB。根因是 x 向
口径 80 mm（2.67λ）相对副面张角太宽；缩到 58 mm 后 H 面回到 −2.99 dB。
这条记为 ISSUE-020，接流水线之前要先定口径。

口径效率 31% 偏低不是数值问题：85 mm 长的喇叭配 80×60 口径，边缘张角
18.6°/16.3°，口面相位偏差已经很大，属于「短而胖」的设计。要提增益得同时
加长喇叭段（`--length`）。

## 踩过的坑（ISSUE-019）

| 现象 | 根因 |
| --- | --- |
| `Prad ≈ 5e-38 W`，方向图是噪点 | `AddRectWaveGuidePort` 的 `a/b` 必须以**米**传，原来传了 mm 数值 22.86/10.16（= 22.86 m 宽的模式，完全不匹配）。官方教程写的是 `a*unit` |
| `S11 ≈ 0 dB`，能量困在腔内缓慢衰减，`Prad≈0` | 喇叭四壁旋转符号写反，喇叭朝内收成封闭角锥。绕 x 轴时 −y 壁取 +ang、+y 壁取 −ang；绕 y 轴时 −x 壁取 −ang、+x 壁取 +ang。两个面手性相反，不能共用一个符号 |
| 两组斜面收不到同一个口径面 | 两个面张角不同时，斜面在旋转前要按 `horn_length/cos(ang)` 加长 |

`horn_pinfeed.py` 末尾有三条物理自检（S11、Dmax、口径效率），数值不对会
直接告警，不用等方向图出来才发现。
