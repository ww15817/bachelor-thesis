# 交接说明

更新时间：2026-10-15 深夜

> **本文是 2026-10-15 的历史记录。** 当前进度、未决问题与下一步都见
> `work/进度台账.md`；本文的「已完成」结论（尤其是环境与结构耦合）已过期。

---

## 一、先说结论：环境全通了

你洗澡这段时间，我把之前卡住的全干完了：

| 组件 | 版本 | 状态 | 验证方式 |
|---|---|---|---|
| WSL | 3.0.3.0 | ✅ | `wsl --version` |
| Ubuntu | 22.04.5 LTS | ✅ | `wsl -l -v` |
| **openEMS** | **v0.37.0** | ✅ | `openEMS --version`，链接 CSXCAD v0.7.0 / hdf5 / VTK 9.1 |
| **FEniCSx** | **dolfinx 0.11.0** | ✅ | 二维 Poisson 求解，误差 1.47e-05 |
| petsc4py | 3.25.6 | ✅ | Windows 上死活装不了的那个，Linux 上一次过 |
| Qt 6.8.3 + MinGW 13.1 | — | ✅ | `work\qt_probe\qt_probe.exe` |
| Windows Python `.venv` | — | ✅ | 13 个科学计算包全可用 |

**你现在什么都不用装了。** 直接进第四节开始干活。

---

## 二、怎么用这套环境

### 开一个 WSL 终端

```powershell
wsl -d Ubuntu-22.04
```

进去之后（环境变量已写进 `~/.bashrc`，新终端自动生效）：

```bash
# openEMS
openEMS --version
which openEMS          # ~/openems-install/bin/openEMS

# FEniCSx
conda activate fenicsx
python -c "import dolfinx; print(dolfinx.__version__)"
```

### 跑项目脚本

项目目录在 WSL 里挂载于 `/mnt/c/...`：

```bash
cd /mnt/c/Users/ww/Documents/Codex/2026-10-09/ai

# GO/PO 部分（Windows 的 .venv 就能跑）
/mnt/c/Users/ww/Documents/Codex/2026-10-09/ai/.venv/Scripts/python.exe examples/cassegrain_20lambda/m3_cassegrain.py

# FEniCSx 部分（要用 conda 环境）
conda activate fenicsx
python src/fenics/check_env.py
```

### 在哪里看文件

| 内容 | Windows 路径 | WSL 路径 |
|---|---|---|
| 项目 | `C:\Users\ww\Documents\Codex\2026-10-09\ai` | `/mnt/c/Users/ww/Documents/Codex/2026-10-09/ai` |
| openEMS 安装 | 不可直接访问 | `~/openems-install` |
| FEniCSx 环境 | 不可直接访问 | `~/miniforge3/envs/fenicsx` |
| openEMS 源码 | 不可直接访问 | `~/openems-src` |

---

## 三、可选：装 AppCSXCAD 图形界面

openEMS 的 **核心计算功能已经全部可用**，唯一没装的是 GUI（AppCSXCAD），
因为它要 `libvtk9-qt-dev`，那个需要 sudo。

如果你想要图形界面看模型（不装也行，我们有 pyvista 画图）：

```bash
sudo apt-get install -y libvtk9-qt-dev
cd ~/openems-src && ./update_openEMS.sh --njobs=16 ~/openems-install
```

要输入 Linux 密码（不显示字符）。

---

## 四、接下来该干的（论文正题）

环境这一步彻底结束了，下面是真正的毕设内容。建议顺序：

> **先把 [`docs/03-工具使用手册.md`](docs/03-工具使用手册.md) 读一遍。**
> 里面有 openEMS / FEniCSx / 自研代码三者的工作流、API 速查、
> 最小可运行例子，以及我踩过的所有真实报错和解决办法。

### 4.1 openEMS 馈源建模

目标：建一个角锥喇叭，导出近场/口径场，喂给 GO/PO 链路。

**我已经实测跑通了官方角锥喇叭教程**（15 GHz，241 万网格，5 分半）：

```
Directivity:          Dmax = 16.9 dBi
Aperture efficiency:  e_a  = 64.6 %
```

结果和输出位置（`/tmp/Horn_Antenna_PinFeed/nf2ff.h5`）都写在
`docs/03-工具使用手册.md` 的 2.4 节。

openEMS 自带教程可以参考：

```bash
mkdir -p ~/openems-tutorials && cd ~/openems-tutorials
cp ~/openems-install/share/openEMS/python/Tutorials/Horn_Antenna.py .
```

产出：`Horn_Antenna.py` 改造成我们算例的馈源，导出表面场数据。

### 4.2 FEniCSx 结构域

目标：反射面面板在自重/风载下的变形 w(x, y)。

先从最小验证开始：均布载荷圆板，与解析解对比。

### 4.3 结构-电磁耦合闭环

把变形 w 转成口径相位扰动 `dphi = 2k w`，回代到 GO/PO 的相位项。

**判据已经有了**：`scripts/validate_po.py` 里的 L3 已经把 Ruze 公式
`G/G0 = exp(-sigma_phi^2)` 验证到 0.001 dB 偏差，缺的只是「真实变形场」。

> 更正（2026-11-04，审查 C-07）：此行原文写的是 `exp(-sigma^2/2)`，与代码不一致。
> 标准 Ruze 是 `exp(-σφ²)`（σφ 为**相位**误差 rms）。写成 `σ²/2` 在同样变形下
> 会少算一倍，是会被评委抓的错。

### 4.4 代码层面的已知待改进点

| 项 | 说明 | 优先级 |
|---|---|---|
| 口径幅度模型 | 现在是 `g(theta)/(r1+r2)` 一阶模型，严格做法要跟踪射线管 Jacobian | 高 |
| 绕射效应 | 目前纯 PO，没有 PTD/UTD，远角副瓣会偏乐观 | 中 |
| `outputs\validate_po.png` 右下角 | 「口径场相位分布」量程设成 ±0.1 rad，图基本全黑，建议调量程或换图 | 低 |
| `examples/openems/quickstart_dipole.py` | 工具链冒烟测试，网格偏粗、物理结果不作数，只用来看链路通不通 | — |

---

## 五、我修过的问题（写论文/答辩可用）

### 卡塞格伦几何（最严重）

`src/common/cassegrain.py` 原来 `z_center = L + a` 是错的，
导致**近焦点算在 z = -6 而不是 +8，偏 14λ**。它自带的
`assert_consistency()` 把它自己判死了，所以这套代码从来没被任何脚本驱动过。

正确关系：

```
c - a = F - L
a = (M-1)(F-L)/2      c = (M+1)(F-L)/2      z_center = L - a
=> 近焦点 = z_center + c = F   （恒成立）
```

另外 **`L = F/M` 不能用来选副面位置**：代进去得到 16.4λ 的副面、
**67% 遮挡**，物理上不可用。`L` 是自由设计参数，控制遮挡比；
默认取 `L = M·F/(M+1)`（馈源落在主面顶点），遮挡约 4.9%。

### 射线追踪

原版把馈源放在抛物面焦点 z=F，按「馈源 -> 主面 -> 副面」追踪，**路径反了**。
正确是「馈源（双曲面远焦点）-> 副面 -> 主面」。

### 光程不变量

原版用 `r1 + r2` 当光程，它本来就不是常数。正确的是 **`r1 + r2 - z_main`**
（出射平行于 z，主面到参考面的距离为 Z0 - z_main）。
改后起伏从 2.93 降到 7.1e-15。

### PO 积分归一化

`far_field_aperture` 原来用被照射函数自己的峰值归一化，导致有/无遮挡
两条曲线用了不同参考，**遮挡增益损失恒为 0 dB**。
加 `normalize=` 参数共用参考后，算出真实的 -0.525 dB。

### PowerShell 脚本

1. 四个 `.ps1` 都没有 UTF-8 BOM。Windows PowerShell 5.1 不把无 BOM 的 UTF-8
   当 UTF-8 读，中文字符会截断字符串引号，脚本直接解析失败。**已全部加 BOM。**
2. `verify-env.ps1` 里 PATH 在 dolfinx 导入**之前**就还原，误报「导入失败」。
3. 同脚本里 `Check 'UFL 版本匹配' { $false }` 是写死的 TODO 占位。
4. `install-wsl-openems.ps1` 最后一行用 `\"` 转义（PowerShell 该用反引号），
   那行被当命令执行了。
5. `fetch_wsl_parallel.ps1` 里 `$parts` 与参数 `$Parts` 冲突
   （**PowerShell 变量大小写不敏感**），且 `-H "Range: ..."` 的空格被拆开
   导致分块全是 0 字节。已改名并改用 `-r`。

### FEniCSx 为什么必须走 WSL

原生 Windows 上两条死结，我都实测确认过：

- `petsc4py` 在 conda-forge 的 **win-64 上不存在**，PyPI 也**没有 Windows 轮子**；
- ffcx 用 CFFI 做表单 JIT，**Windows 上只认 MSVC（cl.exe）**，机器上没装 Visual Studio。

Linux 上这两条都不存在——`petsc4py 3.25.6` 和整个 dolfinx 一次装好。

### WSL 包下载

GitHub 直连只有 25 KB/s（350MB 的包要下 8 小时）。换成 gh-proxy 镜像
**8 连接并行分块**，跑到 4 MB/s，两分半下完。
下载的 MSI 我做过完整校验：大小、SHA256 与 GitHub 官方 digest 一致，
微软数字签名 Valid。

---

## 六、常用命令备忘

```powershell
# WSL 终端
wsl -d Ubuntu-22.04

# 重跑环境自检（Windows 侧）
cd C:\Users\ww\Documents\Codex\2026-10-09\ai
powershell -NoProfile -ExecutionPolicy Bypass -File scripts\verify-env.ps1
```

```bash
# 在 WSL 里
cd /mnt/c/Users/ww/Documents/Codex/2026-10-09/ai

openEMS --version
conda activate fenicsx && python -c "import dolfinx; print(dolfinx.__version__)"

# 一键重装（幂等，已装的会跳过）
bash scripts/bootstrap-ubuntu.sh
```

---

## 七、可以清理的东西（约 2.6 GB）

安装包和临时文件都没用了。我这个会话的删除操作被策略拦住了，你自己删：

```bash
# WSL 里
rm -f ~/miniforge-installer.sh
rm -rf ~/openems-src/*/build        # 编译中间产物
```

```powershell
# Windows 里
Remove-Item 'C:\Users\ww\Documents\Codex\2026-10-09\ai\work\wsl-pkg\wsl.3.0.3.0.x64.msi.part*'
Remove-Item 'C:\Users\ww\Documents\Codex\2026-10-09\ai\work\wsl-pkg\Microsoft.WSL_3.0.3.0_x64_ARM64.msixbundle'
Remove-Item 'C:\Users\ww\Documents\Codex\2026-10-09\ai\work\wsl-pkg\speedtest_*.bin'
Remove-Item 'C:\Users\ww\Documents\Codex\2026-10-09\ai\work\wsl-pkg\rangetest.bin'
```

**`wsl.3.0.3.0.x64.msi` 已经装完了，也可以删。**
