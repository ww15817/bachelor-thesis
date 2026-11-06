# 图表清单

论文用的图放这儿。每张图注明：对应哪一节、哪条命令生成、图里该看什么。

## 已有

| 文件 | 对应章节 | 生成方式 | 看什么 |
|---|---|---|---|
| `m2_validation.png` | 第 4 章·Airy 验证 | `examples/cassegrain_20lambda/validate_dish.py` | 数值 PO 与解析 Airy 重合到 7e-06 |
| `m3_cassegrain.png` | 第 2 章·几何 / 第 4 章·遮挡 | `examples/cassegrain_20lambda/m3_cassegrain.py` | 四联：几何+射线 / 口径幅度 / 方向图 / 遮挡-增益 |
| `cassegrain_20lambda.png` | 第 4 章·主算例方向图 | `run_pipeline.py cases/cassegrain_20lambda.yaml` | 主瓣、第一副瓣、遮挡台阶 |
| `budget_waterfall.png` | 第 4 章·效率预算 | `work/make_figures.py` | 理想 35.96 → 溢出 −4.49 → 锥削+遮挡 −0.22 → 31.24 dBi |
| `q_scan.png` | 第 5 章·馈源优化 | `work/make_figures.py` | 最优 q=47 是平台（47–64），不是单点 |
| `chain_diagram.png` | 第 3 章·集成架构 | `work/make_figures.py` | 三块 + U 形耦合回路；下方注明同一条命令跑完、每段的验证数 |
| `deflection_map.png` | 第 6 章·变形场 | `work/make_figures.py` | 左：FEM 变形 w（0–543 μm）；右：口径相位 4πw/λ（0–0.225 rad） |
| `spill_schematic.png` | 第 5 章·溢出 | `work/make_figures.py` | 按 cos⁸θ 功率密度画，**面积即功率占比**：35.6% 打到副面、64% 溢出 |
| `pattern_compare.png` | 第 6 章·耦合影响 | `work/make_figures.py` | **变形 ×10** 才有可见差别：主瓣变宽、零点被填、副瓣从 −18 抬到 −14 dB。真实 2 mm 板只掉 0.014 dB，画不出差别（那本身是结论） |
| `surface_error.png` | 第 6 章·灵敏度 | `work/make_figures.py` | 随机相位屏下 PO 与 Ruze 一致；确定性变形另有 ISSUE-018 |
| `cassegrain_fd05.png` | 第 4 章·F/D 对照 | `run_pipeline.py cases/cassegrain_fd05.yaml` | fd05 比 20λ 差 1.59 dB |

## 还缺的（按论文顺序）

1. **链路总框图**——openEMS → GO/PO → FEniCSx 三块加耦合箭头。svg 或 ppt 画，不必用代码
2. **效率预算瀑布图**——理想 35.96 → 溢出 −4.49 → 锥削/遮挡 −0.22 → 31.24 dBi。数据在台账 W-0015
3. **q 扫描曲线**——横轴 q、纵轴绝对增益，标出 47 和平台的宽度。数据在 `work/scan_q.py`
4. **溢出发射示意**——馈源张角 vs 副面张角（17.76°）、回来那 64% 去哪了
5. **表面误差灵敏度**——σ_w/λ 对增益损失，加 Ruze 对照。数据在 `work/scan_surface_error.py`
6. **变形场云图**——FEniCSx 的 w(x,y) 或口径相位分布。数据在 `data/feeds/plate_deflection_2mm.npz`
7. **误差账表**（不是图，是表）——每个近似的偏差、方向、来源

## 规矩

- 图只放 `docs/figures/`，别放 `outputs/`（那里是运行产物，会进 `.gitignore`）
- 每张图都要能被一条命令重新生成，命令写在上面这张表里
- 图里不要出现本机路径

## 每张图的假设（写论文时图下要带一句）

| 图 | 基于什么假设 | 哪些部分是示意 |
|---|---|---|
| `m2_validation.png` | 均匀口径、无遮挡 | 无 |
| `m3_cassegrain.png` | D=20λ, F/D=0.4, M=4；副面按几何阴影处理 | 无 |
| `cassegrain_20lambda.png` | 同上 + q=8 解析馈源 | 无 |
| `cassegrain_fd05.png` | 同上，F/D=0.5 | 无 |
| `budget_waterfall.png` | 20λ、q=8；溢出按解析式、口径效率按未遮挡区积分 | 无（柱高即数值） |
| `q_scan.png` | 同上，q 从 4 扫到 96；平台判据 0.10 dB | 无 |
| `surface_error.png` | **随机**相位屏（相关长度远小于口径），12 次实现平均 | 无（确定性变形不适用，见 ISSUE-018） |
| `deflection_map.png` | 方板、简支边施 M=0（偏低 5.4%）、80×80 网格 | 云图是插值，边界一圈有锯齿 |
| `spill_schematic.png` | cos⁸θ 功率密度；副面 17.76°、主面边缘 72.65° | 二维示意图，角度轴不代表空间位置 |
| `pattern_compare.png` | 同 deflection_map，但**变形放大 10 倍**才看得出 | 倍率是人为的，真实值只掉 0.014 dB |
| `chain_diagram.png` | — | 纯示意，方框大小不代表任何比例 |

## 三条不能含糊的

1. **数据图**（前七张）的每个数字都能在仓库里重算；跑法在上一张表里
2. **示意图**（最后两张）不能用来读数字，只能用来讲结构
3. 所有图忠实的是**模型**，不是真实天线。模型本身有已知偏差（板弯曲 5.4%、Ruze 适用条件），图不会藏这些，但也修不掉

