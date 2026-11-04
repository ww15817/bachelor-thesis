# 代码审查 2026-11-04（WorkBuddy）

审查对象：`2d61133`（本轮含 `f4cd0fb` / `94d5af7` / `2d61133` 三个提交）
方法：独立 venv（numpy 2.x + scipy）实跑 `pytest`、两个算例、`scripts/validate_po.py`；关键数字用独立脚本复算，不采信代码里的自述。
约定：**只写本文件，不改 `src/`**；push 前先 `git pull --rebase`。

---

## 本轮已确认修好（4 项，不必再动）

| 原问题 | 现状 | 我的复测 |
|---|---|---|
| `far_field_aperture` 一次性构造 `(n_theta,n_rho,n_phi)` | 已改逐 θ 循环 | 流水线默认规模峰值 **2673 MB → 43 MB**，耗时 3.75→3.76 s；与旧实现 `max\|diff\| = 0.0`（逐位相同） |
| `work/dsh-session/` 2.4 MB 会话日志入库 | 已 `git rm` | `git ls-files` 中已无该目录 |
| 溢出效率完全不在账上 | 已打印诊断 | `0.3557 (-4.49 dB)`，与我独立数值积分 `0.355666` 一致；`1-cos^(q+1)(θe)` 无误 |
| F/D=0.5 对照算例结论方向反了 | 已补绝对增益 | fd05 `29.65 dBi` vs 20λ `31.24 dBi`，**现在正确体现 fd05 更差 1.59 dB** |

顺带核过、都没问题：`feed_half_angle_for_edge` 给的 17.7613° 与「馈源→副面边缘」几何角差 0.000 角分（用这个角做溢出积分是对的）；`57 passed`；L1/L2/L3 全过（L1 偏差 1.4e-5，L3 与 Ruze 差 0.001 dB）。

---

## C-01 | 高 | src/go_po/aperture_po.py:165

现象：`far_field_pattern` 里 `np.exp(1j * np.outer(kx, af.x))` 仍是 `(n_theta × n_points)` 一次成形的大复数组，内存照旧爆。

实测：
```bash
# uniform_aperture(10.0) 默认 12/λ → 57600 点；THETA 4801 点
# phase 数组本身 = 4801 × 57600 × 16 B = 4.42 GB
当前实现：10.22 s，进程峰值 8888 MB
逐 θ 循环： 9.73 s，进程峰值   41 MB，max|diff| = 6.7e-16
```
`python scripts/validate_po.py` 全程峰值 **9993 MB**、138 s。

判据：`po.py` 已按同一手法修好（2673 MB → 43 MB），但 `validate_po.py` 走的是 `aperture_po` 这条路径，所以「内存」这项只修了一半——最容易在 8 GB 笔记本上崩的那个脚本，恰恰没改到。`validate_po.py:43` 的 `THETA` 就是 4801 点，乘 `n_per_lambda=12` 必然 4 GB 起。

建议：照 `po.py` 的循环改法，把 165-166 两行换成
```python
    # 逐 theta 循环，不一次性构造 (n_theta, n_points) 的复数组（同 po.py 的处理）
    we = w * af.e
    field = np.empty(th.size, dtype=complex)
    for i in range(th.size):
        field[i] = np.sum(np.exp(1j * kx[i] * af.x) * we)
```
数值逐位不变，改完 `validate_po.py` 应从 10 GB 掉到几十 MB。

---

## C-02 | 高 | src/antennasim/pipeline.py:173-194

现象：同一屏上 `峰值增益损失 -0.479 dB` 和 `绝对增益 31.24 dBi` **加不起来**。

实测：
```
20λ : 35.96 - 4.49 - 0.479 = 30.99 dBi，程序印 31.24 dBi  → 差 +0.25
fd05: 35.96 - 6.11 - 0.440 = 29.42 dBi，程序印 29.65 dBi  → 差 +0.23
```
复现：跑两个算例，把「理想 − 溢出 − 峰值增益损失」与「绝对增益」比一下，两个算例都差 0.23~0.25 dB。

判据：`peak_db` 是 `20log10(|∫_未遮挡 A dA| / |∫_全口径 A dA|)`，是**场幅比**，隐含「被副面挡掉的那份功率仍然辐射出去了，只是不进主波束」；而 `eta_illum × (1-b)`（`:186-188`）隐含「被挡掉的功率不辐射」。两种归一约定正好差被挡区域那份功率，量级 0.2~0.3 dB。

我把未遮挡区单独做二重积分核过：`eta_ap = 0.9475`（→ -0.234 dB），与预算的 `0.9471` 只差 0.002 dB，说明**预算这一路本身是自洽的**，错的是把它和 `-0.479` 并列摆出来——两者不是同一个坐标系。

建议：二选一，并写进论文。推荐 A（符合方向性定义 `D = 4πU_max/P_rad`，`P_rad` 含全部辐射功率，遮挡散射出去的那份该算）。
**方案 A**：预算统一用 PO 的场幅比折成功率
```python
    # P_rad 取全口径功率（被遮挡部分仍辐射，不进主波束），故 (1-b) 不再单独乘
    eta_taper_block = 10.0 ** (peak_db / 10.0)      # -0.479 dB -> 0.8955
    eta_total = spill_eta * eta_taper_block
```
→ 20λ **30.99 dBi**、fd05 **29.42 dBi**（比现在各低 0.25 / 0.23）。
**方案 B**：保留 `0.3369`，但把 `:174` 改名 `口径场幅损失`，并注明「非增益损失，不可与总效率相加」。

不管选哪个，都要保证屏幕上这三个数能对上，否则答辩时一问就露。

---

## C-03 | 中 | src/antennasim/pipeline.py:199-214

现象：溢出、效率、绝对增益全都只 `say()` 到屏幕，`metrics` / `checks` 里一项都没有。

实测：
```bash
python run_pipeline.py cases/cassegrain_20lambda.yaml --no-figure
# 屏幕有：溢出发射效率 0.3557 / 总效率 0.3369 / 绝对增益 31.24 dBi
python -c "import json;print(list(json.load(open('data/results/cassegrain_20lambda.json'))['metrics']))"
# ['近焦点误差','遮挡面积比','射线命中率','准直误差','光程起伏','艾里自检偏差','峰值增益损失']
# checks: 5 项，全是几何/数值自检，无能量项
```

判据：`data/results/*.json` 是要进论文、进表格、被下游脚本消费的东西。现在 json 里唯一和增益有关的字段还是 `峰值增益损失 -0.479 dB`，谁拿 json 画图/填表都会得到「fd05 比 20λ 好 0.04 dB」的**反向结论**——屏幕对了、存档没对。另外 `fd05` 有 75% 功率溢出发射，`checks` 依然全绿。

建议：在 `metrics` 字典（`:199-207`）里补三项，`checks`（`:209-214`）补一条：
```python
        "溢出发射效率": f"{spill_eta:.4f}",
        "总效率(溢出+锥削+遮挡)": f"{eta_total:.4f}",
        "绝对增益": f"{gain_db:.2f} dBi",
```
```python
        "溢出可接受": spill_eta > 0.5,     # 约 -3 dB，只挡明显病态算例
```
阈值 0.5 是我的建议值，你可以按 fd05 该不该放行自己定；重点是有这条判据，别让 75% 溢出悄悄通过。

---

## C-04 | 中 | tests/（缺）

现象：溢出公式与整条效率预算零测试覆盖。

实测：
```bash
grep -rn "spill\|溢出\|eta_illum\|gain_db" tests/    # 空
```
`57 passed` 与本轮之前数量相同，这三个提交没加任何用例。

判据：这个仓库栽过同类错误——`q_for_edge_taper` 第一版漏因子 2，靠「正算再反解」才发现。溢出式把 `cos^(q+1)` 写成 `cos^q`，q=8 时结果是 **0.323** 而非 0.356（差 0.4 dB），**现有 57 条测试全绿**，只有肉眼看数字才可能发现。

建议：加进 `tests/test_feed.py`，用数值积分当裁判，不依赖被测公式自身：
```python
def test_spill_efficiency_matches_numeric_integral():
    """溢出效率 1-cos^(q+1)(θe) 必须与 cos^q 功率方向图的数值积分一致。"""
    from go_po.feed import cos_power_pattern
    te = np.radians(17.7613)          # 20λ/F0.4 算例的副面张角
    th = np.linspace(0.0, np.pi / 2, 400_001)
    for q in (4.0, 8.0, 47.0):
        g = cos_power_pattern(q)      # 返回的是电压 cos^(q/2)，平方才是功率
        pw = g(th) ** 2 * np.sin(th)
        num = np.trapezoid(pw[th <= te], th[th <= te]) / np.trapezoid(pw, th)
        ana = 1.0 - np.cos(te) ** (q + 1.0)
        assert abs(num - ana) < 1e-6, (q, num, ana)
```
（`np.trapezoid` 需 numpy≥2，兼容 1.x 用 `np.trapz`。我实测 q=8 解析/数值差 1.4e-6、q=47 差 1.1e-6。）

---

## C-05 | 中 | requirements.txt:16

现象：`requirements.txt` 缺 `pyyaml`，而官方安装脚本就是按它装的。

实测：
```bash
grep -i yaml requirements.txt                            # 空
grep -rn "import yaml" src/                              # src/antennasim/config.py:11
sed -n '/dependencies/,/^]/p' pyproject.toml             # 这里写了 "pyyaml"
grep -n "requirements.txt" scripts/setup-python-env.ps1   # :34 确实读 requirements.txt
```

判据：`setup-python-env.ps1` 装完，第一次 `run_pipeline.py` 就 `ModuleNotFoundError: No module named 'yaml'`。`pyproject.toml` 写对了不解决问题——文档和脚本走的都是 `requirements.txt`。这条上一轮提过，仍未修。

建议：在 `# 数据与文件` 段补一行，与 pyproject 对齐：
```
pyyaml>=6.0
```

---

## C-06 | 中 | README.md:25,28,29

现象：目录树里三格指向仓库中不存在的目录。

实测：
```bash
ls -1 src                     # antennasim common fenics go_po —— 没有 openems
git ls-files data | wc -l     # 0
ls -d notebooks               # 不存在
```

判据：README 是评审和答辩老师第一眼看的，指着三个空目录会被读成「规划了但没做」——这恰好是这份毕设现在最敏感的地方（`src/fenics/` 只有 81 行 `check_env.py`，openEMS 只有一个偶极子冒烟测试）。

建议：改成如实描述
```
  fenics/                       FEniCSx 结构有限元与面板变形
examples/cassegrain_20lambda/  统一算例
data/results/                  运行时生成的结果（不入库）
```
删掉 `openems/`、`notebooks/`；若要保留计划，标 `(计划中，尚未提交)`。

---

## C-07 | 低 | HANDOFF.md:3,129

现象：Ruze 公式写成 `exp(-sigma^2/2)`，与代码不一致；文件停在 2026-10-15。

实测：
```bash
sed -n '129p' HANDOFF.md                # G/G0 = exp(-sigma^2/2)
sed -n '161p' scripts/validate_po.py    # 10*log10(exp(-(phi_rms**2)))  → exp(-σ²)
sed -n '3p'  HANDOFF.md                 # 更新时间：2026-10-15 深夜
```

判据：写论文第 3 章时最容易照抄的就是交接文档；`σ²/2` 与 `σ²` 同样变形下差一倍（3 dB 量级），是会被评委抓住的错。同一文件还把「环境全通了」当结论，与 `src/fenics/` 只有 `check_env.py` 的现状对不上。

建议：129 行改成 `G/G0 = exp(-sigma_phi^2)`；文首加一行
`> 本文为 2026-10-15 历史记录，当前进度与未决问题见 work/进度台账.md`。

---

## C-08 | 低 | outputs/ 与 `*.bak-20261104`

现象：生成物、备份文件、Word 锁文件仍在库里，且持续污染工作区。

实测：
```bash
git status --short              # 跑完 validate_po → M outputs/validate_po.txt
git ls-files "*.bak*" | wc -l   # 5
git ls-files outputs/           # validate_po.png/.txt + ~$...docx（Word 锁文件）
```

判据：`.gitignore` 里写了 `*.bak`，但**对已跟踪文件无效**。`outputs/validate_po.txt` 每次运行都变，意味着任何一次正常跑脚本后工作区就是脏的——这正是「push 前先 pull」最容易撞车的地方，`git pull --rebase` 会因为未提交改动直接失败。`~$...docx` 是 Office 崩溃残留，没有内容价值。

建议：
```bash
git rm -r --cached --quiet -- "*.bak-20261104" "*.bak-20261104b" outputs/
printf '\n# 生成物\noutputs/*.txt\noutputs/*.png\n~$*\n' >> .gitignore
```
要给论文用的图，改名移到 `docs/figures/` 并在定稿时单独提交。

---

## C-09 | 低 | src/antennasim/pipeline.py:174

现象：`峰值增益损失` 印三位小数，但该数不随网格收敛，是假精度。

实测：
```
n_rho×n_phi     该值
200 ×320        -0.4786
400 ×640        -0.4898
800 ×1280       -0.4842
1600×2560       -0.4870
```
抖动 ±0.006 dB，无单调趋势。几何遮挡面积比 `4.913%` 单独只能解释 `10log10(1-b) = -0.22 dB`，剩下的来自「被挡的正好是照明最强的中心」，属正常物理效应（该口径的换算见 C-02）。

判据：遮挡台阶 `rho = Rb = 2.2165λ` 落在网格环之间，整个环的面积被整体判给或判走，结果随 `Rb/(R/n_rho)` 的余数振荡，不会随加密收敛。

建议：论文与报告统一写 **-0.49 ± 0.01 dB**（或「约 -0.5 dB」），格式串改 `{peak_db:.2f}`。真要三位小数，得对 `Rb` 所在环做亚网格面积加权，而不是继续加密网格。

---

## 优先级建议

1. **C-01** — 一行循环，10 GB → 40 MB，改完立刻能验证
2. **C-02** — 同屏三个数加不起来，答辩必问
3. **C-03 / C-04** — 存档口径 + 防回归
4. C-05 / C-06 / C-07 / C-08 / C-09 — 新人装得上、面子和卫生

改完在台账里按编号引用即可（`C-01 已修 @<sha>`）。要我复核就再 push，我拉最新版重跑同一套脚本对数字。
