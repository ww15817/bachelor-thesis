# 代码审查 2026-11-04（WorkBuddy）

审查对象：**`12af5df`**（本轮含 `f4cd0fb` / `94d5af7` / `2d61133` / `12af5df`）
修订记录：首版基于 `2d61133`，已由 `17e84b378` 入库；本版对齐 `12af5df`，C-03 缩小范围、新增 C-10。
方法：独立 venv（numpy 2.x + scipy）实跑 `pytest`、两个算例、`work/scan_q.py`、`scripts/validate_po.py`；关键数字用独立脚本复算，不采信代码里的自述。
约定：**只写本文件，不改 `src/`**；push 前先 `git pull --rebase`。

---

## 本轮已确认修好（5 项，不必再动）

| 原问题 | 现状 | 我的复测 |
|---|---|---|
| `far_field_aperture` 一次性构造 `(n_theta,n_rho,n_phi)` | 已改逐 θ 循环 | 峰值 **2673 MB → 43 MB**，耗时 3.75→3.76 s；与旧实现 `max\|diff\| = 0.0`（逐位相同） |
| `work/dsh-session/` 2.4 MB 会话日志入库 | 已 `git rm` | 目录已不在仓库 |
| 溢出效率完全不在账上 | 已打印诊断 | `0.3557 (-4.49 dB)`，与我独立数值积分 `0.355666` 一致；`1-cos^(q+1)(θe)` 无误 |
| F/D=0.5 对照算例结论方向反了 | 已补绝对增益 | fd05 `29.65 dBi` vs 20λ `31.24 dBi`，正确体现 fd05 更差 1.59 dB |
| `metrics` / 汇总表 缺溢出与绝对增益 | `12af5df` 已补 5 个字段 | 汇总表现在会打出「溢出发射效率/照明效率/总效率/绝对增益/理想增益」全 5 项，json 同步 |

顺带核过、都没问题：`feed_half_angle_for_edge` 的 17.7613° 与「馈源→副面边缘」几何角差 0.000 角分（用这个角做溢出积分是对的）；`57 passed`；L1/L2/L3 全过（L1 偏差 1.4e-5，L3 与 Ruze 差 0.001 dB）；q 扫描的最优解是内部极大（q=64 回落到 34.77），不是扫描边界假象。

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

判据：`po.py` 已按同一手法修好（2673 MB → 43 MB），但 `validate_po.py` 走的是 `aperture_po` 这条路径，所以「内存」这项只修了一半——最容易在 8 GB 笔记本上崩的那个脚本，恰恰没改到。`validate_po.py:43` 的 `THETA` 就是 4801 点，乘 `n_per_lambda=12` 必然 4 GB 起。该文件到 `12af5df` 仍未改动（我逐字节比过）。

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

判据：`peak_db` 是 `20log10(|∫_未遮挡 A dA| / |∫_全口径 A dA|)`，是**场幅比**，隐含「被副面挡掉的那份功率仍然辐射出去了，只是不进主波束」；而 `eta_illum × (1-b)`（`:186-188`）隐含「被挡掉的功率不辐射」。两种归一约定正好差被挡区域那份功率（量级 0.2~0.3 dB）。

我把未遮挡区单独做二重积分核过：`eta_ap = 0.9475`（→ -0.234 dB），与预算的 `0.9471` 只差 0.002 dB，说明**预算这一路本身是自洽的**，错的是把它和 `-0.479` 并列摆出来——两者不是同一个坐标系。

**这不是记账洁癖，它会改变设计结论。** 按两种约定分别扫 q（`work/scan_q.py` 的 q 表）：

| q | 溢出 | 峰值增益损失 | 方案 B（现） | 方案 A |
|---|---|---|---|---|
| 8 | -4.49 | -0.479 | 31.24 | 31.00 |
| 47 | -0.44 | -0.723 | **34.85** | 34.80 |
| 64 | -0.19 | -0.846 | 34.77 | **34.93** |

方案 B 最优 q=47，方案 A 最优 q=64。**推荐值 `34.85 dBi / 净增 3.61 dB` 依赖于约定**——定下约定之前，这一节别急着写进论文。另外方案 A 的最优落在扫描上界，说明还得往 q>64 再扫一段才能确认。

建议：二选一，并写进论文。推荐 A（符合方向性定义 `D = 4πU_max/P_rad`，`P_rad` 含全部辐射功率，遮挡散射出去的那份该算）。
**方案 A**：预算统一用 PO 的场幅比折成功率
```python
    # P_rad 取全口径功率（被遮挡部分仍辐射，不进主波束），故 (1-b) 不再单独乘
    eta_taper_block = 10.0 ** (peak_db / 10.0)      # -0.479 dB -> 0.8955
    eta_total = spill_eta * eta_taper_block
```
→ 20λ **30.99 dBi**、fd05 **29.42 dBi**（比现在各低 0.25 / 0.23）。
**方案 B**：保留 `0.3369`，但把 `:174` 改名 `口径场幅损失`，并注明「非增益损失，不可与总效率相加」。

---

## C-03 | 低 | src/antennasim/pipeline.py:214-220

现象：`checks` 五项仍全是几何/数值自检，没有一条管能量分配；75% 功率溢出的算例照样判「通过」。
（`metrics`/汇总表已由 `12af5df` 补齐，这条只剩 `checks`。）

实测：
```bash
python run_pipeline.py cases/cassegrain_fd05.yaml --no-figure | tail -14
# 溢出发射效率 0.2452（即 75% 功率没打到副面）
# ... 状态 通过
```
```python
# pipeline.py:214-220 五项依次为：几何自检 / 射线全部命中 / 准直误差达标 / 等光程达标 / 艾里自检达标
```

判据：这五项验的是「几何和积分器对不对」，全部通过只说明代码没写错，不说明这个馈源配这个副面是合理设计。`fd05` 的 0.2452 在物理上已经很难看，却拿不到任何提示。

建议：`checks` 补一条（0.5 约 -3 dB，够宽，只挡明显病态算例）：
```python
        "溢出可接受": spill_eta > 0.5,
```
阈值你可以按 fd05 该不该放行自己定；重点是有这条判据，别让 75% 溢出悄悄通过。

---

## C-04 | 中 | tests/（缺）

现象：溢出公式与整条效率预算零测试覆盖。

实测：
```bash
grep -rn "spill\|溢出\|eta_illum\|gain_db" tests/    # 空
```
四个提交下来，测试数仍是 `57 passed`，一条没加。

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

判据：`setup-python-env.ps1` 装完，第一次 `run_pipeline.py` 就 `ModuleNotFoundError: No module named 'yaml'`。`pyproject.toml` 写对了不解决问题——文档和脚本走的都是 `requirements.txt`。这条已连续两轮未修。

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
find . -name '~$*'              # outputs/~$...毕业设计规划书.docx（Word 锁文件）
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

## C-10 | 中 | work/scan_q.py:46

现象：q 扫描把基线 `31.24 dBi` 写死在比较式里；同时 `FeedConfig(q=q)` 会丢掉 yaml 里其余的 feed 字段。

实测：
```bash
python work/scan_q.py
# q=47 → 34.85 dBi，末行「相对 q=8（31.24 dBi）净增 +3.61 dB」
sed -n '46p' work/scan_q.py   # print(f"...净增 {best[4] - 31.24:+.2f} dB")
sed -n '33p' work/scan_q.py   # cfg = dataclasses.replace(base, feed=FeedConfig(q=q))
```
跑 `python -W error work/scan_q.py` 会直接崩在第 8 行：`SyntaxWarning: invalid escape sequence '\.'`（docstring 里的 `.\venv\Scripts`）。

判据：这行的基线正是 C-02 要动的那个数。一旦按方案 A 改，q=8 变成 31.00、q=47 变成 34.80，脚本仍印 `34.80 - 31.24 = +3.56`，而真实净增是 **+3.80**——「优化净收益」这个要进论文的数会静默错 0.24 dB，且没人会发现，因为脚本自己看起来是自洽的。

`FeedConfig(q=q)` 的问题同理：`dataclasses.replace` 只替换 feed 整体，`edge_taper_db` 等字段被重置成默认值。当前 `cassegrain_20lambda.yaml` 只写了 `q`，暂时无碍，但哪天配置里加了 `edge_taper_db`，这条路径会悄悄绕过 `:142-148` 那个偏差告警。

建议：
```python
    # 基线从本次扫描里取，不写字面量（C-02 改了口径后这里自动跟着变）
    base_q8 = next(g for (qq, _s, _i, _t, g) in rows if qq == 8.0)
    print(f"相对 q=8（{base_q8:.2f} dBi）净增 {best[4] - base_q8:+.2f} dB")
```
```python
    cfg = dataclasses.replace(
        base, feed=dataclasses.replace(base.feed, q=q)   # 只改 q，保留其余 feed 字段
    )
```
docstring 改 `r"""..."""`（或把反斜杠写双）。另外把 `Q_LIST` 上界开到 96：方案 A 的最优落在 q=64 这个扫描边界上，不确认一下分不清是极大还是还没扫到头。

---

## 优先级建议

1. **C-02** — 先定归一约定，否则 q 扫描和绝对增益都是浮的（最优 q 会在 47/64 之间摇摆）
2. **C-01** — 一行循环，10 GB → 40 MB，改完立刻能验证
3. **C-10 / C-03 / C-04** — 让 q 扫描的数字站得住
4. C-05 / C-06 / C-07 / C-08 / C-09 — 新人装得上、面子和卫生

改完在台账里按编号引用即可（`C-01 已修 @<sha>`）。要我复核就再 push，我拉最新版重跑同一套脚本对数字。

---

# 复审 2026-11-05（基线 `2fb0ef0`）

审查对象：`2fb0ef08c`（修复提交 `111182fa9`，台账 W-0015）。
方法同前：独立 venv 实跑 `pytest`、两个算例、`work/scan_q.py`、`scripts/validate_po.py`；关键结论用**解析可判的判例**（均匀口径+遮挡）交叉验证，网格 200×320 / 400×640 / 800×1280 三档查收敛。

> **先说一句要紧的：C-02 的方案 A 是我给错的。** 你照着实现了，实现也忠实，但式子本身把遮挡算了两遍，见 C-11。旧值 `31.24 dBi` 本来就是对的。

## 本轮已确认修好（8 项）

| 编号 | 现状 | 我的复测 |
|---|---|---|
| C-01 | `aperture_po.py` 已改逐 θ 循环 | `validate_po.py` 全程峰值 **9993 MB → 129 MB**，138 s → 113.8 s，L1/L2/L3 全过；与建议代码逐行一致 |
| C-04 | 补 `test_spill_efficiency_matches_numeric_integral` | 58 passed；把 `cos^(q+1)` 写成 `cos^q` 会精确挂掉（q=8 得 0.323 vs 0.356），是真判据不是恒真 |
| C-05 | `requirements.txt:9` 补 `pyyaml>=6.0` | ✓ |
| C-06 | README 目录树改成实际结构，并注明 `src/openems/`、`notebooks/` 尚未提交 | ✓ |
| C-07 | HANDOFF 加了 `exp(-σφ²)` 的更正注 | ✓（保留原文+更正注，可接受） |
| C-08 | `outputs/` 全部出库（含 Word 锁文件），两份 docx 归到 `docs/` | 远程树 `outputs` tracked = 0 ✓；`.bak` 漏了，见 C-15 |
| C-09 | 只印两位小数并注明 ±0.01 dB | ✓ |
| C-10 | 基线不再写死、`dataclasses.replace` 保留 feed 字段、`Q_LIST` 开到 96、docstring 改 `r"""` | ✓ 四小项全中 |

---

## C-11 | 高 | src/antennasim/pipeline.py:199-200

现象：`eta_taper_block = 10.0 ** (peak_db / 10.0)` 把遮挡算了两遍——`peak_db` 已经是「未遮挡/全口径」的场幅比，它再平方一次又乘回去一份遮挡。绝对增益因此从 31.24 掉到 31.00（fd05 从 29.65 掉到 29.42）。

实测：拿**解析可判**的均匀口径+遮挡做判例（真值恒为 `η = 1-b`，b=4.913%）：
```
真值  1-b                      = 0.950873  -> -0.219 dB
精确  |∫un|²/(A·∫un²)          = 0.950938  -> -0.218 dB   ✓
现实现 (|∫un|/|∫all|)²         = 0.904283  -> -0.437 dB   ✗ 正好是两倍 dB
```
真实 q=8 锥削口径、三档网格（数值稳定到 ±0.001）：
```
精确  0.94819/0.94709/0.94764 -> 31.24 dBi（三档都是 31.24）
现实现 0.89565/0.89334/0.89450 -> 31.00 / 30.98 / 30.99
差 -0.25 dB（与 1-b 的量级一致，说明就是多乘了一份遮挡）
```

判据：方向性 `D = 4π|∫A dA|² / (λ² ∫|A|² dA)`，分母是 `∫|A|²`，不是 `(∫A)²/A`。`10**(peak_db/10)` 只在 A 均匀时才勉强对上，而均匀时它给的又是 `(1-b)²` —— 教科书遮挡口径是 `1-b`，差一倍。换句话说：**场幅比的平方不是效率，效率必须回到二重积分里算。**

建议：别从场幅比反推，直接积分未遮挡区：
```python
    amp_on_grid = np.interp(grid_rho, rho_tab, amp_tab)
    amp_un = np.where(grid_rho >= block_radius, amp_on_grid, 0.0)
    # 遮挡口径的方向性：分子只取未遮挡区，分母 ∫|A|² 也只取未遮挡区
    eta_taper_block = float(np.sum(amp_un * grid_dA)) ** 2 / (
        aperture_area * float(np.sum(amp_un ** 2 * grid_dA))
    )
    eta_total = spill_eta * eta_taper_block
```
→ 20λ **31.24 dBi**、fd05 **29.65 dBi**（= 旧值）。
想看分解就印成「锥削 × 遮挡」：`0.9960 (-0.017 dB) × 0.9513 (-0.217 dB，解析 10log10(1-b) = -0.219) = 0.9475 (-0.234 dB)`，三者自洽，答辩时好讲。
并且 `-0.479 dB` 必须标注为「20log10 场幅比」，它不是增益损失——真写成损失就是 C-12 那个错。

---

## C-12 | 高 | docs/01-技术路线.md:113

现象：验证表里写「副面遮挡增益损失 | `~20log10(1-eps)` | -0.479 dB（eps=4.913%，理论 -0.4376 dB）」——用 20log10 算出来的量，却挂在「增益」名下。

实测：
```
10*log10(1-0.04913) = -0.219 dB   （功率/增益，教科书遮挡口径 G/G0 = 1-b）
20*log10(1-0.04913) = -0.437 dB   （场幅比，不是增益）
数值 -0.479 dB                     （场幅比 + 锥削贡献 0.04 dB）
```

判据：这是整条 C-02 混乱的源头，而且它还活着。论文里若照抄「遮挡增益损失 -0.479 dB」，读者按 `D = 4πU_max/P_rad` 复算会差 0.26 dB，且差得没道理（多了个平方）。上一轮我把这个混用判成「归一约定之争」，是我当时没回到定义——抱歉。

建议：拆成两行，把「场幅」和「功率」分开写：
```
| 副面遮挡（场幅比） | 20log10(1-ε) | -0.437 dB（数值 -0.479，差 0.04 来自锥削） |
| 副面遮挡（功率）   | 10log10(1-ε) | -0.219 dB（数值 -0.234） |
```

---

## C-13 | 中 | work/scan_q.py:6 与 work/进度台账.md（W-0015）

现象：台账 W-0015 记「方案 A 下 q 最优 64 / 34.93 dBi」，这个结论建立在 C-11 那个多算一份遮挡的公式上。

实测：同一套口径表，两种公式扫 q（800×1280，精确列三档网格一致）：
```
 q   溢出dB   场幅dB   精确G    现方案A G   差
 8   -4.49   -0.486   31.24    30.99    -0.25
16   -2.49   -0.532   33.20    32.94    -0.26
32   -0.97   -0.631   34.57    34.37    -0.21
47   -0.44   -0.734   34.89    34.79    -0.10
64   -0.19   -0.859   34.83    34.92    +0.09
96   -0.04   -1.121   34.28    34.80    +0.52
```
精确公式最优 **q≈47 → 34.89 dBi**（q=8 为 31.24，净增 +3.65 dB）；q=47 与 q=64 只差 **0.06 dB**，是个平台，不是一个点。

判据：误差随 q **变号**（低 q 低估 0.25、高 q 高估 0.52），所以它不只是整体偏 0.25 dB，而是把最优从 47 推到了 64。上轮我说「C-02 会改变设计结论」，方向说反了——真正的原因是公式错，不是约定选错。

建议：先按 C-11 改公式再重跑；台账 W-0015 补一条更正：「精确公式下 q=47 与 64 处 34.89 / 34.83 dBi，属 0.06 dB 平台，基线 q=8 = 31.24 dBi，净增 +3.65 dB」。论文那一节按**平台**写（「q 取 40~70 均可，建议 50 左右」），别给单点最优——0.06 dB 在模型不确定度内，答辩被问「凭什么正好 47」答不上来。

---

## C-14 | 中 | README.md:91 / README-antennasim.md:26,73,75 / docs/04-环境配置任务清单-v2.md:247 / docs/04-环境配置任务清单.md:232 / examples/cassegrain_20lambda/m3_cassegrain.py:139-140,217 / work/进度台账.md:368

现象：指标改名（`峰值增益损失` → `口径场幅损失`）只改了 `pipeline.py` 和 `tests/`，6 处文档与示例仍在用旧名和旧值；台账 368 行还写着「`总效率 = 溢出 × 照明 × 遮挡`」，与新代码（`总效率 = 溢出 × 场幅比²`）直接矛盾。

实测：
```bash
grep -rn "峰值增益损失" README.md README-antennasim.md docs/ examples/
# README.md:91  -0.479 dB
# README-antennasim.md:26 / :73 / :75
# docs/04-环境配置任务清单-v2.md:247  -0.479 dB
# docs/04-环境配置任务清单.md:232      -0.525 dB（更旧）
# m3_cassegrain.py:140 / :217
```

判据：`m3_cassegrain.py:139` 的 `g_blocked - g_clear` 是两个 `20log10` 场幅之差，正是 C-12 那个错，而它还被画成了子图的 y 轴标签（`:217`），会直接进答辩 PPT。名称错、数值对不上新口径、台账公式错，三样凑在一起，论文写的时候必然抄到错的那份。

建议：统一改成「口径场幅损失（20log10 场幅比，非增益损失）」；m3 里补一句功率等效：
```python
    print(f"[PO] 遮挡造成的口径场幅损失 = {g_blocked - g_clear:.3f} dB"
          f"（场幅比；功率等效约 {(g_blocked - g_clear) / 2:.3f} dB）")
```
轴标签同步改成「口径场幅损失 [dB]」。台账 368 改成「`总效率 = 溢出 × 口径(锥削×遮挡)`」。

---

## C-15 | 中 | 5 个 `*.bak-20261104`

现象：C-08 只把 `outputs/` 出库了，5 个备份文件仍在远程树里。

实测（GitHub git tree，119 个 blob）：
```
src/go_po/po.py.bak-20261104
src/go_po/po.py.bak-20261104b
src/go_po/aperture_po.py.bak-20261104
src/antennasim/pipeline.py.bak-20261104
examples/cassegrain_20lambda/m3_cassegrain.py.bak-20261104
.gitignore:38  *.bak          <-- 匹配不到 .bak-20261104
```

判据：两处叠加——(1) `.gitignore` 写的 `*.bak` 只匹配以 `.bak` **结尾**的文件，这里结尾是日期，得 `*.bak-*`；(2) 即便模式写对，gitignore 也拦不住已入库文件，必须 `git rm --cached`。上一轮 `work/dsh-session/` 就是栽在 (2)，这次 (1)(2) 一起来了。

建议：
```bash
git rm --cached src/go_po/po.py.bak-20261104 src/go_po/po.py.bak-20261104b \
              src/go_po/aperture_po.py.bak-20261104 src/antennasim/pipeline.py.bak-20261104 \
              examples/cassegrain_20lambda/m3_cassegrain.py.bak-20261104
# .gitignore 里把 *.bak 改成
*.bak*
```

---

## C-16 | 低 | tests/test_feed.py:201,231 与 requirements.txt:2

现象：测试用了 `np.trapezoid`（numpy **2.0** 才有），`requirements.txt` 写的是 `numpy>=1.26`。

实测：`np.trapezoid` 在 1.26 上不存在，只有 `np.trapz`（2.0 起弃用）。装到 1.26.x 的环境里 `pytest` 直接 `AttributeError`。

判据：README 承诺「照 requirements 装就能跑」，但版本下界允许一个跑不了的组合。`np.trapezoid` 出现在 `:201`（旧测试）和 `:231`（本轮新增），说明早就默认 2.x 了，只是没写进清单。

建议：`requirements.txt:2` 改 `numpy>=2.0`（更干脆，代码里已按 2.x 写），或测试里 `getattr(np, "trapezoid", np.trapz)` 兜底。

---

## C-17 | 低 | tests/test_pipeline.py:1

现象：这个文件带 UTF-8 BOM（本轮只改 3 行，被编辑器写回去了）。

实测：`head -c 3 tests/test_pipeline.py | od -An -tx1` → `ef bb bf`。
`src/` 与 `tests/` 里只剩它一个带 BOM（`work/` 下另有 10 个历史脚本带）。

判据：Python 能吃 BOM，但 `git diff` 第一行永远显示为改动（本轮 diff 就是这样），合并时容易冲突；仓库之前专门清过 BOM，这是回潮。

建议：编辑器切「UTF-8 无 BOM」保存。清已入库的那个：
```bash
python -c "p='tests/test_pipeline.py'; s=open(p,encoding='utf-8-sig').read(); open(p,'w',encoding='utf-8',newline='\n').write(s)"
```

---

## 优先级（本轮）

1. **C-11** — 5 行改动，31.00 → 31.24，且它决定 C-13 的全部数字
2. **C-12** — 论文那行表格，改个 log 系数，别把场幅比当增益写进「验证」章
3. **C-13** — 改完 C-11 重跑；结论按 34.89 dBi / 平台写，台账 W-0015 加更正
4. **C-14 / C-15** — 改名传播 + 备份出库，属于「老师第一眼」那一类
5. C-16 / C-17 — 版本下界与 BOM

改完仍按编号在台账里引用即可。
