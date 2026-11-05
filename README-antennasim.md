# antennasim 使用说明

把 openEMS、自研 GO/PO、FEniCSx 串成一条可复现链路的小框架。
本文档说明四个工程特性怎么用、以及**复试时怎么讲**。

---

## 一条命令跑完整个算例

```powershell
# 在仓库根目录下
.\.venv\Scripts\python.exe run_pipeline.py cases\cassegrain_20lambda.yaml
```

输出分三段：五个阶段的进度、校验指标汇总、生成的文件路径。

```
[1/5] 建模几何
      主面 D=20λ  F=8λ  M=4  L=6.4000λ
      副面口径 4.4329λ   遮挡比 4.913%
[2/5] 射线追踪
      命中 1872/1872   准直误差 4.44e-16   光程起伏 7.11e-15
[3/5] 构造口径场
      口径边缘锥削 -1.91 dB  (q=8)
      溢出发射效率 0.3557  (-4.49 dB)   副面张角 17.76°
[4/5] 远场积分
      均匀口径自检  max|数值-解析| = 7.01e-06
      口径场幅损失 -0.48 dB（遮挡台阶在网格之间，只到 ±0.01 dB）
      口径效率 0.9482（锥削 0.9960 × 遮挡 0.9509，诊断用）
        → 总效率 0.3372  绝对增益 31.24 dBi  (理想 35.96 dBi)
[5/5] 保存结果
```

常用参数：

```powershell
# 只输出最终汇总
run_pipeline.py cases\cassegrain_20lambda.yaml --quiet

# 跳过出图（跑批量扫描时更快）
run_pipeline.py cases\cassegrain_20lambda.yaml --no-figure

# 输出目录相对指定根目录解析
run_pipeline.py cases\cassegrain_20lambda.yaml --project-root D:\run1
```

**复试时可以这样讲**：
> "整条链路是配置驱动的，一条命令从建模跑到出图。跑一次大约 4 秒，
> 每个环节都会打印校验指标，任何一项不达标会在我看到结果之前就暴露出来。"

---

## 配置文件：改参数不用改代码

算例参数全部放在 `cases/*.yaml`：

```yaml
case:
  name: cassegrain_20lambda

geometry:
  diameter: 20.0        # 主面口径，单位波长
  f_over_d: 0.4         # 焦径比
  magnification: 4.0    # 双曲面放大率
  sub_vertex: null      # null 表示用默认 M*F/(M+1)

feed:
  q: 8.0                # 馈源方向图 cos^q(theta)

sampling:
  n_radial: 26          # 射线采样密度
  n_rho: 200            # 口径积分采样密度
```

**演示参数扫描**——仓库里已经有两个算例：

| 算例 | F/D | 遮挡比 | 口径场幅损失（场幅比） | 绝对增益 |
|---|---|---|---|---|
| `cassegrain_20lambda.yaml` | 0.4 | 4.913 % | −0.48 dB | **31.24 dBi** |
| `cassegrain_fd05.yaml` | 0.5 | 4.551 % | −0.44 dB | 29.65 dBi |

注意「口径场幅损失」是 `20log10` 的场幅比，**不是增益损失**；功率等效约 −0.22 dB。真正进预算的是绝对增益那一列。

两个算例的差别**只有配置文件里的一行**，代码一个字没动。

```powershell
run_pipeline.py cases\cassegrain_20lambda.yaml --quiet
run_pipeline.py cases\cassegrain_fd05.yaml --quiet
```

配置对象的读写也能直接调：

```python
from antennasim.config import load_config

cfg = load_config("cases/cassegrain_20lambda.yaml")
print(cfg.focal_length, cfg.blocking_ratio)   # 派生量自动算
cfg.dump("cases/my_case.yaml")                # 存回文件
```

非法的配置会在读取阶段就报错，不会带着错参数跑完再出问题：

```python
load_config("bad.yaml")
# ConfigError: geometry.magnification 必须 > 1，当前 0.5
```

**复试时可以这样讲**：
> "参数和代码分离。做一个参数扫描只需要生成不同的配置文件，
> 不需要动代码，也就不会引入人为错误。配置在读取阶段就会校验，
> 非法取值直接报错。"

---

## 测试：把物理判据变成自动回归

```powershell
.\.venv\Scripts\python.exe -m pytest tests -q
# 33 passed in 27.76s
```

测试不是走过场，每条都对应一个**物理判据**：

| 测试文件 | 验证什么 | 判据来源 |
|---|---|---|
| `test_geometry.py` | 双曲面近焦点是否精确落在抛物面焦点 | 卡塞格伦共焦条件 |
| `test_rays.py` | 出射方向是否平行于视轴 | 几何光学准直 |
| `test_rays.py` | 光程 r1+r2−z 是否为常数 | 等光程原理 |
| `test_po.py` | 数值积分是否等于解析艾里图案 | 均匀圆口径解析解 |
| `test_config.py` | 非法配置是否被拒绝 | 输入校验 |
| `test_cli.py` | 命令行与配置的优先级是否正确 | 接口契约 |
| `test_pipeline.py` | 同配置两次运行结果是否一致 | 可复现性 |

挑几条看：

```
tests/test_geometry.py::test_near_focus_coincides_with_parabola_focus PASSED
tests/test_rays.py::test_collimation_error_is_at_machine_precision PASSED
tests/test_rays.py::test_optical_path_is_constant PASSED
tests/test_po.py::test_uniform_aperture_matches_airy PASSED
tests/test_pipeline.py::test_pipeline_is_reproducible PASSED
```

**复试时可以这样讲**：
> "我把验证判据写成了自动化测试。核心三条是：几何自检、射线追踪的两个
> 物理不变量、以及口径积分与解析解比对。目前 33 条测试全过。
> 改代码之后跑一遍 pytest，任何一条不达标立刻能发现。"

这一条是最能体现"科研习惯"的，多数本科生没有。

---

## 包结构：代码怎么组织

```
antennasim/
├── run_pipeline.py          一条命令入口
├── pyproject.toml           工程元数据与工具配置
├── cases/                   算例配置文件
│   ├── cassegrain_20lambda.yaml
│   └── cassegrain_fd05.yaml
├── tests/                   自动化测试（33 条）
├── src/
│   ├── antennasim/          本框架：配置、流水线、命令行
│   │   ├── config.py
│   │   ├── pipeline.py
│   │   └── cli.py
│   ├── common/              公共几何与常量
│   │   └── cassegrain.py
│   └── go_po/               自研数值内核
│       ├── rays.py          射线追踪
│       ├── po.py            物理光学积分
│       └── feed.py          馈源方向图
└── data/results/            结果输出（npz + json + png）
```

分层原则：**框架负责编排，内核负责计算**。
`antennasim` 不碰数值细节，`go_po` 不关心怎么配置和输出。

**复试时可以这样讲**：
> "代码分成三层：几何与数值内核、流水线编排、命令行接口。
> 每层职责单一，所以能单独测试，也方便换掉其中一层。"

---

## 结果文件

每次运行产出三个文件（在 `data/results/`）：

| 文件 | 内容 |
|---|---|
| `<tag>.npz` | 数值结果：方向图、口径幅度、关键几何量 |
| `<tag>.json` | 校验指标 + 完整配置快照 |
| `<tag>.png` | 两张图：口径幅度分布、远场方向图 |

JSON 里同时存了配置，所以**每条结果都能追溯到产生它的参数**——
这是可复现性的基础。

---

## 已知待办

- `data/results/cassegrain_fd05.png` 是修复配置优先级 bug **之前**留下的旧文件，
  当前配置写的是 `figure: false`，重跑不会再生成。可以手动删除。
- 口径幅度还是一阶模型 `g(θ)/(r1+r2)`，严格做法需要跟踪射线管 Jacobian。
- 结构域（FEniCSx）与双域耦合尚未接入流水线，目前只覆盖射频侧。
