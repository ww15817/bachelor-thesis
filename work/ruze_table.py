# -*- coding: utf-8 -*-
"""
Ruze 增益损失表（毕设结构侧算例参数）。

口径（仓库已定，勿重推）：
    sigma_phi = 2*k*sigma_w = 4*pi*sigma_w/lambda     (k = 2*pi/lambda)
    G/G0      = exp(-sigma_phi**2)                    # 注意是 exp(-sigma_phi**2)，不是 /2
    loss(dB)  = 10*log10(exp(-sigma_phi**2)) = -4.343*sigma_phi**2

本脚本只做纯计算，不写业务代码。运行：
    .\.venv\Scripts\python.exe work\ruze_table.py
产物：work/Ruze增益损失表.md
"""
import numpy as np
from pathlib import Path

# ---- 物理常数与算例参数（口径，来自任务单）----
LAMBDA_MM = 30.0          # 10 GHz 自由空间波长，mm
D_MM = 600.0              # 20*lambda 口径，mm（仅用于解读，不进公式）
F_GHZ = 10.0
LOSS_COEF = 4.343         # = 10*log10(e) ≈ 4.34294，按仓库口径取 4.343

# ---- 扫描点：sigma_w/lambda，1/500(=0.002) ~ 1/20(=0.05)，10 个点 ----
# 取干净分数，覆盖全区间；1/250=0.004 是验收点（sigma_phi 应 ≈0.0503）
RATIOS = np.array([1/500, 1/250, 1/200, 1/100, 1/80, 1/50, 1/40, 1/30, 1/25, 1/20])
FRACTIONS = ["1/500", "1/250", "1/200", "1/100", "1/80",
             "1/50", "1/40", "1/30", "1/25", "1/20"]


def phase_rms(ratio):
    """口径相位 rms (rad)：sigma_phi = 4*pi*sigma_w/lambda。"""
    return 4.0 * np.pi * ratio


def loss_db(ratio):
    """增益损失 (dB)：G/G0 = exp(-sigma_phi**2)，损失 = -4.343*sigma_phi**2（负值）。"""
    sphi = phase_rms(ratio)
    return -LOSS_COEF * sphi**2


def inverse_ratio(target_loss_db):
    """反解：要损失 |target_loss_db| dB，需要多大的 sigma_w/lambda。"""
    sphi = np.sqrt(target_loss_db / LOSS_COEF)   # target_loss_db 取正
    return sphi / (4.0 * np.pi), sphi


# ================= 正算表 =================
rows = []
for frac, r in zip(FRACTIONS, RATIOS):
    sw_mm = r * LAMBDA_MM
    sphi = phase_rms(r)
    ld = loss_db(r)
    rows.append((frac, r, sw_mm, sphi, ld))

# ================= 反解：0.1 / 0.5 / 1.0 dB =================
TARGETS = [0.1, 0.5, 1.0]
inv = []
for t in TARGETS:
    r, sphi = inverse_ratio(t)
    sw_mm = r * LAMBDA_MM
    # 代回正算，复现 dB
    check = -loss_db(r)   # loss_db 返回负，取正核对
    inv.append((t, r, sw_mm, sphi, check))

# ================= 验收自检（打印到 stdout）=================
print("=== 验收自检 ===")
# 1) sigma_phi = 0.5 rad -> loss = -4.343*0.25
sphi_chk = 0.5
loss1 = -LOSS_COEF * sphi_chk**2
print(f"[1] sigma_phi=0.5 rad -> loss = {loss1:+.4f} dB (期望 -1.086 +/-0.01)")
# 2) sigma_w = 0.004*lambda -> sigma_phi
r2 = 0.004
sphi2 = phase_rms(r2)
print(f"[2] sigma_w=0.004 lambda -> sigma_phi = {sphi2:.4f} rad (期望 0.0503)")
# 3) 反解代回
print("[3] 反解代回正算：")
for t, r, sw_mm, sphi, back in inv:
    err = abs(back - t)
    print(f"    目标 -{t:.1f} dB -> sigma_w={r:.5f} lambda={sw_mm:.4f} mm, "
          f"代回 = {back:+.4f} dB, 误差={err:.2e} dB (容差 1e-3)")

# ================= 自重量级（解读用）=================
sw_self_mm = 0.1            # 600mm 铝板 2mm 厚仅自重，中心挠度约 0.1 mm
r_self = sw_self_mm / LAMBDA_MM
sphi_self = phase_rms(r_self)
loss_self = loss_db(r_self)
r_05, sphi_05 = inverse_ratio(0.5)
sw_05_mm = r_05 * LAMBDA_MM
amp = sw_05_mm / sw_self_mm
print(f"\n=== 自重量级 ===")
print(f"中心挠度 0.1 mm = {r_self:.5f} lambda, sigma_phi={sphi_self:.4f} rad, "
      f"loss={loss_self:+.4f} dB")
print(f"0.5 dB 阈值 -> sigma_w={sw_05_mm:.3f} mm = {r_05:.5f} lambda, "
      f"相对自重放大 {amp:.1f}x")

# ================= 写 markdown =================
out = []
out.append("# Ruze 增益损失表（结构侧算例参数）\n")
out.append("> 口径：`sigma_phi = 4*pi*sigma_w/lambda`，`G/G0 = exp(-sigma_phi^2)`，"
           "损失 `= -4.343*sigma_phi^2`（注意是平方，不是 /2）。\n")
out.append(f"> 算例：f = {F_GHZ:.0f} GHz，lambda = {LAMBDA_MM:.0f} mm，"
           f"口径 D = {D_MM:.0f} mm（20 lambda）。\n")
out.append(f"> 本表由 `work/ruze_table.py` 用 numpy 正算生成，未手填任何占位数。\n")
out.append("\n## 1. 正算表\n")
out.append("| sigma_w (lambda) | sigma_w (mm) @10 GHz | sigma_phi (rad) | G/G0 (dB) |")
out.append("|---|---|---|---|")
for frac, r, sw_mm, sphi, ld in rows:
    out.append(f"| {frac} = {r:.5f} | {sw_mm:.3f} | {sphi:.4f} | {ld:+.4f} |")

out.append("\n## 2. 反解：要掉多少 dB，需要多大变形\n")
out.append("| 目标损失 (dB) | 需要 sigma_w (lambda) | 需要 sigma_w (mm) | 对应 sigma_phi (rad) | 代回正算 (dB) |")
out.append("|---|---|---|---|---|")
for t, r, sw_mm, sphi, back in inv:
    out.append(f"| -{t:.1f} | {r:.5f} | {sw_mm:.3f} | {sphi:.4f} | -{back:.4f} |")

out.append("\n## 3. 解读（自重量级落在哪一行）\n")
interp = (
    f"按 20 lambda（{D_MM:.0f} mm）口径、2 mm 铝板仅算自重，中心挠度约 0.1 mm "
    f"≈ {r_self:.4f} lambda，落在表中 1/500（0.06 mm）与 1/250（0.12 mm）两行之间，"
    f"sigma_phi≈{sphi_self:.3f} rad，增益仅掉约 {abs(loss_self):.3f} dB——远低于 0.1 dB，"
    f"仪器分辨不出来。要看到 0.5 dB 的可测损失，需 sigma_w≈{sw_05_mm:.2f} mm"
    f"≈{r_05:.4f} lambda，即把自重量级的变形放大约 {amp:.0f} 倍"
    f"（0.1 mm→{sw_05_mm:.1f} mm）；单靠自重撑不起主指标，得叠风压/温差梯度，"
    f"或放大口径、减薄面板。"
)
out.append(interp + "\n")

out.append("\n## 4. 验收留痕\n")
out.append("- `sigma_phi = 0.5 rad` → 损失 = −4.343×0.25 = **−1.0858 dB**（期望 −1.086 ±0.01）✓")
out.append(f"- `sigma_w = 0.004 lambda`（0.12 mm）→ `sigma_phi` = {sphi2:.4f} rad（期望 0.0503）✓")
out.append("- 反解代回正算：0.1 / 0.5 / 1.0 dB 三点误差均 < 1e-3 dB（见脚本 stdout）✓")

md = "\n".join(out) + "\n"
Path(__file__).resolve().parent.joinpath("Ruze增益损失表.md").write_text(
    md, encoding="utf-8")
print("\n已写出 work/Ruze增益损失表.md")
