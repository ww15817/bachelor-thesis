r"""扫馈源 q，看总增益怎么变。

回答 2026-11-04 外部审查提的问题：q=8 时馈源波束比副面张角宽，64% 功率溢出。
把 q 加大能让溢出降下来，但边缘锥削变深会损失口径效率——要找两者平衡的那一点。

归一约定按审查 C-02 的方案 A（方向性 D = 4πU_max/P_rad，被挡的功率仍算辐射），
所以最优 q 可能落在扫描上界，Q_LIST 开到 96 才能确认是极大还是没扫到头。

用法（项目根目录下）：

    .\.venv\Scripts\python.exe work\scan_q.py
"""

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np

from antennasim.config import load_config
from antennasim.pipeline import run_pipeline

Q_LIST = (4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 47.0, 64.0, 80.0, 96.0)


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    base = load_config(root / "cases" / "cassegrain_20lambda.yaml")

    print(f"{'q':>5} {'溢出':>7} {'照明':>7} {'总效率':>8} {'绝对增益':>9}")
    print("-" * 42)
    rows = []
    for q in Q_LIST:
        # 只改 q，保留 feed 里其余字段（审查 C-10：直接替换 FeedConfig 会把
        # edge_taper_db 等重置成默认值，绕过流水线里的偏差告警）
        cfg = dataclasses.replace(
            base, feed=dataclasses.replace(base.feed, q=q)
        )
        res = run_pipeline(cfg, project_root=root, make_figure=False, verbose=False)
        m = res.metrics
        spill = float(m["溢出发射效率"])
        illum = float(m["照明效率"])
        total = float(m["总效率"])
        gain = float(m["绝对增益"].split()[0])
        rows.append((q, spill, illum, total, gain))
        print(f"{q:5.1f} {spill:7.4f} {illum:7.4f} {total:8.4f} {gain:9.2f}")

    best = max(rows, key=lambda r: r[4])
    # 基线从本次扫描里取，不写死字面量（审查 C-10：改了口径这里会自动跟着变）
    base_q8 = next(g for (qq, _s, _i, _t, g) in rows if qq == 8.0)
    print()
    print(f"最优 q = {best[0]:g}，总效率 {best[3]:.4f}，绝对增益 {best[4]:.2f} dBi")
    print(f"相对 q=8（{base_q8:.2f} dBi）净增 {best[4] - base_q8:+.2f} dB")
    if best[0] == Q_LIST[-1]:
        print("⚠️ 最优落在扫描上界，还不能确认是极大——把 Q_LIST 再往上开")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
