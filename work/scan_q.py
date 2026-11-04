"""扫馈源 q，看总增益怎么变。

回答 2026-11-04 外部审查提的问题：q=8 时馈源波束比副面张角宽，64% 功率溢出。
把 q 加大能让溢出降下来，但边缘锥削变深会损失照明效率——要找到两者平衡的那一点。

用法（项目根目录下）：

    .\.venv\Scripts\python.exe work\scan_q.py
"""

import dataclasses
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np

from antennasim.config import FeedConfig, load_config
from antennasim.pipeline import run_pipeline

Q_LIST = (4.0, 6.0, 8.0, 12.0, 16.0, 24.0, 32.0, 47.0, 64.0)


def main() -> int:
    root = Path(__file__).resolve().parent.parent
    base = load_config(root / "cases" / "cassegrain_20lambda.yaml")

    print(f"{'q':>5} {'溢出':>7} {'照明':>7} {'总效率':>8} {'绝对增益':>9}")
    print("-" * 42)
    rows = []
    for q in Q_LIST:
        cfg = dataclasses.replace(base, feed=FeedConfig(q=q))
        res = run_pipeline(cfg, project_root=root, make_figure=False, verbose=False)
        m = res.metrics
        spill = float(m["溢出发射效率"])
        illum = float(m["照明效率"])
        total = float(m["总效率"])
        gain = float(m["绝对增益"].split()[0])
        rows.append((q, spill, illum, total, gain))
        print(f"{q:5.1f} {spill:7.4f} {illum:7.4f} {total:8.4f} {gain:9.2f}")

    best = max(rows, key=lambda r: r[4])
    print()
    print(f"最优 q = {best[0]:g}，总效率 {best[3]:.4f}，绝对增益 {best[4]:.2f} dBi")
    print(f"相对 q=8（31.24 dBi）净增 {best[4] - 31.24:+.2f} dB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
