"""扫表面误差，看 20λ 口径下能看出什么。

回答 ISSUE-003：面板变形要多大，增益才看得出来掉、副瓣才明显抬高。

这是**参数化研究**：给定表面误差 rms（以波长为单位），看远场怎么变。
不绑定具体结构（自重/风载/温度各算各的），因为 20λ 口径配自重的量级已经算过，
只有 0.01 dB 级，撑不起主指标。

用法（项目根目录下）：

    .\.venv\Scripts\python.exe work\scan_surface_error.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np

from go_po.aperture_po import axial_amplitude, far_field_pattern, tapered_aperture

R = 10.0                     # 20λ 口径半径（λ 归一化）
K = 2.0 * np.pi
EDGE_TAPER_DB = -1.83        # q=8 时双反射面口径的实际边缘锥削
N_PER_LAMBDA = 12.0
TRIALS = 16
SIGMA_LIST = (0.0, 0.002, 0.005, 0.010, 0.020, 0.040, 0.060)

THETA_DEG = np.linspace(0.0, 12.0, 481)


def first_sidelobe_db(pat_db: np.ndarray) -> tuple[float, float]:
    """第一副瓣电平（dB）和它的位置（度）。"""
    i0 = int(np.argmax(pat_db))
    j = i0
    while j < len(pat_db) - 1 and pat_db[j + 1] <= pat_db[j]:
        j += 1
    k2 = j
    while k2 < len(pat_db) - 1 and pat_db[k2 + 1] >= pat_db[k2]:
        k2 += 1
    return float(pat_db[k2]), float(THETA_DEG[k2])


def run() -> int:
    rng = np.random.default_rng(20261104)

    base = tapered_aperture(R, edge_taper_db=EDGE_TAPER_DB, n_per_lambda=N_PER_LAMBDA)
    ref = axial_amplitude(base)
    r = np.hypot(base.x, base.y)

    rows = []
    for sigma in SIGMA_LIST:
        if sigma == 0.0:
            pat_db = 20.0 * np.log10(np.clip(
                far_field_pattern(base, THETA_DEG, normalize=ref), 1e-12, None))
            gain_db = 0.0
            sll_db, sll_deg = first_sidelobe_db(pat_db)
            rows.append((sigma, 0.0, gain_db, sll_db, sll_deg))
            continue

        gains, slls, degs = [], [], []
        for _ in range(TRIALS):
            dphi = 2.0 * K * sigma * rng.normal(size=r.shape)
            af = tapered_aperture(
                R, edge_taper_db=EDGE_TAPER_DB, phase_error=dphi,
                n_per_lambda=N_PER_LAMBDA,
            )
            g = float(far_field_pattern(af, np.array([0.0]), normalize=ref)[0])
            gains.append(20.0 * np.log10(max(g, 1e-12)))

            pat_db = 20.0 * np.log10(np.clip(
                far_field_pattern(af, THETA_DEG, normalize=ref), 1e-12, None))
            sll, deg = first_sidelobe_db(pat_db)
            slls.append(sll)
            degs.append(deg)

        rows.append((sigma, 2.0 * K * sigma, float(np.mean(gains)),
                     float(np.mean(slls)), float(np.mean(degs))))

    print(f"{'sigma/lam':>10} {'sigma_phi':>10} {'gain_dB':>9} {'SLL_dB':>8} {'SLL_deg':>8}")
    print("-" * 52)
    for sigma, sphi, gain, sll, deg in rows:
        print(f"{sigma:10.4f} {sphi:10.4f} {gain:9.3f} {sll:8.2f} {deg:8.2f}")

    print()
    # 反解：要掉 0.5 / 1.0 dB，需要多大的 sigma
    for target in (0.5, 1.0, 2.0):
        sphi_needed = np.sqrt(target / 4.343)
        sigma_needed = sphi_needed / (2.0 * K)
        print(f"掉 {target:.1f} dB 需要 sigma = {sigma_needed:.4f} lam "
              f"({sigma_needed * 30.0:.3f} mm @10GHz)")

    return 0


if __name__ == "__main__":
    raise SystemExit(run())
