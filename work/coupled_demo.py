"""双域耦合的第一条完整闭环：结构变形 → 口径相位 → 远场增益。

用**解析的圆板挠度**（不是 FEM，FEM 那边还在 ISSUE-017），把变形折成口径相位，
喂给现有的 PO 积分器，看增益怎么掉，再对 Ruze。

这样做的意义：把「结构 → 电磁」这条路先整段跑通、量级对上，
以后把解析挠度换成 FEniCSx 的变形场时，下游一行都不用改。

跑法（项目根目录下）：

    .\.venv\Scripts\python.exe work\coupled_demo.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

import numpy as np

from go_po.aperture_po import axial_amplitude, far_field_pattern, tapered_aperture

# ---- 几何与物理 -------------------------------------------------------
R_NORM = 10.0              # 口径半径（λ）
EDGE_TAPER_DB = -1.83      # q=8 时双反射面的实际边缘锥削
LAMBDA_M = 0.030           # 10 GHz
E = 70e9                   # 铝
NU = 0.3
RHO = 2700.0
G = 9.81

THETA_DEG = np.linspace(0.0, 8.0, 321)


def plate_deflection(r_m: np.ndarray, a_m: float, t_m: float, load_pa: float) -> np.ndarray:
    """简支圆板在均布载荷下的挠度（解析解）。

        w(r) = q/(64D) · (a² - r²) · [ (5+ν)/(1+ν)·a² - r² ]

    中心处退化为 ``q a⁴/(64D)·(5+ν)/(1+ν)``，可直接核对。
    """
    D = E * t_m**3 / (12.0 * (1.0 - NU**2))
    a2 = a_m**2
    return load_pa / (64.0 * D) * (a2 - r_m**2) * (
        (5.0 + NU) / (1.0 + NU) * a2 - r_m**2
    )


def main() -> int:
    a_m = R_NORM * LAMBDA_M           # 300 mm

    print(f"{'厚度 (mm)':>9} {'w_max (μm)':>11} {'σ_w/λ':>9} "
          f"{'σ_φ (rad)':>10} {'PO 损失':>9} {'Ruze':>9}")
    print("-" * 64)

    for t_mm in (2.0, 5.0, 10.0, 20.0):
        t_m = t_mm / 1000.0
        load = RHO * G * t_m          # 自重面载荷

        # 口径采样点（λ 单位）
        af0 = tapered_aperture(R_NORM, edge_taper_db=EDGE_TAPER_DB, n_per_lambda=12.0)
        ref = axial_amplitude(af0)
        r_norm = np.hypot(af0.x, af0.y)
        r_m = r_norm * LAMBDA_M

        w_m = plate_deflection(r_m, a_m, t_m, load)
        sigma_w_lambda = float(np.std(w_m)) / LAMBDA_M

        # 变形 → 口径相位：dφ = 2kw = 4π w/λ
        dphi = 4.0 * np.pi * (w_m / LAMBDA_M)

        af = tapered_aperture(
            R_NORM, edge_taper_db=EDGE_TAPER_DB, phase_error=dphi,
            n_per_lambda=12.0,
        )
        g = float(far_field_pattern(af, np.array([0.0]), normalize=ref)[0])
        loss_po = 20.0 * np.log10(max(g, 1e-12))

        sigma_phi = 4.0 * np.pi * sigma_w_lambda
        loss_ruze = 10.0 * np.log10(np.exp(-(sigma_phi**2)))

        print(f"{t_mm:9.1f} {np.max(np.abs(w_m))*1e6:11.2f} {sigma_w_lambda:9.5f} "
              f"{sigma_phi:10.4f} {loss_po:8.3f}dB {loss_ruze:8.3f}dB")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
