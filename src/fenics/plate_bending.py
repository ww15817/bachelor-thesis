"""薄板弯曲：用两个串联泊松方程解 Kirchhoff 板，对 Navier 解析级数校验。

背景
----
反射面的面板可以当薄板。弯曲方程是四阶的 ``∇⁴w = q/D``，直接解需要 C¹ 连续单元，
而 FEniCSx 没有现成的板/壳单元（这是 ISSUE-004）。办法是拆成两个二阶方程（混合形式）：

    M = -D ∇²w          （弯矩和）
    ∇²M = -q            （平衡）

两个都是泊松，C0 拉格朗日单元就能解。简支边的边界条件是 ``w = 0`` 且 ``M = 0``，
所以两个方程用同一套 Dirichlet 条件。

校验
----
方形简支板（边长 2a）中心挠度的 Navier 级数解：

    w = 16 q a⁴ / (π⁶ D) · Σ_{m,n odd} 1 / (m n (m² + n²)²)

薄板理论只在 ``t/a`` 很小时成立；这里顺带算一下 ``t/a = 1/20`` 的厚板偏离量。

跑法（WSL 里的 conda 环境）：

    wsl -d Ubuntu-22.04 -e bash -lc "source ~/miniforge3/etc/profile.d/conda.sh && \
        conda activate fenicsx && python /mnt/c/.../src/fenics/plate_bending.py"
"""

import numpy as np
from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem
from mpi4py import MPI
import ufl


def navier_center(a: float, D: float, q: float, terms: int = 199) -> float:
    """方形简支板中心挠度的 Navier 级数解。

    板边长 ``L = 2a``（这里 ``a`` 是半边长），中心挠度：

        w = 16 q L⁴ / (π⁶ D) · Σ_{m,n odd} 1 / (m n (m² + n²)²)

    注意 L⁴ 不是 a⁴——写成 a⁴ 会少 16 倍。
    """
    L = 2.0 * a
    s = 0.0
    for m in range(1, terms + 1, 2):
        for n in range(1, terms + 1, 2):
            s += 1.0 / (m * n * (m * m + n * n) ** 2)
    return 16.0 * q * L**4 / (np.pi**6 * D) * s


def solve_plate(a: float, t: float, E: float, nu: float, q: float, n_cell: int):
    """解 2a×2a 简支方板，返回 (w_center, 网格上的 w 函数, D)。"""
    D = E * t**3 / (12.0 * (1.0 - nu**2))

    msh = mesh.create_rectangle(
        MPI.COMM_WORLD, [[-a, -a], [a, a]], [n_cell, n_cell],
        cell_type=mesh.CellType.triangle,
    )
    V = fem.functionspace(msh, ("Lagrange", 1))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    def on_edge(x):
        return np.isclose(np.abs(x[0]), a) | np.isclose(np.abs(x[1]), a)

    fdim = msh.topology.dim - 1
    facets = mesh.locate_entities_boundary(msh, fdim, on_edge)
    dofs = fem.locate_dofs_topological(V, fdim, facets)
    bc = fem.dirichletbc(fem.Constant(msh, 0.0), dofs, V)

    bilinear = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx

    # 第一步：∇²M = -q
    M_h = LinearProblem(
        bilinear, fem.Constant(msh, q) * v * ufl.dx, bcs=[bc],
        petsc_options_prefix="plate_M",
    ).solve()

    # 第二步：∇²w = -M/D
    w_h = LinearProblem(
        bilinear, M_h / D * v * ufl.dx, bcs=[bc],
        petsc_options_prefix="plate_w",
    ).solve()

    # 中心点取值：不依赖 eval（0.11 里对任意点取值容易定位失败），
    # 直接取离 (0,0) 最近的节点值。
    coords = V.tabulate_dof_coordinates()
    radial = np.hypot(coords[:, 0], coords[:, 1])
    nearest = int(np.argmin(radial))
    w_center = float(w_h.x.array[nearest])
    return w_center, D


def export_deflection(out_path: str, a: float = 0.3, t: float = 0.002,
                      E: float = 70e9, nu: float = 0.3, n_cell: int = 80) -> str:
    """解一次板，把变形场导成 ``x, y, w``（米）存 npz，供耦合链路读取。

    这样「结构 → 相位 → 远场」里的第一步就可以换成有限元的真实变形场，
    而不是解析挠度。注意 ISSUE-017：简支边条件用的是近似，w 偏低约 5.4%。
    """
    q = 2700.0 * 9.81 * t          # 自重面载荷
    w_center, D = solve_plate(a, t, E, nu, q, n_cell)

    # 重新解一次只为拿到函数值（solve_plate 只返回中心点）
    msh = mesh.create_rectangle(
        MPI.COMM_WORLD, [[-a, -a], [a, a]], [n_cell, n_cell],
        cell_type=mesh.CellType.triangle,
    )
    V = fem.functionspace(msh, ("Lagrange", 1))
    u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

    def on_edge(x):
        return np.isclose(np.abs(x[0]), a) | np.isclose(np.abs(x[1]), a)

    fdim = msh.topology.dim - 1
    facets = mesh.locate_entities_boundary(msh, fdim, on_edge)
    dofs = fem.locate_dofs_topological(V, fdim, facets)
    bc = fem.dirichletbc(fem.Constant(msh, 0.0), dofs, V)
    bilinear = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx

    M_h = LinearProblem(
        bilinear, fem.Constant(msh, q) * v * ufl.dx, bcs=[bc],
        petsc_options_prefix="export_M",
    ).solve()
    w_h = LinearProblem(
        bilinear, M_h / D * v * ufl.dx, bcs=[bc],
        petsc_options_prefix="export_w",
    ).solve()

    coords = V.tabulate_dof_coordinates()
    np.savez(out_path, x=coords[:, 0], y=coords[:, 1],
             w=w_h.x.array, a=a, t=t, D=D, w_center=w_center)
    return out_path


def main() -> int:
    # 600 mm 口径的反射面板：半径 300 mm，铝，2 mm 厚，只算自重
    a = 0.3
    E = 70e9
    nu = 0.3
    rho = 2700.0
    g = 9.81

    print(f"{'t (mm)':>7} {'t/a':>7} {'网格':>6} {'数值 w':>12} {'Navier':>12} {'相对差':>9}")
    print("-" * 62)
    for t_mm in (2.0, 5.0, 15.0):
        t = t_mm / 1000.0
        q = rho * g * t                     # 自重面载荷 Pa
        for n_cell in (40, 80):
            w_num, D = solve_plate(a, t, E, nu, q, n_cell)
            w_ref = navier_center(a, D, q)
            rel = (w_num - w_ref) / w_ref
            print(f"{t_mm:7.1f} {t/a:7.4f} {n_cell:6d} "
                  f"{w_num*1e6:11.4f}μm {w_ref*1e6:11.4f}μm {rel:+9.2%}")
    # 导出一份变形场供耦合链路用
    out = export_deflection("data/feeds/plate_deflection_2mm.npz", t=0.002)
    print(f"变形场已导出: {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
