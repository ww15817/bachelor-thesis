"""圆板弯曲：用 gmsh 建圆盘，双泊松串联解 Kirchhoff 板，对解析解校验。

为什么要有这个
--------------
`plate_bending.py` 用的是**方板**（`create_rectangle`），而反射面是**圆口径**。
方板和圆板的挠度不一样，耦合链路里直接拿方板去对圆口径的相位，会引一个
几何形状的偏差进去。这里补上圆板，让结构侧和射频侧的口径一致。

解析解（简支圆板、均布载荷）：

    w(r) = q/(64D) · (a² − r²) · [ (5+ν)/(1+ν)·a² − r² ]
    w_max = q a⁴/(64D) · (5+ν)/(1+ν)

同 `plate_bending.py` 一样，简支边施的是 `M=0`（弯矩和）而不是 `M_nn=0`，
所以会偏低一点（ISSUE-017），这里顺带把偏差量出来。

跑法（WSL 里的 conda 环境）：

    wsl -d Ubuntu-22.04 -e bash -lc "source ~/miniforge3/etc/profile.d/conda.sh && \
        conda activate fenicsx && python /mnt/c/.../src/fenics/plate_disk.py"
"""

import numpy as np
from dolfinx import fem, mesh
from dolfinx.fem.petsc import LinearProblem
from dolfinx.io import gmsh as gmshio
from mpi4py import MPI
import gmsh
import ufl


def make_disk(a: float, h: float, tag: str = "disk"):
    """用 gmsh 建圆盘网格，读进 dolfinx。"""
    gmsh.initialize()
    gmsh.model.add(tag)
    gmsh.model.occ.addDisk(0.0, 0.0, 0.0, a, a)
    gmsh.model.occ.synchronize()
    # dolfinx 的 gmsh 读取器要求有 physical groups，不然报
    # "No 'physical groups' in gmsh mesh"
    surfaces = gmsh.model.getEntities(2)
    gmsh.model.addPhysicalGroup(2, [s[1] for s in surfaces], 1)
    curves = gmsh.model.getEntities(1)
    gmsh.model.addPhysicalGroup(1, [c[1] for c in curves], 2)
    gmsh.option.setNumber("Mesh.MeshSizeMax", h)
    gmsh.model.mesh.generate(2)
    path = f"/tmp/{tag}.msh"
    gmsh.write(path)
    gmsh.finalize()

    md = gmshio.read_from_msh(path, MPI.COMM_WORLD, rank=0, gdim=2)
    return md.mesh, md.facet_tags


def solid_simply_supported_w(a, t, E, nu, q):
    """简支圆板中心挠度的解析值。"""
    D = E * t**3 / (12.0 * (1.0 - nu**2))
    return q * a**4 / (64.0 * D) * (5.0 + nu) / (1.0 + nu), D


def solve(a=0.3, t=0.002, E=70e9, nu=0.3, h=0.02, n_runs=2):
    q = 2700.0 * 9.81 * t
    w_ref, D = solid_simply_supported_w(a, t, E, nu, q)

    print(f"{'网格尺寸':>10} {'节点数':>8} {'数值 w':>12} {'解析 w':>12} {'相对差':>9}")
    print("-" * 58)
    for i in range(n_runs):
        msh, _ = make_disk(a, h / (2 ** i), tag=f"disk{i}")

        V = fem.functionspace(msh, ("Lagrange", 1))
        u, v = ufl.TrialFunction(V), ufl.TestFunction(V)

        # 圆盘边界：所有外边界点
        fdim = msh.topology.dim - 1
        facets = mesh.exterior_facet_indices(msh.topology)
        dofs = fem.locate_dofs_topological(V, fdim, facets)
        bc = fem.dirichletbc(fem.Constant(msh, 0.0), dofs, V)

        bilinear = ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx
        M_h = LinearProblem(bilinear, fem.Constant(msh, q) * v * ufl.dx, bcs=[bc],
                            petsc_options_prefix=f"diskM{i}").solve()
        w_h = LinearProblem(bilinear, M_h / D * v * ufl.dx, bcs=[bc],
                            petsc_options_prefix=f"diskW{i}").solve()

        coords = V.tabulate_dof_coordinates()
        radial = np.hypot(coords[:, 0], coords[:, 1])
        w_center = float(w_h.x.array[int(np.argmin(radial))])
        n_nodes = V.dofmap.index_map.size_local
        rel = (w_center - w_ref) / w_ref
        print(f"{h/(2**i):10.4f} {n_nodes:8d} {w_center*1e6:11.3f}μm "
              f"{w_ref*1e6:11.3f}μm {rel:+9.2%}")
    return 0


def main() -> int:
    return solve()


if __name__ == "__main__":
    raise SystemExit(main())
