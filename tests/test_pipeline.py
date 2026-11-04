"""端到端流水线：配置文件进，结果文件出。"""

import json

import numpy as np

from antennasim.pipeline import run_pipeline


def test_pipeline_runs_and_reports_ok(cfg, tmp_path):
    result = run_pipeline(cfg, project_root=tmp_path, make_figure=False, verbose=False)
    assert result.ok, f"流水线未通过校验: {result.metrics}"
    assert result.seconds > 0


def test_pipeline_writes_outputs(cfg, tmp_path):
    result = run_pipeline(cfg, project_root=tmp_path, make_figure=False, verbose=False)
    npz_path = tmp_path / cfg.output.directory / f"{cfg.output.tag}.npz"
    json_path = tmp_path / cfg.output.directory / f"{cfg.output.tag}.json"
    assert npz_path.exists()
    assert json_path.exists()

    meta = json.loads(json_path.read_text(encoding="utf-8"))
    assert meta["case"] == cfg.name
    assert all(meta["checks"].values()), meta["checks"]


def test_pipeline_gain_accounts_for_blockage(cfg, tmp_path):
    """遮挡变大，峰值增益必须按被遮掉的口径功率下降。

    旧版跑两遍比较 metrics 字典：``run_pipeline`` 是纯确定性函数（没有随机量、
    没有时间戳进指标），两遍必然相同，等于没测。这里换物理判据：

    * 视轴增益损失应等于均匀口径的解析值 ``20log10(1-b)``，b = 遮挡面积比。
      实际照明有锥削，但遮挡区（L=6.4 时 rho<2.2λ）内幅度只比中心低约 1.5 dB，
      对损失的修正远小于 0.1 dB（实测 0.041 dB）。
    * 把副面顶点从 6.4λ 移到 6.0λ，遮挡 4.9% -> 7.7%，增益必须跟着变差。
      把遮挡归一化掉、或返回常数指标的实现过不了这两关。
    """
    import dataclasses

    base = run_pipeline(cfg, project_root=tmp_path, make_figure=False, verbose=False)
    blocked_cfg = dataclasses.replace(
        cfg, geometry=dataclasses.replace(cfg.geometry, sub_vertex=6.0)
    )
    blocked = run_pipeline(
        blocked_cfg, project_root=tmp_path, make_figure=False, verbose=False
    )
    assert base.ok and blocked.ok

    ratio_base = float(base.metrics["遮挡面积比"].split()[0]) / 100.0
    ratio_blocked = float(blocked.metrics["遮挡面积比"].split()[0]) / 100.0
    loss_base = float(base.metrics["峰值增益损失"].split()[0])
    loss_blocked = float(blocked.metrics["峰值增益损失"].split()[0])

    assert ratio_blocked > ratio_base

    for ratio, loss in ((ratio_base, loss_base), (ratio_blocked, loss_blocked)):
        expected = 20.0 * np.log10(1.0 - ratio)
        assert abs(loss - expected) < 0.1, (
            f"遮挡损失 {loss:.3f} dB 与解析值 {expected:.3f} dB 不符"
        )

    assert loss_blocked < loss_base


def test_pipeline_rejects_inconsistent_geometry(cfg, tmp_path):
    """把副面位置改到非法范围，应在建模阶段就报错。"""
    import dataclasses

    import pytest

    from antennasim.config import ConfigError

    bad_geom = dataclasses.replace(cfg.geometry, sub_vertex=20.0)
    bad_cfg = dataclasses.replace(cfg, geometry=bad_geom)
    with pytest.raises(ConfigError):
        bad_cfg.validate()
