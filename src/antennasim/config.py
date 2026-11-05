"""配置文件读取与校验。

把算例参数从代码里抽出来，放在 YAML 文件里。改参数不需要改代码。
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

import yaml

__all__ = ["CaseConfig", "load_config", "ConfigError"]


class ConfigError(ValueError):
    """配置文件缺字段或取值非法。"""


@dataclass(frozen=True)
class GeometryConfig:
    diameter: float = 20.0
    f_over_d: float = 0.4
    magnification: float = 4.0
    sub_vertex: float | None = None   # None 表示取默认 M*F/(M+1)


@dataclass(frozen=True)
class FeedConfig:
    q: float = 8.0
    # 期望的口径边缘锥削（dB），可为 None。
    #
    # ⚠️ 它**不参与计算**：边缘锥削是结果，由 q、几何和双反射面的扩散因子共同决定
    # （q=8、20λ 算例的实际值是约 -1.9 dB，不是 -10 dB）。
    # 流水线会把实际算出的值跟它比，差超过 1 dB 就告警；
    # 想真做到某个锥削值，应该去反解 q，而不是在这里写一个数。
    edge_taper_db: float | None = None

    # ---- 馈源来源（P2，2026-11-05）--------------------------------------
    # analytic：用解析 cos^q（q 决定波束宽度），原来的行为
    # npz     ：从 openEMS 全波算出的方向图文件里读（需要 file 与 plane）
    source: str = "analytic"
    file: str | None = None
    # 取哪个面的方向图：E 面是 phi=90°（TE10 的 E 沿 y），H 面是 phi=0°（x-z 面）。
    # ⚠️ 矩形喇叭两个面的锥削差很多（ISSUE-020：80×60 口径差 3.8 dB），
    # 而流水线假设轴对称，所以这两个面会给出明显不同的增益，别混着说。
    plane: str = "E"



@dataclass(frozen=True)
class SamplingConfig:
    n_radial: int = 26
    n_azimuth: int = 72
    n_rho: int = 200
    n_phi: int = 320


@dataclass(frozen=True)
class OutputConfig:
    directory: str = "data/results"
    tag: str = "cassegrain"
    figure: bool = True


@dataclass(frozen=True)
class CaseConfig:
    name: str
    description: str = ""
    geometry: GeometryConfig = field(default_factory=GeometryConfig)
    feed: FeedConfig = field(default_factory=FeedConfig)
    sampling: SamplingConfig = field(default_factory=SamplingConfig)
    output: OutputConfig = field(default_factory=OutputConfig)

    # ---------------------------------------------------------- 派生量
    @property
    def focal_length(self) -> float:
        return self.geometry.f_over_d * self.geometry.diameter

    @property
    def sub_vertex(self) -> float:
        """副面顶点位置；未指定时取 M*F/(M+1)。"""
        if self.geometry.sub_vertex is not None:
            return self.geometry.sub_vertex
        m = self.geometry.magnification
        return m * self.focal_length / (m + 1.0)

    def validate(self) -> None:
        g = self.geometry
        if g.diameter <= 0:
            raise ConfigError(f"geometry.diameter 必须为正，当前 {g.diameter}")
        if not 0.0 < g.f_over_d < 2.0:
            raise ConfigError(f"geometry.f_over_d 应在 (0, 2)，当前 {g.f_over_d}")
        if g.magnification <= 1.0:
            raise ConfigError(f"geometry.magnification 必须 > 1，当前 {g.magnification}")
        if not 0.0 < self.sub_vertex < self.focal_length:
            raise ConfigError(
                f"副面顶点须满足 0 < L < F({self.focal_length:.3f})，"
                f"当前 L={self.sub_vertex:.3f}"
            )
        s = self.sampling
        if s.n_radial < 4 or s.n_azimuth < 8:
            raise ConfigError("采样密度过低：n_radial >= 4，n_azimuth >= 8")
        if s.n_rho < 16 or s.n_phi < 16:
            raise ConfigError("口径采样过低：n_rho >= 16，n_phi >= 16")

    def as_dict(self) -> dict:
        """导出为与配置文件一致的结构，可原样写回 YAML。"""
        return {
            "case": {"name": self.name, "description": self.description},
            "geometry": asdict(self.geometry),
            "feed": asdict(self.feed),
            "sampling": asdict(self.sampling),
            "output": asdict(self.output),
        }

    def dump(self, path: str | Path) -> Path:
        """把当前配置写回 YAML 文件（改完参数可直接存档）。"""
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(
            yaml.safe_dump(self.as_dict(), allow_unicode=True, sort_keys=False),
            encoding="utf-8",
        )
        return p


def _section(raw: dict, key: str) -> dict:
    value = raw.get(key, {})
    if not isinstance(value, dict):
        raise ConfigError(f"配置项 {key} 必须是映射（键值对），当前是 {type(value).__name__}")
    return value


def load_config(path: str | Path) -> CaseConfig:
    """读取 YAML 配置文件并校验。"""
    p = Path(path)
    if not p.exists():
        raise ConfigError(f"找不到配置文件：{p}")

    try:
        raw = yaml.safe_load(p.read_text(encoding="utf-8")) or {}
    except yaml.YAMLError as exc:
        raise ConfigError(f"YAML 解析失败：{exc}") from exc

    if not isinstance(raw, dict):
        raise ConfigError("配置文件顶层必须是键值对")

    case = _section(raw, "case")
    name = case.get("name")
    if not name:
        raise ConfigError("配置项 case.name 缺失")

    try:
        cfg = CaseConfig(
            name=str(name),
            description=str(case.get("description", "")),
            geometry=GeometryConfig(**_section(raw, "geometry")),
            feed=FeedConfig(**_section(raw, "feed")),
            sampling=SamplingConfig(**_section(raw, "sampling")),
            output=OutputConfig(**_section(raw, "output")),
        )
    except TypeError as exc:
        # 多打一个字（比如 feed: {qq: 8}）以前会直接 traceback 出来，
        # 看起来像程序坏了；统一转成 ConfigError，命令行会报「配置错误」并返回 2。
        raise ConfigError(f"配置里有无法识别的字段：{exc}") from exc
    cfg.validate()
    return cfg
