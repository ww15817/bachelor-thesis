"""openEMS 运行环境引导（Windows 原生，无需 WSL）。

openEMS 的 Python 绑定（CSXCAD / openEMS 两个 wheel）依赖同目录下的
CSXCAD.dll 等二进制。Python 3.8 之后 PATH 不再参与扩展模块的 DLL 搜索，
必须显式 add_dll_directory，因此统一在这里做一次。

用法：脚本开头 `import _env` 即可，需先设环境变量 OPENEMS_ROOT 指向
解压后的 openEMS 目录（含 openEMS.exe / CSXCAD.dll）。

    set OPENEMS_ROOT=D:\\openEMS
    python horn_pinfeed.py

如果没装 openEMS，import 时会给出明确的提示，而不是一句
"DLL load failed"。
"""
import os
import sys

EMS_ROOT = os.environ.get("OPENEMS_ROOT", "")

if EMS_ROOT and os.path.isdir(EMS_ROOT):
    os.add_dll_directory(EMS_ROOT)
    os.environ["PATH"] = EMS_ROOT + os.pathsep + os.environ.get("PATH", "")
    os.environ.setdefault("OPENEMS_BIN", os.path.join(EMS_ROOT, "openEMS.exe"))
else:
    sys.stderr.write(
        "[_env] 未设置 OPENEMS_ROOT 或目录不存在。\n"
        "       请从 https://github.com/thliebig/openEMS-Project/releases 下载\n"
        "       openEMS_x64_v0.37.0_msvc.zip，解压后把 OPENEMS_ROOT 指向里面的\n"
        "       openEMS 目录，并 pip install 该目录 python/ 下的两个 wheel。\n"
    )

try:
    import matplotlib
    matplotlib.use("Agg")     # 无显示环境下画图不弹窗
except ImportError:
    pass
