"""向后兼容层 - 已弃用

此模块保留仅为兼容旧代码。所有功能已迁移至：
- lib.chemistry: 化学计算
- lib.database: 数据库操作
- lib.viewer3d: 3D 查看器
- lib.version: 版本信息

请使用 `from lib import ...` 导入所需功能。
"""

import warnings

warnings.warn(
    "lib.ctlcore 已弃用，请使用 lib.chemistry / lib.database / lib.viewer3d",
    DeprecationWarning,
    stacklevel=2,
)

# 向后兼容导出
from lib.version import CORE_VERSION, VERSION
from lib.database import load_smiles_database, get_smiles_options, get_database_info
from lib.chemistry import (
    formula_to_bondline,
    overlay_force_field,
    analyze_molecule,
    get_chemical_info,
)
from lib.viewer3d import MoleculeViewer


def core_version():
    return CORE_VERSION


def show_chemical_info(smiles):
    info = get_chemical_info(smiles)
    return "\n".join(f"{k}: {v}" for k, v in info.items())
