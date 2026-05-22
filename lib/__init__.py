"""Chem2Line 核心库

提供化学式→键线式转换、数据库管理、3D 分子查看等功能。
"""

from lib.version import VERSION, CORE_VERSION, DEVELOPER, DATE, APP_NAME, REPOSITORY
from lib.config import load_config, save_config
from lib.i18n import load_language, get_text
from lib.database import load_smiles_database, get_smiles_options, get_database_info
from lib.chemistry import (
    formula_to_bondline,
    formula_to_svg,
    overlay_force_field,
    analyze_molecule,
    get_chemical_info,
)
from lib.history import load_history, save_history, add_to_history, clear_history
from lib.viewer3d import MoleculeViewer

__all__ = [
    # 版本
    "VERSION", "CORE_VERSION", "DEVELOPER", "DATE", "APP_NAME", "REPOSITORY",
    # 配置
    "load_config", "save_config",
    # 国际化
    "load_language", "get_text",
    # 数据库
    "load_smiles_database", "get_smiles_options", "get_database_info",
    # 化学
    "formula_to_bondline", "formula_to_svg", "overlay_force_field",
    "analyze_molecule", "get_chemical_info",
    # 历史
    "load_history", "save_history", "add_to_history", "clear_history",
    # 3D 查看器
    "MoleculeViewer",
]
