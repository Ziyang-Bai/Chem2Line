"""向后兼容层 - 已弃用

此模块保留仅为兼容旧代码。所有功能已迁移至：
- lib.config: 配置管理
- lib.i18n: 国际化
- lib.history: 历史记录
- lib.chemistry: 化学计算
- lib.database: 数据库操作
- lib.viewer3d: 3D 查看器
- lib.version: 版本信息

请使用 `from lib import ...` 导入所需功能。
"""

import warnings

warnings.warn(
    "lib.engine 已弃用，请使用 lib 包中的各子模块",
    DeprecationWarning,
    stacklevel=2,
)

# 向后兼容导出
from lib.version import VERSION, DEVELOPER, DATE, CORE_VERSION
from lib.config import load_config, save_config
from lib.i18n import load_language
from lib.database import load_smiles_database, get_smiles_options, get_database_info
from lib.chemistry import formula_to_bondline as formula_to_image, analyze_molecule
from lib.history import load_history, save_history, add_to_history
from lib.viewer3d import MoleculeViewer

DATABASE_PATH = "lib/db/default_database.xml"
