# engine.py
import os
import sys
import time
import datetime
import logging
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional, Union

# 基础导入
import numpy as np
import tkinter as tk
from PIL import Image
from tkinter import messagebox, Canvas, Frame

# 尝试导入RDKit
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
except ImportError as e:
    logging.error(f"无法导入RDKit: {e}")
    Chem = None
    AllChem = None
    Draw = None

# 检查RDKit是否可用
if None in (Chem, AllChem, Draw):
    raise ImportError("RDKit未正确安装，请确保安装了RDKit及其所有依赖")

__all__ = [
    'VERSION', 'DEVELOPER', 'DATE', 'CORE_VERSION', 'DATABASE_PATH',
    'load_config', 'save_config', 'load_language',
    'load_history', 'save_history', 'add_to_history',
    'get_smiles_options', 'formula_to_image',
    'analyze_molecule', 'show_3d_viewer',
    'search_smiles_database', 'convert_sdf_to_smiles',
    'MoleculeViewer'
]

# 版本信息
VERSION = "1.4.1"
DEVELOPER = "Ziyang-Bai"
DATE = datetime.date(2025, 6, 9).strftime("%Y-%m-%d")

# 路径常量
HISTORY_FILE = "lib/history.xml"
DATABASE_PATH = "lib/db/default_database.xml"
CONFIG_FILE = "lib/config/config.xml"
CORE_VERSION = VERSION

VERSION = "1.4.1"
DEVELOPER = "Ziyang-Bai"
DATE = "2025-01-01"

# 全局变量
HISTORY_FILE = "lib/history.xml"
DATABASE_PATH = "lib/db/default_database.xml"

def get_version():
    """返回软件版本号"""
    return VERSION

def get_developer():
    """返回开发者信息"""
    return DEVELOPER

def get_date():
    """返回发布日期"""
    return DATE

def get_default_database_path():
    """返回默认数据库路径"""
    return DATABASE_PATH

def load_history():
    """
    从XML文件加载历史记录
    返回: list 历史记录列表
    """
    if not os.path.exists(HISTORY_FILE) or os.path.getsize(HISTORY_FILE) == 0:
        return []
    try:
        tree = ET.parse(HISTORY_FILE)
        root = tree.getroot()
        history = []
        for entry in root.findall("entry"):
            smiles_elem = entry.find("smiles")
            timestamp_elem = entry.find("timestamp")
            input_text_elem = entry.find("input_text")
            
            # 安全地获取文本值
            smiles = smiles_elem.text if smiles_elem is not None else None
            timestamp = timestamp_elem.text if timestamp_elem is not None else None
            input_text = input_text_elem.text if input_text_elem is not None else None
            
            if all([smiles, timestamp, input_text]):  # 确保所有值都不为None
                history.append({"smiles": smiles, "timestamp": timestamp, "input_text": input_text})
        return history
    except (ET.ParseError, AttributeError):
        return []

def save_history(history):
    """
    将历史记录保存到XML文件
    :param history: list 历史记录列表
    """
    root = ET.Element("history")
    for entry in history:
        if all(entry.get(key) for key in ["smiles", "timestamp", "input_text"]):
            entry_element = ET.SubElement(root, "entry")
            smiles_element = ET.SubElement(entry_element, "smiles")
            smiles_element.text = str(entry["smiles"])
            timestamp_element = ET.SubElement(entry_element, "timestamp")
            timestamp_element.text = str(entry["timestamp"])
            input_text_element = ET.SubElement(entry_element, "input_text")
            input_text_element.text = str(entry["input_text"])
    tree = ET.ElementTree(root)
    tree.write(HISTORY_FILE)

def add_to_history(smiles, input_text, history):
    """
    添加新记录到历史记录中
    :param smiles: str SMILES字符串
    :param input_text: str 输入文本
    :param history: list 历史记录列表
    :return: list 更新后的历史记录列表
    """
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    if smiles not in [entry["smiles"] for entry in history]:
        history.append({"smiles": smiles, "timestamp": timestamp, "input_text": input_text})
        save_history(history)
    return history

def formula_to_image(smiles, image_type="PNG"):
    """
    根据 SMILES 生成键线式图像
    :param smiles: str SMILES 字符串
    :param image_type: str 图像类型 ("PNG" 或 "SVG")
    :return: 图像对象 (PIL.Image 或 SVG 字符串)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("无效的 SMILES 字符串")

    if image_type.upper() == "SVG":
        return Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(300, 300), useSVG=True)
    elif image_type.upper() == "PNG":
        return Draw.MolToImage(mol, size=(300, 300))
    else:
        raise ValueError("不支持的图像类型，请使用 'PNG' 或 'SVG'")

def get_smiles_options(formula, smiles_dict):
    """
    获取对应化学式的所有可能 SMILES
    :param formula: str 输入的化学式或 SMILES
    :param smiles_dict: dict 化学式和 SMILES 的映射字典
    :return: list 包含所有可能 SMILES 的列表
    """
    matches = []
    for key, values in smiles_dict.items():
        if isinstance(values, list):
            if key == formula or formula in values:
                matches.extend(values)
        elif isinstance(values, dict):
            if formula == values.get("formula") or formula == values.get("name"):
                matches.append(values.get("smiles"))
            
    if not matches:
        # 如果没有找到匹配，且输入可能是 SMILES
        try:
            mol = Chem.MolFromSmiles(formula)
            if mol is not None:
                matches.append(formula)
        except:
            pass
            
    return matches

def search_smiles_database(query, smiles_dict, search_type="name"):
    """
    在数据库中查询化学物质
    :param query: str 查询字符串(名称或化学式)
    :param smiles_dict: dict SMILES数据库
    :param search_type: str 搜索类型("name"或"formula")
    :return: list 匹配的化学物质列表 [(name, formula, smiles), ...]
    """
    results = []
    query = query.lower()
    for smiles, info in smiles_dict.items():
        if isinstance(info, dict):
            if query in info.get(search_type, "").lower():
                results.append((info.get("name", "N/A"), info.get("formula", "N/A"), smiles))
        elif isinstance(info, list):
            for formula in info:
                if query in formula.lower():
                    results.append(("N/A", formula, smiles))
    return results

def convert_sdf_to_smiles(sdf_file, output_file):
    """
    将SDF文件转换为SMILES格式
    :param sdf_file: str SDF文件路径
    :param output_file: str 输出文件路径
    """
    with open(sdf_file, "r") as sdf, open(output_file, "w") as output:
        suppl = Chem.SDMolSupplier(sdf_file)
        for mol in suppl:
            if mol is not None:
                smiles = Chem.MolToSmiles(mol)
                output.write(f"{smiles}\n")

class MoleculeViewer:
    """3D分子查看器基类"""
    def __init__(self, root, smiles="C", model_type="ball_and_stick", lang_dict=None):
        self.root = root
        self.lang_dict = lang_dict or {}
        self.root.title(self.lang_dict.get("molecule_viewer_3d", "3D Molecule Viewer"))
        self.canvas = Canvas(root, width=600, height=600, bg="black")
        self.canvas.pack()

        # 视角参数
        self.angle_x = self.angle_y = self.angle_z = 0
        self.scale = 100         # 缩放范围从50到500
        self.perspective_d = 5   # 透视参数范围从1到20

        # 摄像机位置
        self.cam_offset_x = self.cam_offset_y = 0

        # 鼠标拖动状态
        self.last_x = self.last_y = 0
        self.last_pan_x = self.last_pan_y = 0

        # 模型类型
        self.model_type = model_type
        self.setup()
        self.load_molecule(smiles)

    def setup(self):
        """初始化控件和事件绑定"""
        self.create_controls()
        self.setup_mouse_controls()

    def create_controls(self):
        """创建控制按钮"""
        control_frame = tk.Frame(self.root)
        control_frame.pack(pady=10)

        buttons = [
            ("reset_view", "重置视角", self.reset_view, 0, 1),
            ("zoom_in", "缩放+", self.zoom_in, 1, 0),
            ("zoom_out", "缩小", self.zoom_out, 1, 2),
            ("perspective_increase", "透视+", lambda: self.adjust_perspective(1), 2, 0),
            ("perspective_decrease", "透视-", lambda: self.adjust_perspective(-1), 2, 2)
        ]

        for key, default, command, row, col in buttons:
            text = self.lang_dict.get(key) or default
            tk.Button(control_frame, text=text, command=command).grid(row=row, column=col, padx=5)

    def load_molecule(self, smiles):
        """加载分子"""
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("无效的SMILES字符串")
        
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        conf = mol.GetConformer()
        self.atoms = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
        self.atom_colors = [(self.get_atom_color(atom.GetSymbol()), atom.GetSymbol()) 
                           for atom in mol.GetAtoms()]
        self.bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble())
                     for bond in mol.GetBonds()]
        self.render_molecule()

    def get_atom_color(self, symbol):
        """获取原子颜色"""
        colors = {
            "H": "#FFFFFF", "C": "#00FFFF", "N": "#0000FF", 
            "O": "#FF0000", "F": "#00FF00", "Cl": "#00FF00",
            "Br": "#8B0000", "I": "#9400D3", "S": "#FFFF00", 
            "P": "#FFA500"
        }
        return colors.get(symbol, "#FF69B4")

    def get_vdw_radius(self, symbol):
        """获取范德华半径"""
        radii = {
            "H": 1.2, "C": 1.7, "N": 1.55, "O": 1.52,
            "F": 1.47, "Cl": 1.75, "Br": 1.85, "I": 1.98,
            "S": 1.8, "P": 1.8
        }
        return radii.get(symbol, 1.5)

    # MoleculeViewer类的方法
    def setup_mouse_controls(self):
        """设置鼠标控制"""
        self.canvas.bind("<Button-1>", self.start_drag)
        self.canvas.bind("<B1-Motion>", self.drag)
        self.canvas.bind("<Button-3>", self.start_pan)
        self.canvas.bind("<B3-Motion>", self.pan)
        self.canvas.bind("<MouseWheel>", self.on_zoom_wheel)
        self.canvas.bind("<Control-MouseWheel>", self.on_perspective_wheel)

    def start_drag(self, event):
        """开始旋转"""
        self.last_x = event.x
        self.last_y = event.y

    def drag(self, event):
        """旋转处理"""
        dx = event.x - self.last_x
        dy = event.y - self.last_y
        self.angle_y += dx * 0.5
        self.angle_x -= dy * 0.5
        self.last_x = event.x
        self.last_y = event.y
        self.render_molecule()

    def start_pan(self, event):
        """开始平移"""
        self.last_pan_x = event.x
        self.last_pan_y = event.y

    def pan(self, event):
        """平移处理"""
        dx = event.x - self.last_pan_x
        dy = event.y - self.last_pan_y
        self.cam_offset_x += dx
        self.cam_offset_y += dy
        self.last_pan_x = event.x
        self.last_pan_y = event.y
        self.render_molecule()

    def on_zoom_wheel(self, event):
        """缩放处理"""
        if event.delta > 0:
            self.scale = min(1000, self.scale + 10)
        else:
            self.scale = max(10, self.scale - 10)
        self.render_molecule()

    def on_perspective_wheel(self, event):
        """透视调整"""
        if event.delta > 0:
            self.perspective_d = min(20, self.perspective_d + 0.5)
        else:
            self.perspective_d = max(1, self.perspective_d - 0.5)
        self.render_molecule()

    def reset_view(self):
        """重置视图"""
        self.angle_x = self.angle_y = self.angle_z = 0
        self.scale = 100
        self.perspective_d = 5
        self.cam_offset_x = self.cam_offset_y = 0
        self.render_molecule()

    def zoom_in(self):
        """放大"""
        self.scale = min(500, self.scale + 10)
        self.render_molecule()

    def zoom_out(self):
        """缩小"""
        self.scale = max(50, self.scale - 10)
        self.render_molecule()

    def adjust_perspective(self, delta):
        """调整透视"""
        self.perspective_d = max(1, min(20, self.perspective_d + delta))
        self.render_molecule()

    def render_molecule(self):
        """渲染分子"""
        self.canvas.delete("all")
        width, height = 600, 600
        
        # 获取旋转矩阵
        rotation_matrix = self.get_rotation_matrix()
        rotated_coords = np.dot(self.atoms, rotation_matrix.T)
        
        # 根据模型类型调整缩放
        if self.model_type == "space_filling":
            scale_factor = 0.2
            rotated_coords *= scale_factor
            
        # 绘制分子键
        if self.model_type == "ball_and_stick":
            self._render_bonds(rotated_coords, width, height)
            
        # 绘制原子
        self._render_atoms(rotated_coords, width, height)
        
        # 绘制方向立方体
        self.draw_direction_cube(rotation_matrix)

    def _render_bonds(self, rotated_coords, width, height):
        """渲染分子键"""
        sorted_bonds = sorted(self.bonds, key=lambda b: (rotated_coords[b[0], 2] + rotated_coords[b[1], 2]) / 2)
        for begin_idx, end_idx, bond_type in sorted_bonds:
            z1 = rotated_coords[begin_idx, 2]
            z2 = rotated_coords[end_idx, 2]
            
            # 透视投影计算
            x1 = rotated_coords[begin_idx, 0] * self.perspective_d / (self.perspective_d - z1)
            y1 = rotated_coords[begin_idx, 1] * self.perspective_d / (self.perspective_d - z1)
            x2 = rotated_coords[end_idx, 0] * self.perspective_d / (self.perspective_d - z2)
            y2 = rotated_coords[end_idx, 1] * self.perspective_d / (self.perspective_d - z2)
            
            # 转换到画布坐标
            cx1 = x1 * self.scale + width/2 + self.cam_offset_x
            cy1 = -y1 * self.scale + height/2 + self.cam_offset_y
            cx2 = x2 * self.scale + width/2 + self.cam_offset_x
            cy2 = -y2 * self.scale + height/2 + self.cam_offset_y
            
            # 根据键类型绘制
            if bond_type == 1:  # 单键
                self.canvas.create_line(cx1, cy1, cx2, cy2, fill="#404040", width=2)
            elif bond_type == 2:  # 双键
                offset = 3
                self.canvas.create_line(cx1-offset, cy1-offset, cx2-offset, cy2-offset, fill="#404040", width=2)
                self.canvas.create_line(cx1+offset, cy1+offset, cx2+offset, cy2+offset, fill="#404040", width=2)
            elif bond_type == 3:  # 三键
                offset = 4
                self.canvas.create_line(cx1, cy1, cx2, cy2, fill="#404040", width=2)
                self.canvas.create_line(cx1-offset, cy1-offset, cx2-offset, cy2-offset, fill="#404040", width=2)
                self.canvas.create_line(cx1+offset, cy1+offset, cx2+offset, cy2+offset, fill="#404040", width=2)

    def _render_atoms(self, rotated_coords, width, height):
        """渲染原子"""
        sorted_indices = np.argsort(rotated_coords[:, 2])
        for i in sorted_indices:
            x, y, z = rotated_coords[i]
            
            # 透视投影
            x_proj = x * self.perspective_d / (self.perspective_d - z)
            y_proj = y * self.perspective_d / (self.perspective_d - z)
            
            # 转换到画布坐标
            cx = x_proj * self.scale + width/2 + self.cam_offset_x
            cy = -y_proj * self.scale + height/2 + self.cam_offset_y
            
            # 根据模型类型设置原子大小
            if self.model_type == "ball_and_stick":
                atom_size = max(3, 8 * (self.perspective_d / (self.perspective_d - z)))
            else:  # space_filling
                symbol = self.atom_colors[i][1]
                atom_size = self.get_vdw_radius(symbol) * self.scale / 10
            
            # 绘制原子
            atom_color = self.atom_colors[i][0]
            self.canvas.create_oval(
                cx - atom_size, cy - atom_size,
                cx + atom_size, cy + atom_size,
                fill=atom_color, outline=""
            )

    def draw_direction_cube(self, rotation_matrix):
        """绘制方向立方体"""
        size = 50
        margin = 10
        center_x = 600 - size - margin
        center_y = margin + size

        # 定义立方体顶点
        cube_vertices = np.array([
            [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],  # 后面
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]       # 前面
        ]) * (size / 2)

        # 应用旋转
        rotated_vertices = np.dot(cube_vertices, rotation_matrix.T)

        # 投影到2D
        projected_vertices = []
        for vertex in rotated_vertices:
            x_proj = vertex[0] + center_x
            y_proj = -vertex[1] + center_y
            projected_vertices.append((x_proj, y_proj))

        # 定义立方体边
        cube_edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),  # 后面
            (4, 5), (5, 6), (6, 7), (7, 4),  # 前面
            (0, 4), (1, 5), (2, 6), (3, 7)   # 连接边
        ]

        # 绘制边
        for start, end in cube_edges:
            x1, y1 = projected_vertices[start]
            x2, y2 = projected_vertices[end]
            self.canvas.create_line(x1, y1, x2, y2, fill="#FFFFFF")

    def get_rotation_matrix(self):
        """获取旋转矩阵"""
        rx = np.radians(self.angle_x)
        ry = np.radians(self.angle_y)
        rz = np.radians(self.angle_z)

        # 轴旋转矩阵
        Rx = np.array([
            [1, 0, 0],
            [0, np.cos(rx), -np.sin(rx)],
            [0, np.sin(rx), np.cos(rx)]
        ])
        
        Ry = np.array([
            [np.cos(ry), 0, np.sin(ry)],
            [0, 1, 0],
            [-np.sin(ry), 0, np.cos(ry)]
        ])
        
        Rz = np.array([
            [np.cos(rz), -np.sin(rz), 0],
            [np.sin(rz), np.cos(rz), 0],
            [0, 0, 1]
        ])

        return Rz @ Rx @ Ry

def analyze_molecule(smiles, lang_dict=None):
    """
    分析分子性质
    :param smiles: str SMILES字符串
    :param lang_dict: dict 语言字典
    :return: dict 分子性质字典
    """
    if lang_dict is None:
        lang_dict = {}
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(lang_dict.get('error_invalid_smiles', '无效的SMILES字符串'))

        # 基本属性
        properties = {
            lang_dict.get("num_atoms", "原子数"): mol.GetNumAtoms(),
            lang_dict.get("num_bonds", "键数"): mol.GetNumBonds(),
            lang_dict.get("num_rings", "环数"): len(mol.GetRingInfo().AtomRings())
        }

        # 统计原子类型
        atom_types = {}
        for atom in mol.GetAtoms():
            symbol = atom.GetSymbol()
            atom_types[symbol] = atom_types.get(symbol, 0) + 1
        
        properties[lang_dict.get("atom_types", "原子类型统计")] = ", ".join(
            f"{symbol}: {count}" for symbol, count in atom_types.items()
        )

        # 统计键类型
        bond_types = {1: 0, 2: 0, 3: 0, 1.5: 0}  # 单键、双键、三键、芳香键
        for bond in mol.GetBonds():
            bond_type = bond.GetBondTypeAsDouble()
            bond_types[bond_type] = bond_types.get(bond_type, 0) + 1
        
        properties[lang_dict.get("bond_types", "键类型统计")] = ", ".join(
            f"{name}: {count}" for type_val, count in bond_types.items() 
            for name in [["单键", "双键", "三键", "芳香键"][int(type_val*2-2)] if type_val in [1, 2, 3] else "芳香键"]
            if count > 0
        )

        return properties
    except Exception as e:
        raise RuntimeError(f"{lang_dict.get('error_occurred', '发生错误')}: {e}")

def load_config():
    """
    从XML文件加载配置
    返回: dict 配置字典
    """
    if not os.path.exists(CONFIG_FILE):
        # 默认配置
        return {
            "language": "en_us",
            "default_database": DATABASE_PATH,
            "available_languages": ["en_us", "zh_cn", "fr_fr", "de_de"]
        }
    
    try:
        tree = ET.parse(CONFIG_FILE)
        root = tree.getroot()
        config = {}
        
        # 读取语言设置
        language_elem = root.find("language")
        if language_elem is not None:
            config["language"] = language_elem.text
            
        # 读取数据库设置
        database_elem = root.find("default_database")
        if database_elem is not None:
            config["default_database"] = database_elem.text
        else:
            config["default_database"] = DATABASE_PATH
            
        # 读取可用语言列表
        available_languages = []
        langs_elem = root.find("available_languages")
        if langs_elem is not None:
            for lang_elem in langs_elem.findall("language"):
                if lang_elem.text:
                    available_languages.append(lang_elem.text)
        
        if not available_languages:
            available_languages = ["en_us", "zh_cn", "fr_fr", "de_de"]
        config["available_languages"] = available_languages
        
        return config
    except ET.ParseError:
        return {
            "language": "en_us",
            "default_database": DATABASE_PATH,
            "available_languages": ["en_us", "zh_cn", "fr_fr", "de_de"]
        }

def save_config(config):
    """
    将配置保存到XML文件
    :param config: dict 配置字典
    """
    root = ET.Element("config")
    
    # 保存语言设置
    if "language" in config:
        language_elem = ET.SubElement(root, "language")
        language_elem.text = str(config["language"])
    
    # 保存数据库设置
    if "default_database" in config:
        database_elem = ET.SubElement(root, "default_database")
        database_elem.text = str(config["default_database"])
    
    # 保存可用语言列表
    if "available_languages" in config:
        langs_elem = ET.SubElement(root, "available_languages")
        for lang in config["available_languages"]:
            lang_elem = ET.SubElement(langs_elem, "language")
            lang_elem.text = str(lang)
    
    # 保存到文件
    tree = ET.ElementTree(root)
    tree.write(CONFIG_FILE, encoding="utf-8", xml_declaration=True)

def load_language(lang_file):
    """
    从XML文件加载语言包
    :param lang_file: str 语言文件路径
    返回: dict 语言字典
    """
    if not os.path.exists(lang_file):
        logging.error(f"找不到语言文件: {lang_file}")
        return {}
    
    try:
        tree = ET.parse(lang_file)
        root = tree.getroot()
        lang_dict = {}
        
        # 遍历所有元素，将标签名作为键，文本内容作为值
        for elem in root.iter():
            if elem.tag != "language" and elem.text:  # 跳过根元素
                lang_dict[elem.tag] = elem.text.strip()
        
        return lang_dict
    except ET.ParseError as e:
        logging.error(f"解析语言文件时出错: {e}")
        return {}

def show_3d_viewer(smiles, model_type="ball_and_stick", lang_dict=None):
    """
    显示3D分子查看器
    :param smiles: str SMILES字符串
    :param model_type: str 模型类型 ("ball_and_stick" 或 "space_filling")
    :param lang_dict: dict 语言字典
    """
    root = tk.Tk()
    viewer = MoleculeViewer(root, smiles, model_type, lang_dict)
    root.mainloop()

