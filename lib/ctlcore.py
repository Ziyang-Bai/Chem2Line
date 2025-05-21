# ctlcore.py
import os
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import tkinter as tk
from tkinter import Canvas
import numpy as np

CORE_VERSION = "1.2.3"

def core_version():
    return CORE_VERSION

def load_smiles_database(file="default_database.xml"):
    """
    加载化学式和 SMILES 映射的数据库
    :param file: XML 文件路径
    :return: 字典形式的化学式和 SMILES 映射
    """
    if not os.path.exists(file):
        raise FileNotFoundError(f"找不到数据库文件: {file}")

    tree = ET.parse(file)
    root = tree.getroot()

    smiles_dict = {}
    for compound in root.findall("compound"):
        name = compound.find("formula").text
        smiles = compound.find("smiles").text
        if name in smiles_dict:
            smiles_dict[name].append(smiles)
        else:
            smiles_dict[name] = [smiles]

    return smiles_dict  # 这一行应当在 for 循环之后执行

def get_smiles_options(formula, smiles_dict):
    """
    获取对应化学式的所有可能 SMILES
    :param formula: 输入的化学式或 SMILES
    :param smiles_dict: 化学式和 SMILES 的映射字典
    :return: 包含所有可能 SMILES 的列表
    """
    matches = []
    for key, values in smiles_dict.items():
        if key == formula or formula in values:  # 支持直接输入 SMILES
            matches.extend(values)
    if not matches:
        matches.append(formula)  # 如果没有找到匹配的化学式，直接返回输入的 SMILES
    return matches

def formula_to_bondline(smiles, image_type="PNG"):
    """
    根据 SMILES 生成键线式图像。
    :param smiles: SMILES 字符串
    :param image_type: 图像类型 ("PNG" 或 "SVG")
    :return: 图像对象 (PIL.Image 或 SVG 字符串)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("Invalid SMILES string")

    if image_type.upper() == "SVG":
        return Draw.MolsToGridImage([mol], molsPerRow=1, subImgSize=(300, 300), useSVG=True)
    elif image_type.upper() == "PNG":
        return Draw.MolToImage(mol, size=(300, 300))
    else:
        raise ValueError("Unsupported image type. Use 'PNG' or 'SVG'.")

def formula_to_structure(smiles):
    """
    从 SMILES 生成可能的结构式图像
    :param smiles: SMILES 表示
    :return: 图像对象
    """
    try:
        # 从 SMILES 生成分子
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"无法从 SMILES {smiles} 推导结构，请提供正确的 SMILES 表示形式")

        # 优化分子结构
        AllChem.Compute2DCoords(mol)

        # 生成分子图像，设置 kekulize=True 以生成结构式
        img = Draw.MolToImage(mol, kekulize=True)
        return img

    except Exception as e:
        raise RuntimeError(f"发生错误: {e}")

def overlay_force_field(mol):
    """
        叠加力场显示
    :param mol: 分子对象
    :return: 力场叠加图像
    """
    try:
        # 计算力场
        AllChem.Compute2DCoords(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        img = Draw.MolToImage(mol, kekulize=True)
        return img
    except Exception as e:
        raise RuntimeError(f"发生错误: {e}")

def get_database_info(file="default_database.xml"):
    """
    获取数据库的基本信息
    :param file: XML 文件路径
    :return: 字典形式的信息（大小、简介、发布者、元素量）
    """
    if not os.path.exists(file):
        return {"大小": "未知", "简介": "无", "发布者": "未知", "元素量": 0}

    tree = ET.parse(file)
    root = tree.getroot()
    filename = file
    description = root.find("description")
    publisher = root.find("publisher")
    compounds = root.findall("compound")

    info = {
        "文件名": filename,
        "大小": f"{os.path.getsize(file) / 1024:.2f} KB",
        "简介": description.text if description is not None else "无",
        "发布者": publisher.text if publisher is not None else "未知",
        "元素量": len(compounds)
    }
    return info

def analyze_molecule(smiles, lang_dict):
    """
    分析化学分子的性质
    :param smiles: SMILES 表示
    :param lang_dict: 语言字典
    :return: 分子性质的字典
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"{lang_dict.get('error_invalid_smiles', '无法从 SMILES 推导结构，请提供正确的 SMILES 表示形式')}")

        # 计算分子性质
        mol_weight = Descriptors.MolWt(mol)
        num_atoms = mol.GetNumAtoms()
        num_bonds = mol.GetNumBonds()
        num_rings = Chem.GetSSSR(mol)
        formula = Chem.rdMolDescriptors.CalcMolFormula(mol)

        properties = {
            lang_dict.get("molecular_weight", "分子量"): mol_weight,
            lang_dict.get("num_atoms", "原子数"): num_atoms,
            lang_dict.get("num_bonds", "键数"): num_bonds,
            lang_dict.get("num_rings", "环数"): num_rings,
            lang_dict.get("molecular_formula", "分子式"): formula
        }

        return properties

    except Exception as e:
        raise RuntimeError(f"{lang_dict.get('error_occurred', '发生错误')}: {e}")

def get_chemical_info(smiles):
    """
    获取更多化学信息
    :param smiles: SMILES 表示
    :return: 化学信息的字典
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("无效的SMILES字符串")

        info = {
            "分子量": Descriptors.MolWt(mol),
            "LogP": Descriptors.MolLogP(mol),
            "氢键受体数量": Descriptors.NumHAcceptors(mol),
            "氢键供体数量": Descriptors.NumHDonors(mol),
            "旋转键数量": Descriptors.NumRotatableBonds(mol),
            "拓扑极性表面积": Descriptors.TPSA(mol)
        }
        return info

    except Exception as e:
        raise RuntimeError(f"获取化学信息时发生错误: {e}")

def show_chemical_info(smiles):
    """
    显示化学信息
    :param smiles: SMILES 表示
    :return: 化学信息的字典
    """
    try:
        info = get_chemical_info(smiles)
        info_str = "\n".join([f"{key}: {value}" for key, value in info.items()])
        return info_str
    except Exception as e:
        raise RuntimeError(f"获取化学信息时发生错误: {e}")

def search_database(query, smiles_dict):
    """
    在数据库中查询化学物质
    :param query: 查询字符串，可以是名称或化学式
    :param smiles_dict: 已加载的 SMILES 数据库
    :return: 匹配的化学物质列表，每个元素为 (名称, 化学式, SMILES)
    """
    results = []
    for smiles, info in smiles_dict.items():
        if query.lower() in info.get("name", "").lower() or query.lower() in info.get("formula", "").lower():
            results.append((info.get("name", "N/A"), info.get("formula", "N/A"), smiles))
    return results

class MoleculeViewer:
    def __init__(self, root, smiles="C", lang_dict=None):
        self.root = root
        self.lang_dict = lang_dict or {}
        self.root.title(self.lang_dict.get("3d_molecule_viewer", "3D Molecule Viewer"))
        self.canvas = Canvas(root, width=600, height=600, bg="black")
        self.canvas.pack()

        # 初始化视角参数
        self.angle_x = 0
        self.angle_y = 0
        self.angle_z = 0
        self.scale = 100         # 缩放范围从50到500
        self.perspective_d = 5   # 透视参数范围从1到20

        # 摄像机平移偏移量
        self.cam_offset_x = 0
        self.cam_offset_y = 0

        # 鼠标拖动相关变量（旋转）
        self.last_x = 0
        self.last_y = 0

        # 鼠标拖动相关变量（平移）
        self.last_pan_x = 0
        self.last_pan_y = 0

        # 创建控件
        self.create_controls()
        
        # 加载分子
        self.load_molecule(smiles)
        self.setup_mouse_controls()
        self.render_molecule()

    def create_controls(self):
        control_frame = tk.Frame(self.root)
        control_frame.pack(pady=10)

        # 旋转/视图控制按钮
        tk.Button(control_frame, text=self.lang_dict.get("reset_view", "重置视角"), command=self.reset_view).grid(row=0, column=1, padx=5)
        tk.Button(control_frame, text=self.lang_dict.get("zoom_in", "缩放+"), command=self.zoom_in).grid(row=1, column=0, padx=5)
        tk.Button(control_frame, text=self.lang_dict.get("zoom_out", "缩放-"), command=self.zoom_out).grid(row=1, column=2, padx=5)
        tk.Button(control_frame, text=self.lang_dict.get("perspective_increase", "透视+"), command=lambda: self.adjust_perspective(1)).grid(row=2, column=0, padx=5)
        tk.Button(control_frame, text=self.lang_dict.get("perspective_decrease", "透视-"), command=lambda: self.adjust_perspective(-1)).grid(row=2, column=2, padx=5)

    def setup_mouse_controls(self):
        # 左键拖动旋转视角
        self.canvas.bind("<Button-1>", self.start_drag)
        self.canvas.bind("<B1-Motion>", self.drag)
        # 右键拖动平移摄像机位置
        self.canvas.bind("<Button-3>", self.start_pan)
        self.canvas.bind("<B3-Motion>", self.pan)
        # 鼠标滚轮缩放（不按Ctrl时）
        self.canvas.bind("<MouseWheel>", self.on_zoom_wheel)
        # Ctrl+鼠标滚轮调整透视参数
        self.canvas.bind("<Control-MouseWheel>", self.on_perspective_wheel)

    def load_molecule(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("无效的SMILES字符串")
        
        # 生成3D坐标
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # 获取原子坐标和颜色
        conf = mol.GetConformer()
        self.atoms = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
        self.atom_colors = [self.get_atom_color(atom.GetSymbol()) for atom in mol.GetAtoms()]
        
        # 获取键信息，包含类型
        self.bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble())
                     for bond in mol.GetBonds()]

    def get_atom_color(self, symbol):
        color_map = {
            "H": "#FFFFFF",   # 白色
            "C": "#00FFFF",   # 青色
            "N": "#0000FF",   # 蓝色
            "O": "#FF0000",   # 红色
            "F": "#00FF00",   # 绿色
            "Cl": "#00FF00",  # 绿色
            "Br": "#8B0000",  # 深红色
            "I": "#9400D3",   # 紫色
            "S": "#FFFF00",   # 黄色
            "P": "#FFA500"    # 橙色
        }
        return color_map.get(symbol, "#FF69B4")  # 默认粉色

    def reset_view(self):
        self.angle_x = self.angle_y = self.angle_z = 0
        self.scale = 100
        self.perspective_d = 5
        self.cam_offset_x = 0
        self.cam_offset_y = 0
        self.render_molecule()

    def zoom_in(self):
        self.scale = min(500, self.scale + 10)
        self.render_molecule()

    def zoom_out(self):
        self.scale = max(50, self.scale - 10)
        self.render_molecule()

    def adjust_perspective(self, delta):
        # 每次改变1个单位，范围设定在1到20之间
        self.perspective_d = max(1, min(20, self.perspective_d + delta))
        self.render_molecule()

    # 鼠标左键拖动旋转
    def start_drag(self, event):
        self.last_x = event.x
        self.last_y = event.y

    def drag(self, event):
        dx = event.x - self.last_x
        dy = event.y - self.last_y
        self.angle_y += dx * 0.5
        self.angle_x -= dy * 0.5
        self.last_x = event.x
        self.last_y = event.y
        self.render_molecule()

    # 鼠标右键拖动平移摄像机位置
    def start_pan(self, event):
        self.last_pan_x = event.x
        self.last_pan_y = event.y

    def pan(self, event):
        dx = event.x - self.last_pan_x
        dy = event.y - self.last_pan_y
        self.cam_offset_x += dx
        self.cam_offset_y += dy
        self.last_pan_x = event.x
        self.last_pan_y = event.y
        self.render_molecule()

    # 鼠标滚轮缩放（不按Ctrl时）
    def on_zoom_wheel(self, event):
        if event.delta > 0:
            self.scale = min(1000, self.scale + 10)
        else:
            self.scale = max(10, self.scale - 10)
        self.render_molecule()

    # Ctrl+鼠标滚轮调整透视参数
    def on_perspective_wheel(self, event):
        if event.delta > 0:
            self.perspective_d = min(20, self.perspective_d + 0.5)
        else:
            self.perspective_d = max(1, self.perspective_d - 0.5)
        self.render_molecule()

    def get_rotation_matrix(self):
        # 转换为弧度
        rx = np.radians(self.angle_x)
        ry = np.radians(self.angle_y)
        rz = np.radians(self.angle_z)

        # 各轴旋转矩阵
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

        # 组合旋转矩阵：Y -> X -> Z
        return Rz @ Rx @ Ry

    def render_molecule(self):
        self.canvas.delete("all")
        width, height = 600, 600
        
        # 应用旋转矩阵
        rotation_matrix = self.get_rotation_matrix()
        rotated_coords = np.dot(self.atoms, rotation_matrix.T)
        
        # 先对键按照平均深度排序（远处先绘制）
        sorted_bonds = sorted(
            self.bonds,
            key=lambda bond: (rotated_coords[bond[0], 2] + rotated_coords[bond[1], 2]) / 2
        )
        for bond in sorted_bonds:
            i, j, bond_type = bond
            z1 = rotated_coords[i, 2]
            z2 = rotated_coords[j, 2]
            
            # 计算透视投影
            x1 = rotated_coords[i, 0] * self.perspective_d / (self.perspective_d - z1)
            y1 = rotated_coords[i, 1] * self.perspective_d / (self.perspective_d - z1)
            x2 = rotated_coords[j, 0] * self.perspective_d / (self.perspective_d - z2)
            y2 = rotated_coords[j, 1] * self.perspective_d / (self.perspective_d - z2)
            
            # 转换到画布坐标（加入平移偏移量）
            cx1 = x1 * self.scale + width/2 + self.cam_offset_x
            cy1 = -y1 * self.scale + height/2 + self.cam_offset_y
            cx2 = x2 * self.scale + width/2 + self.cam_offset_x
            cy2 = -y2 * self.scale + height/2 + self.cam_offset_y
            
            if bond_type == 2:  # 双键，画两条平行线
                # 计算垂直方向
                dx = cx2 - cx1
                dy = cy2 - cy1
                length = (dx**2 + dy**2) ** 0.5
                if length == 0:
                    offset_x, offset_y = 0, 0
                else:
                    offset_x = -dy / length * 4  # 4像素偏移
                    offset_y = dx / length * 4
                self.canvas.create_line(cx1 + offset_x, cy1 + offset_y, cx2 + offset_x, cy2 + offset_y, fill="#404040", width=2)
                self.canvas.create_line(cx1 - offset_x, cy1 - offset_y, cx2 - offset_x, cy2 - offset_y, fill="#404040", width=2)
            else:
                self.canvas.create_line(cx1, cy1, cx2, cy2, fill="#404040", width=2)
        
        # 对原子按照 z 坐标升序排序（远处先绘制，近处后绘制）
        sorted_indices = np.argsort(rotated_coords[:, 2])
        for i in sorted_indices:
            x, y, z = rotated_coords[i]
            
            # 透视投影
            x_proj = x * self.perspective_d / (self.perspective_d - z)
            y_proj = y * self.perspective_d / (self.perspective_d - z)
            
            # 转换到画布坐标（加入平移偏移量）
            cx = x_proj * self.scale + width/2 + self.cam_offset_x
            cy = -y_proj * self.scale + height/2 + self.cam_offset_y
            
            # 根据原子类型调整大小
            atom_size = max(3, 8 * (self.perspective_d / (self.perspective_d - z)))
            self.canvas.create_oval(
                cx - atom_size, cy - atom_size,
                cx + atom_size, cy + atom_size,
                fill=self.atom_colors[i], outline=""
            )

def show_3d_viewer(smiles):
    root = tk.Tk()
    viewer = MoleculeViewer(root, smiles)
    root.mainloop()

if __name__ == "__main__":
    print("核心版本：",core_version())
    print("本模块为Chem2Line的核心模块，无法直接运行，请使用ctlgui.")