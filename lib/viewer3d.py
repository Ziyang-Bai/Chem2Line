"""3D 分子查看器模块

提供基于 Tkinter Canvas 的 3D 分子可视化功能，
支持球棍模型和比例模型两种显示方式。
"""

import tkinter as tk
from tkinter import Canvas
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

# 原子颜色映射（CPK 配色方案）
ATOM_COLORS = {
    "H": "#FFFFFF",
    "C": "#00FFFF",
    "N": "#0000FF",
    "O": "#FF0000",
    "F": "#00FF00",
    "Cl": "#00FF00",
    "Br": "#8B0000",
    "I": "#9400D3",
    "S": "#FFFF00",
    "P": "#FFA500",
}
DEFAULT_ATOM_COLOR = "#FF69B4"

# 范德华半径（单位：Å）
VDW_RADII = {
    "H": 1.2,
    "C": 1.7,
    "N": 1.55,
    "O": 1.52,
    "F": 1.47,
    "Cl": 1.75,
    "Br": 1.85,
    "I": 1.98,
    "S": 1.8,
    "P": 1.8,
}
DEFAULT_VDW_RADIUS = 1.5


def _get_atom_color(symbol: str) -> str:
    """获取原子的 CPK 颜色"""
    return ATOM_COLORS.get(symbol, DEFAULT_ATOM_COLOR)


def _get_vdw_radius(symbol: str) -> float:
    """获取原子的范德华半径"""
    return VDW_RADII.get(symbol, DEFAULT_VDW_RADIUS)


def _rotation_matrix(angle_x: float, angle_y: float, angle_z: float) -> np.ndarray:
    """
    计算组合旋转矩阵 (Y → X → Z 顺序)。

    :param angle_x: 绕 X 轴旋转角度（度）
    :param angle_y: 绕 Y 轴旋转角度（度）
    :param angle_z: 绕 Z 轴旋转角度（度）
    :return: 3x3 旋转矩阵
    """
    rx = np.radians(angle_x)
    ry = np.radians(angle_y)
    rz = np.radians(angle_z)

    Rx = np.array([
        [1, 0, 0],
        [0, np.cos(rx), -np.sin(rx)],
        [0, np.sin(rx), np.cos(rx)],
    ])
    Ry = np.array([
        [np.cos(ry), 0, np.sin(ry)],
        [0, 1, 0],
        [-np.sin(ry), 0, np.cos(ry)],
    ])
    Rz = np.array([
        [np.cos(rz), -np.sin(rz), 0],
        [np.sin(rz), np.cos(rz), 0],
        [0, 0, 1],
    ])

    return Rz @ Rx @ Ry


class MoleculeViewer:
    """
    基于 Tkinter Canvas 的 3D 分子查看器。

    支持：
    - 左键拖动旋转
    - 右键拖动平移
    - 滚轮缩放
    - Ctrl+滚轮调整透视
    - 球棍模型 / 比例模型切换
    """

    CANVAS_SIZE = 600

    def __init__(self, root: tk.Toplevel, smiles: str = "C",
                 model_type: str = "ball_and_stick", lang_dict: dict = None):
        """
        :param root: Tkinter 父窗口
        :param smiles: 要显示的分子 SMILES
        :param model_type: 模型类型 ("ball_and_stick" 或 "space_filling")
        :param lang_dict: 语言字典
        """
        self.root = root
        self.lang_dict = lang_dict or {}
        self.model_type = model_type

        self.root.title(self.lang_dict.get("molecule_viewer_3d", "3D Molecule Viewer"))
        self.canvas = Canvas(
            root, width=self.CANVAS_SIZE, height=self.CANVAS_SIZE, bg="black"
        )
        self.canvas.pack()

        # 视角参数
        self.angle_x = 0.0
        self.angle_y = 0.0
        self.angle_z = 0.0
        self.scale = 100.0
        self.perspective_d = 5.0

        # 摄像机平移
        self.cam_offset_x = 0.0
        self.cam_offset_y = 0.0

        # 鼠标状态
        self._last_x = 0
        self._last_y = 0
        self._last_pan_x = 0
        self._last_pan_y = 0

        # 分子数据
        self.atoms = None
        self.atom_info = []  # [(color, symbol), ...]
        self.bonds = []  # [(begin_idx, end_idx, bond_order), ...]

        self._create_controls()
        self._bind_mouse()
        self._load_molecule(smiles)
        self._render()

    # ─── 控件 ───────────────────────────────────────────────

    def _create_controls(self):
        frame = tk.Frame(self.root)
        frame.pack(pady=10)

        def _t(key, default):
            return self.lang_dict.get(key) or default

        tk.Button(frame, text=_t("reset_view", "重置视角"),
                  command=self._reset_view).grid(row=0, column=1, padx=5)
        tk.Button(frame, text=_t("zoom_in", "缩放+"),
                  command=self._zoom_in).grid(row=1, column=0, padx=5)
        tk.Button(frame, text=_t("zoom_out", "缩放-"),
                  command=self._zoom_out).grid(row=1, column=2, padx=5)
        tk.Button(frame, text=_t("perspective_increase", "透视+"),
                  command=lambda: self._adjust_perspective(1)).grid(row=2, column=0, padx=5)
        tk.Button(frame, text=_t("perspective_decrease", "透视-"),
                  command=lambda: self._adjust_perspective(-1)).grid(row=2, column=2, padx=5)

    def _bind_mouse(self):
        self.canvas.bind("<Button-1>", self._on_drag_start)
        self.canvas.bind("<B1-Motion>", self._on_drag)
        self.canvas.bind("<Button-3>", self._on_pan_start)
        self.canvas.bind("<B3-Motion>", self._on_pan)
        self.canvas.bind("<MouseWheel>", self._on_zoom_wheel)
        self.canvas.bind("<Control-MouseWheel>", self._on_perspective_wheel)

    # ─── 分子加载 ─────────────────────────────────────────────

    def _load_molecule(self, smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"无效的 SMILES 字符串: {smiles}")

        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)

        conf = mol.GetConformer()
        self.atoms = np.array([
            conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())
        ])
        self.atom_info = [
            (_get_atom_color(atom.GetSymbol()), atom.GetSymbol())
            for atom in mol.GetAtoms()
        ]
        self.bonds = [
            (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx(), bond.GetBondTypeAsDouble())
            for bond in mol.GetBonds()
        ]

    # ─── 视角控制 ─────────────────────────────────────────────

    def _reset_view(self):
        self.angle_x = self.angle_y = self.angle_z = 0.0
        self.scale = 100.0
        self.perspective_d = 5.0
        self.cam_offset_x = self.cam_offset_y = 0.0
        self._render()

    def _zoom_in(self):
        self.scale = min(500, self.scale + 10)
        self._render()

    def _zoom_out(self):
        self.scale = max(50, self.scale - 10)
        self._render()

    def _adjust_perspective(self, delta: float):
        self.perspective_d = max(1, min(20, self.perspective_d + delta))
        self._render()

    # ─── 鼠标事件 ─────────────────────────────────────────────

    def _on_drag_start(self, event):
        self._last_x = event.x
        self._last_y = event.y

    def _on_drag(self, event):
        dx = event.x - self._last_x
        dy = event.y - self._last_y
        self.angle_y += dx * 0.5
        self.angle_x -= dy * 0.5
        self._last_x = event.x
        self._last_y = event.y
        self._render()

    def _on_pan_start(self, event):
        self._last_pan_x = event.x
        self._last_pan_y = event.y

    def _on_pan(self, event):
        dx = event.x - self._last_pan_x
        dy = event.y - self._last_pan_y
        self.cam_offset_x += dx
        self.cam_offset_y += dy
        self._last_pan_x = event.x
        self._last_pan_y = event.y
        self._render()

    def _on_zoom_wheel(self, event):
        delta = 10 if event.delta > 0 else -10
        self.scale = max(10, min(1000, self.scale + delta))
        self._render()

    def _on_perspective_wheel(self, event):
        delta = 0.5 if event.delta > 0 else -0.5
        self.perspective_d = max(1, min(20, self.perspective_d + delta))
        self._render()

    # ─── 渲染 ─────────────────────────────────────────────────

    def _project(self, x: float, y: float, z: float):
        """透视投影 + 画布坐标转换"""
        factor = self.perspective_d / (self.perspective_d - z)
        cx = x * factor * self.scale + self.CANVAS_SIZE / 2 + self.cam_offset_x
        cy = -y * factor * self.scale + self.CANVAS_SIZE / 2 + self.cam_offset_y
        return cx, cy, factor

    def _render(self):
        self.canvas.delete("all")
        if self.atoms is None:
            return

        rot = _rotation_matrix(self.angle_x, self.angle_y, self.angle_z)
        coords = self.atoms @ rot.T

        # 比例模型缩放
        if self.model_type == "space_filling":
            coords = coords * 0.2

        # 绘制键（球棍模型）
        if self.model_type == "ball_and_stick":
            self._render_bonds(coords)

        # 绘制原子
        self._render_atoms(coords)

        # 方向立方体
        self._render_direction_cube(rot)

    def _render_bonds(self, coords: np.ndarray):
        sorted_bonds = sorted(
            self.bonds, key=lambda b: (coords[b[0], 2] + coords[b[1], 2]) / 2
        )
        for begin, end, bond_type in sorted_bonds:
            cx1, cy1, _ = self._project(*coords[begin])
            cx2, cy2, _ = self._project(*coords[end])

            if bond_type == 2:
                dx = cx2 - cx1
                dy = cy2 - cy1
                length = max((dx ** 2 + dy ** 2) ** 0.5, 1e-6)
                ox = -dy / length * 4
                oy = dx / length * 4
                self.canvas.create_line(cx1 + ox, cy1 + oy, cx2 + ox, cy2 + oy,
                                        fill="#404040", width=2)
                self.canvas.create_line(cx1 - ox, cy1 - oy, cx2 - ox, cy2 - oy,
                                        fill="#404040", width=2)
            elif bond_type == 3:
                dx = cx2 - cx1
                dy = cy2 - cy1
                length = max((dx ** 2 + dy ** 2) ** 0.5, 1e-6)
                ox = -dy / length * 4
                oy = dx / length * 4
                self.canvas.create_line(cx1, cy1, cx2, cy2, fill="#404040", width=2)
                self.canvas.create_line(cx1 + ox, cy1 + oy, cx2 + ox, cy2 + oy,
                                        fill="#404040", width=2)
                self.canvas.create_line(cx1 - ox, cy1 - oy, cx2 - ox, cy2 - oy,
                                        fill="#404040", width=2)
            else:
                self.canvas.create_line(cx1, cy1, cx2, cy2, fill="#404040", width=2)

    def _render_atoms(self, coords: np.ndarray):
        sorted_indices = np.argsort(coords[:, 2])
        for i in sorted_indices:
            x, y, z = coords[i]
            cx, cy, factor = self._project(x, y, z)

            if self.model_type == "ball_and_stick":
                atom_size = max(3, 8 * factor)
            else:
                symbol = self.atom_info[i][1]
                atom_size = _get_vdw_radius(symbol) * self.scale / 10

            color = self.atom_info[i][0]
            self.canvas.create_oval(
                cx - atom_size, cy - atom_size,
                cx + atom_size, cy + atom_size,
                fill=color, outline="",
            )

    def _render_direction_cube(self, rotation_matrix: np.ndarray):
        """在右上角绘制方向指示立方体"""
        size = 50
        margin = 10
        center_x = self.CANVAS_SIZE - size - margin
        center_y = margin + size

        vertices = np.array([
            [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1],
        ]) * (size / 2)

        rotated = vertices @ rotation_matrix.T
        projected = [(v[0] + center_x, -v[1] + center_y) for v in rotated]

        edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),
            (4, 5), (5, 6), (6, 7), (7, 4),
            (0, 4), (1, 5), (2, 6), (3, 7),
        ]
        for start, end in edges:
            x1, y1 = projected[start]
            x2, y2 = projected[end]
            self.canvas.create_line(x1, y1, x2, y2, fill="#FFFFFF")
