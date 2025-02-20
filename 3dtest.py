import tkinter as tk
from tkinter import Canvas
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

class MoleculeViewer:
    def __init__(self, root, smiles="C"):
        self.root = root
        self.root.title("3D Molecule Viewer")
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
        tk.Button(control_frame, text="重置视角", command=self.reset_view).grid(row=0, column=1, padx=5)
        tk.Button(control_frame, text="缩放+", command=self.zoom_in).grid(row=1, column=0, padx=5)
        tk.Button(control_frame, text="缩放-", command=self.zoom_out).grid(row=1, column=2, padx=5)
        tk.Button(control_frame, text="透视+", command=lambda: self.adjust_perspective(1)).grid(row=2, column=0, padx=5)
        tk.Button(control_frame, text="透视-", command=lambda: self.adjust_perspective(-1)).grid(row=2, column=2, padx=5)

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
        
        # 获取键信息
        self.bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) 
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
            i, j = bond
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

# 创建并运行程序
if __name__ == "__main__":
    root = tk.Tk()
    viewer = MoleculeViewer(root, "CCC1=C2C(C(=C(C2=NC3=C1C=C(C=C3C=C(C)C)C=O)C)C)CCC(=O)OC4C(C(C(C(O4)C)O)O)O")
    root.mainloop()
