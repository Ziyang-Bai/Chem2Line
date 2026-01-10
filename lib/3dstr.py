import tkinter as tk
from tkinter import Canvas
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem

class MoleculeViewer:
    def __init__(self, root, smiles="C", model_type="ball_and_stick"):
        """
        初始化 MoleculeViewer
        :param root: Tkinter root
        :param smiles: SMILES str
        :param model_type: modle ("ball_and_stick" 或 "space_filling")
        """
        self.root = root
        self.root.title("3D Molecule Viewer")
        self.canvas = Canvas(root, width=600, height=600, bg="black")
        self.canvas.pack()

        # init sight
        self.angle_x = 0
        self.angle_y = 0
        self.angle_z = 0
        self.scale = 100         # zoom 50 500
        self.perspective_d = 5   # pre 0 20

        # cam offsets
        self.cam_offset_x = 0
        self.cam_offset_y = 0

        # mouse rot
        self.last_x = 0
        self.last_y = 0

        # mouse trans
        self.last_pan_x = 0
        self.last_pan_y = 0

        # mdoel
        self.model_type = model_type

        # controls
        self.create_controls()
        
        # load moles
        self.load_molecule(smiles)
        self.setup_mouse_controls()
        self.render_molecule()

    def create_controls(self):
        control_frame = tk.Frame(self.root)
        control_frame.pack(pady=10)

        # controls
        tk.Button(control_frame, text="重置视角", command=self.reset_view).grid(row=0, column=1, padx=5)
        tk.Button(control_frame, text="缩放+", command=self.zoom_in).grid(row=1, column=0, padx=5)
        tk.Button(control_frame, text="缩放-", command=self.zoom_out).grid(row=1, column=2, padx=5)
        tk.Button(control_frame, text="透视+", command=lambda: self.adjust_perspective(1)).grid(row=2, column=0, padx=5)
        tk.Button(control_frame, text="透视-", command=lambda: self.adjust_perspective(-1)).grid(row=2, column=2, padx=5)

    def setup_mouse_controls(self):
        # left rot
        self.canvas.bind("<Button-1>", self.start_drag)
        self.canvas.bind("<B1-Motion>", self.drag)
        # right cam trans
        self.canvas.bind("<Button-3>", self.start_pan)
        self.canvas.bind("<B3-Motion>", self.pan)
        # wheel whithout control
        self.canvas.bind("<MouseWheel>", self.on_zoom_wheel)
        # wheel with contorl
        self.canvas.bind("<Control-MouseWheel>", self.on_perspective_wheel)

    def load_molecule(self, smiles):
        if smiles == "teapot":
            self.show_teapot = True
            self.load_teapot()
            return
        self.show_teapot = False
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("无效的SMILES字符串")
        
        # 3d coord
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        AllChem.MMFFOptimizeMolecule(mol)
        
        # dm
        self.mol = mol
        
        # coord colour
        conf = mol.GetConformer()
        self.atoms = np.array([conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())])
        self.atom_colors = [self.get_atom_color(atom.GetSymbol(), index=i) for i, atom in enumerate(mol.GetAtoms())]
        
        # bonds
        self.bonds = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()) 
                     for bond in mol.GetBonds()]

        # tea
        self.load_teapot()

    def get_atom_color(self, symbol, index=None):
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
        base_color = color_map.get(symbol, "#FF69B4")  # default pink

        # dyn
        if index is not None:
            adjustment = (index % 10) * 0.05  # bright ness
            r, g, b = int(base_color[1:3], 16), int(base_color[3:5], 16), int(base_color[5:7], 16)
            r = min(255, int(r * (1 - adjustment)))
            g = min(255, int(g * (1 - adjustment)))
            b = min(255, int(b * (1 - adjustment)))
            return f"#{r:02X}{g:02X}{b:02X}"

        return base_color

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
        self.perspective_d = max(1, min(20, self.perspective_d + delta))
        self.render_molecule()

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

    def on_zoom_wheel(self, event):
        if event.delta > 0:
            self.scale = min(1000, self.scale + 10)
        else:
            self.scale = max(10, self.scale - 10)
        self.render_molecule()

    def on_perspective_wheel(self, event):
        if event.delta > 0:
            self.perspective_d = min(20, self.perspective_d + 0.5)
        else:
            self.perspective_d = max(1, self.perspective_d - 0.5)
        self.render_molecule()

    def get_rotation_matrix(self):
        rx = np.radians(self.angle_x)
        ry = np.radians(self.angle_y)
        rz = np.radians(self.angle_z)

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

    def get_vdw_radius(self, symbol):
        """获取原子的范德华半径Å"""
        vdw_radii = {
            "H": 1.2,  # 氢
            "C": 1.7,  # 碳
            "N": 1.55, # 氮
            "O": 1.52, # 氧
            "F": 1.47, # 氟
            "Cl": 1.75, # 氯
            "Br": 1.85, # 溴
            "I": 1.98, # 碘
            "S": 1.8,  # 硫
            "P": 1.8   # 磷
        }
        return vdw_radii.get(symbol, 1.5)  # 默认值

    def load_teapot(self):
        body = []
        for t in np.linspace(0, 2*np.pi, 60, endpoint=False):
            x = 1.5 * np.cos(t)
            y = 1.0 * np.sin(t)
            z = 2.4 + 0.18 * np.sin(3*t)
            body.append([x, y, z])
        lid = []
        for t in np.linspace(0, 2*np.pi, 30, endpoint=False):
            x = 0.7 * np.cos(t)
            y = 0.7 * np.sin(t)
            z = 2.8 + 0.09 * np.cos(2*t)
            lid.append([x, y, z])
        spout = []
        for t in np.linspace(0, 1, 20):
            x = 1.5 + 1.2*t + 0.2*np.sin(2*np.pi*t)
            y = 0.25*np.sin(np.pi*t)
            z = 2.5 + 0.3*t + 0.08*np.cos(3*t)
            spout.append([x, y, z])
        handle = []
        for t in np.linspace(-np.pi/2, 3*np.pi/2, 30):
            x = -1.5 + 0.8 * np.cos(t)
            y = 0.95 * np.sin(t)
            z = 2.5 + 0.12 * np.cos(2*t)
            handle.append([x, y, z])
        knob = []
        for t in np.linspace(0, 2*np.pi, 15, endpoint=False):
            x = 0.18 * np.cos(t)
            y = 0.18 * np.sin(t)
            z = 2.95 + 0.05 * np.sin(2*t)
            knob.append([x, y, z])
        self.teapot_lines = [body, lid, spout, handle, knob]

    def render_teapot(self):
        self.canvas.delete("all")
        width, height = 600, 600
        rotation_matrix = self.get_rotation_matrix()
        scale = self.scale * 1.2
        cx, cy = width/2 + self.cam_offset_x, height/2 + self.cam_offset_y
        for line in self.teapot_lines:
            pts = np.dot(np.array(line), rotation_matrix.T)
            proj = []
            for v in pts:
                z = v[2]
                x_proj = v[0] * self.perspective_d / (self.perspective_d - z)
                y_proj = v[1] * self.perspective_d / (self.perspective_d - z)
                px = x_proj * scale + cx
                py = -y_proj * scale + cy
                proj.append((px, py))
            for i in range(len(proj)-1):
                self.canvas.create_line(proj[i][0], proj[i][1], proj[i+1][0], proj[i+1][1], fill="#FFD700", width=3)
            if line is self.teapot_lines[0] or line is self.teapot_lines[1] or line is self.teapot_lines[4]:
                self.canvas.create_line(proj[-1][0], proj[-1][1], proj[0][0], proj[0][1], fill="#FFD700", width=3)
        self.draw_direction_cube(rotation_matrix)

    def render_molecule(self):
        if hasattr(self, 'show_teapot') and self.show_teapot:
            self.render_teapot()
            return
        self.canvas.delete("all")
        width, height = 600, 600
        
        rotation_matrix = self.get_rotation_matrix()
        rotated_coords = np.dot(self.atoms, rotation_matrix.T)
        
        if self.model_type == "space_filling":
            scale_factor = 0.2
            rotated_coords *= scale_factor


        if self.model_type == "ball_and_stick":
            sorted_bonds = sorted(
                self.bonds,
                key=lambda bond: (rotated_coords[bond[0], 2] + rotated_coords[bond[1], 2]) / 2
            )
            for bond in sorted_bonds:
                i, j = bond
                z1 = rotated_coords[i, 2]
                z2 = rotated_coords[j, 2]
                
                x1 = rotated_coords[i, 0] * self.perspective_d / (self.perspective_d - z1)
                y1 = rotated_coords[i, 1] * self.perspective_d / (self.perspective_d - z1)
                x2 = rotated_coords[j, 0] * self.perspective_d / (self.perspective_d - z2)
                y2 = rotated_coords[j, 1] * self.perspective_d / (self.perspective_d - z2)
                
                cx1 = x1 * self.scale + width/2 + self.cam_offset_x
                cy1 = -y1 * self.scale + height/2 + self.cam_offset_y
                cx2 = x2 * self.scale + width/2 + self.cam_offset_x
                cy2 = -y2 * self.scale + height/2 + self.cam_offset_y
                
                bond_type = self.get_bond_type(i, j)
                if bond_type == 1:
                    self.canvas.create_line(cx1, cy1, cx2, cy2, fill="#404040", width=2)
                elif bond_type == 2:
                    offset = 3
                    self.canvas.create_line(cx1 - offset, cy1 - offset, cx2 - offset, cy2 - offset, fill="#404040", width=2)
                    self.canvas.create_line(cx1 + offset, cy1 + offset, cx2 + offset, cy2 + offset, fill="#404040", width=2)
                elif bond_type == 3:
                    offset = 4
                    self.canvas.create_line(cx1, cy1, cx2, cy2, fill="#404040", width=2)
                    self.canvas.create_line(cx1 - offset, cy1 - offset, cx2 - offset, cy2 - offset, fill="#404040", width=2)
                    self.canvas.create_line(cx1 + offset, cy1 + offset, cx2 + offset, cy2 + offset, fill="#404040", width=2)

        sorted_indices = np.argsort(rotated_coords[:, 2])
        for i in sorted_indices:
            x, y, z = rotated_coords[i]
            
            x_proj = x * self.perspective_d / (self.perspective_d - z)
            y_proj = y * self.perspective_d / (self.perspective_d - z)
            
            cx = x_proj * self.scale + width/2 + self.cam_offset_x
            cy = -y_proj * self.scale + height/2 + self.cam_offset_y
            
            if self.model_type == "ball_and_stick":
                atom_size = max(3, 8 * (self.perspective_d / (self.perspective_d - z)))
            elif self.model_type == "space_filling":
                atom_size = self.get_vdw_radius(self.mol.GetAtomWithIdx(int(i)).GetSymbol()) * self.scale / 10
            
            atom_color = self.get_atom_color(self.mol.GetAtomWithIdx(int(i)).GetSymbol(), index=i)
            
            self.canvas.create_oval(
                cx - atom_size, cy - atom_size,
                cx + atom_size, cy + atom_size,
                fill=atom_color, outline=""
            )

        self.draw_direction_cube(rotation_matrix)

    def draw_direction_cube(self, rotation_matrix):
        """立方体"""
        size = 50 
        margin = 10  # 间距
        center_x = 600 - size - margin
        center_y = margin + size

        # 局部坐标系
        cube_vertices = np.array([
            [-1, -1, -1], [1, -1, -1], [1, 1, -1], [-1, 1, -1],  # 后面
            [-1, -1, 1], [1, -1, 1], [1, 1, 1], [-1, 1, 1]       # 前面
        ]) * (size / 2)

        # 旋转矩阵
        rotated_vertices = np.dot(cube_vertices, rotation_matrix.T)

        # 投影
        projected_vertices = []
        for vertex in rotated_vertices:
            x_proj = vertex[0] + center_x
            y_proj = -vertex[1] + center_y
            projected_vertices.append((x_proj, y_proj))

        # 立方体边
        cube_edges = [
            (0, 1), (1, 2), (2, 3), (3, 0),  # 后面
            (4, 5), (5, 6), (6, 7), (7, 4),  # 前面
            (0, 4), (1, 5), (2, 6), (3, 7)   # 连接
        ]

        # 绘制
        for edge in cube_edges:
            start, end = edge
            x1, y1 = projected_vertices[start]
            x2, y2 = projected_vertices[end]
            self.canvas.create_line(x1, y1, x2, y2, fill="#FFFFFF")

        

    def get_bond_type(self, atom1_idx, atom2_idx):
        """获取键的类型"""
        for bond in self.mol.GetBonds():
            if (bond.GetBeginAtomIdx() == atom1_idx and bond.GetEndAtomIdx() == atom2_idx) or \
               (bond.GetBeginAtomIdx() == atom2_idx and bond.GetEndAtomIdx() == atom1_idx):
                return bond.GetBondTypeAsDouble()
        return 1


if __name__ == "__main__":
    root = tk.Tk()
    # 你可以在这里切换为teapot或SMILES
    #viewer = MoleculeViewer(root, "teapot", "space_filling")
    viewer = MoleculeViewer(root, "C1=CC=C(C=C1)C=O", "ball_and_stick")
    root.mainloop()
