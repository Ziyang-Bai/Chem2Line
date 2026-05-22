"""独立 3D 分子查看器演示脚本

可直接运行此文件查看 3D 分子效果。
使用 lib.viewer3d.MoleculeViewer 作为核心渲染引擎。
"""

import tkinter as tk
from lib.viewer3d import MoleculeViewer


def main():
    root = tk.Tk()
    # 示例：维生素 D3 的 SMILES
    smiles = "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CCC/C2=C\\C=C/3\\C[C@H](CCC3=C)O)C"
    MoleculeViewer(root, smiles, model_type="ball_and_stick")
    root.mainloop()


if __name__ == "__main__":
    main()
