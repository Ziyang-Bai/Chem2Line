from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt

# 1. 加载分子
smiles = "CCO"  # 乙醇
mol = Chem.MolFromSmiles(smiles)

# 2. 生成 2D 结构
mol_2d = Chem.Mol(mol)
AllChem.Compute2DCoords(mol_2d)

# 3. 生成 3D 结构并优化
mol_3d = Chem.AddHs(mol)  # 添加氢原子
AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())  # 生成 3D 坐标
AllChem.UFFOptimizeMolecule(mol_3d)  # 使用 UFF 力场优化

# 4. 叠加显示
fig, axes = plt.subplots(1, 2, figsize=(8, 4))

# 4.1 显示 2D 结构
img_2d = Draw.MolToImage(mol_2d)
axes[0].imshow(img_2d)
axes[0].axis("off")
axes[0].set_title("2D 结构")

# 4.2 显示 3D 结构
img_3d = Draw.MolToImage(mol_3d)
axes[1].imshow(img_3d)
axes[1].axis("off")
axes[1].set_title("3D 结构（力场优化）")

plt.show()
