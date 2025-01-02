import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import matplotlib.pyplot as plt
import numpy as np

def load_smiles_database(file="smiles_database.xml"):
    """
    加载化学式和 SMILES 映射的数据库
    :param file: XML 文件路径
    :return: 字典形式的化学式和 SMILES 映射
    """
    tree = ET.parse(file)
    root = tree.getroot()
    
    smiles_dict = {}
    for compound in root.findall("compound"):
        name = compound.find("formula").text
        smiles = compound.find("smiles").text
        smiles_dict[name] = smiles
    
    return smiles_dict

def formula_to_bondline(formula, smiles_dict):
    """
    从化学式生成可能的键线式图像，并显示在弹出的窗口中
    :param formula: 化学式字符串（例如：C2H6O）
    :param smiles_dict: 化学式与 SMILES 的映射字典
    :return: None
    """
    try:
        # 查找化学式对应的 SMILES
        if formula in smiles_dict:
            formula = smiles_dict[formula]
        else:
            print(f"没有找到 {formula} 的 SMILES 表示。")
            return
        
        # 从 SMILES 生成分子
        mol = Chem.MolFromSmiles(formula)
        if mol is None:
            print(f"无法从 SMILES {formula} 推导结构，请提供正确的 SMILES 表示形式")
            return
        
        # 优化分子结构
        AllChem.Compute2DCoords(mol)
        
        # 生成分子图像
        img = Draw.MolToImage(mol)
        
        # 使用 matplotlib 显示图像
        plt.figure(figsize=(6, 6))
        plt.imshow(img)
        plt.axis('off')  # 不显示坐标轴
        plt.show()  # 弹出窗口并保持
    
    except Exception as e:
        print(f"发生错误: {e}")

# 示例调用
if __name__ == "__main__":
    smiles_dict = load_smiles_database("smiles_database.xml")
    chemical_formula = input("请输入化学式（例如 C2H6O）：")
    formula_to_bondline(chemical_formula, smiles_dict)
    print("请查看弹出的窗口。")
