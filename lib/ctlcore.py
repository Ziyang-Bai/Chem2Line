# ctlcore.py
import os
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
CORE_VERSION = "1.1"
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
    return matches

def formula_to_bondline(smiles):
    """
    从 SMILES 生成可能的键线式图像
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

        # 生成分子图像
        img = Draw.MolToImage(mol)
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

if __name__ == "__main__":
    print("核心版本：",core_version())
    print("本模块为Chem2Line的核心模块，无法直接运行，请使用ctlgui.")