"""化学计算模块：键线式生成、分子分析、力场优化"""

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors
from PIL import Image
from typing import Optional


def formula_to_bondline(smiles: str, size: tuple = (300, 300)) -> Image.Image:
    """
    根据 SMILES 生成键线式 PNG 图像。

    :param smiles: SMILES 字符串
    :param size: 图像尺寸 (宽, 高)
    :return: PIL Image 对象
    :raises ValueError: SMILES 无效时抛出
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无效的 SMILES 字符串: {smiles}")
    return Draw.MolToImage(mol, size=size)


def formula_to_svg(smiles: str, file_path: str) -> None:
    """
    根据 SMILES 生成键线式 SVG 文件。

    :param smiles: SMILES 字符串
    :param file_path: 输出 SVG 文件路径
    :raises ValueError: SMILES 无效时抛出
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无效的 SMILES 字符串: {smiles}")
    Draw.MolToFile(mol, file_path, imageType="svg")


def overlay_force_field(smiles: str) -> Image.Image:
    """
    对分子进行力场优化后生成图像。

    :param smiles: SMILES 字符串
    :return: PIL Image 对象
    :raises ValueError: SMILES 无效时抛出
    :raises RuntimeError: 力场优化失败时抛出
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无效的 SMILES 字符串: {smiles}")

    mol = Chem.AddHs(mol)
    try:
        AllChem.Compute2DCoords(mol)
        AllChem.MMFFOptimizeMolecule(mol)
    except Exception as e:
        raise RuntimeError(f"力场优化失败: {e}")

    return Draw.MolToImage(mol, kekulize=True)


def analyze_molecule(smiles: str, lang_dict: Optional[dict] = None) -> dict:
    """
    分析分子的基本性质。

    :param smiles: SMILES 字符串
    :param lang_dict: 语言字典（用于本地化属性名称）
    :return: 分子性质字典
    :raises ValueError: SMILES 无效时抛出
    """
    if lang_dict is None:
        lang_dict = {}

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(
            lang_dict.get("error_invalid_smiles", "无效的 SMILES 字符串")
        )

    return {
        lang_dict.get("molecular_weight", "分子量"): round(Descriptors.MolWt(mol), 4),
        lang_dict.get("num_atoms", "原子数"): mol.GetNumAtoms(),
        lang_dict.get("num_bonds", "键数"): mol.GetNumBonds(),
        lang_dict.get("num_rings", "环数"): len(Chem.GetSSSR(mol)),
        lang_dict.get("molecular_formula", "分子式"): rdMolDescriptors.CalcMolFormula(mol),
    }


def get_chemical_info(smiles: str) -> dict:
    """
    获取分子的详细化学信息（药物化学相关描述符）。

    :param smiles: SMILES 字符串
    :return: 化学信息字典
    :raises ValueError: SMILES 无效时抛出
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"无效的 SMILES 字符串: {smiles}")

    return {
        "分子量": round(Descriptors.MolWt(mol), 4),
        "LogP": round(Descriptors.MolLogP(mol), 4),
        "氢键受体数量": Descriptors.NumHAcceptors(mol),
        "氢键供体数量": Descriptors.NumHDonors(mol),
        "旋转键数量": Descriptors.NumRotatableBonds(mol),
        "拓扑极性表面积": round(Descriptors.TPSA(mol), 4),
    }
