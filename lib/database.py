"""SMILES 数据库加载与查询模块"""

import os
import xml.etree.ElementTree as ET
from typing import Dict, List, Optional


def load_smiles_database(file_path: str) -> Dict[str, List[str]]:
    """
    加载 XML 格式的 SMILES 数据库。

    :param file_path: XML 数据库文件路径
    :return: 字典，key 为化学式，value 为对应的 SMILES 列表
    :raises FileNotFoundError: 文件不存在时抛出
    :raises ET.ParseError: XML 解析失败时抛出
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"找不到数据库文件: {file_path}")

    tree = ET.parse(file_path)
    root = tree.getroot()

    smiles_dict: Dict[str, List[str]] = {}
    for compound in root.findall("compound"):
        formula_elem = compound.find("formula")
        smiles_elem = compound.find("smiles")
        if formula_elem is None or smiles_elem is None:
            continue
        formula = formula_elem.text
        smiles = smiles_elem.text
        if not formula or not smiles:
            continue

        if formula in smiles_dict:
            smiles_dict[formula].append(smiles)
        else:
            smiles_dict[formula] = [smiles]

    return smiles_dict


def get_smiles_options(formula: str, smiles_dict: Dict[str, List[str]]) -> List[str]:
    """
    获取对应化学式的所有可能 SMILES 表示。
    支持直接输入 SMILES 字符串。

    :param formula: 输入的化学式或 SMILES
    :param smiles_dict: 化学式→SMILES 映射字典
    :return: 匹配的 SMILES 列表
    """
    matches = []
    for key, values in smiles_dict.items():
        if key == formula or formula in values:
            matches.extend(values)

    if not matches:
        # 没有在数据库中找到，将输入本身作为 SMILES 返回
        matches.append(formula)

    return matches


def get_database_info(file_path: str) -> dict:
    """
    获取数据库的元信息。

    :param file_path: XML 数据库文件路径
    :return: 包含文件名、大小、简介、发布者、元素量的字典
    """
    if not os.path.exists(file_path):
        return {"文件名": file_path, "大小": "未知", "简介": "无", "发布者": "未知", "元素量": 0}

    tree = ET.parse(file_path)
    root = tree.getroot()

    description = root.find("description")
    publisher = root.find("publisher")
    compounds = root.findall("compound")

    return {
        "文件名": file_path,
        "大小": f"{os.path.getsize(file_path) / 1024:.2f} KB",
        "简介": description.text if description is not None else "无",
        "发布者": publisher.text if publisher is not None else "未知",
        "元素量": len(compounds),
    }
