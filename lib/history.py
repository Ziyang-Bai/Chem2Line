"""历史记录管理模块"""

import os
import time
import xml.etree.ElementTree as ET
from typing import List

HISTORY_FILE = "lib/history.xml"


def load_history() -> List[dict]:
    """
    从 XML 文件加载历史记录。

    :return: 历史记录列表，每项为 {"smiles": ..., "timestamp": ..., "input_text": ...}
    """
    if not os.path.exists(HISTORY_FILE) or os.path.getsize(HISTORY_FILE) == 0:
        return []

    try:
        tree = ET.parse(HISTORY_FILE)
        root = tree.getroot()
        history = []
        for entry in root.findall("entry"):
            smiles_elem = entry.find("smiles")
            timestamp_elem = entry.find("timestamp")
            input_text_elem = entry.find("input_text")

            smiles = smiles_elem.text if smiles_elem is not None else None
            timestamp = timestamp_elem.text if timestamp_elem is not None else None
            input_text = input_text_elem.text if input_text_elem is not None else None

            if all([smiles, timestamp, input_text]):
                history.append({
                    "smiles": smiles,
                    "timestamp": timestamp,
                    "input_text": input_text,
                })
        return history
    except (ET.ParseError, OSError):
        return []


def save_history(history: List[dict]) -> None:
    """
    将历史记录保存到 XML 文件。

    :param history: 历史记录列表
    """
    root = ET.Element("history")
    for entry in history:
        if not all(entry.get(k) for k in ("smiles", "timestamp", "input_text")):
            continue
        entry_elem = ET.SubElement(root, "entry")
        for key in ("smiles", "timestamp", "input_text"):
            child = ET.SubElement(entry_elem, key)
            child.text = str(entry[key])

    tree = ET.ElementTree(root)
    tree.write(HISTORY_FILE, encoding="utf-8", xml_declaration=True)


def add_to_history(smiles: str, input_text: str, history: List[dict]) -> List[dict]:
    """
    添加一条记录到历史记录中（去重）。

    :param smiles: SMILES 字符串
    :param input_text: 用户输入的原始文本
    :param history: 当前历史记录列表
    :return: 更新后的历史记录列表
    """
    if smiles in [entry["smiles"] for entry in history]:
        return history

    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    history.append({
        "smiles": smiles,
        "timestamp": timestamp,
        "input_text": input_text,
    })
    save_history(history)
    return history


def clear_history() -> List[dict]:
    """
    清空历史记录。

    :return: 空列表
    """
    history = []
    save_history(history)
    return history
