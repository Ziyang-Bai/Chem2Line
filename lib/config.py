"""配置文件管理模块"""

import os
import xml.etree.ElementTree as ET

CONFIG_FILE = "lib/config/config.xml"

DEFAULT_CONFIG = {
    "language": "zh_cn",
    "default_database": "lib/db/default_database.xml",
    "available_languages": ["zh_cn", "en_us", "fr_fr", "de_de"],
    "record_history": True,
    "default_3d_model": "ball_and_stick",
}


def load_config() -> dict:
    """
    从 XML 文件加载配置。
    如果文件不存在或解析失败，返回默认配置。
    """
    if not os.path.exists(CONFIG_FILE):
        return DEFAULT_CONFIG.copy()

    try:
        tree = ET.parse(CONFIG_FILE)
        root = tree.getroot()
        config = {}

        for child in root:
            if child.tag == "available_languages":
                config["available_languages"] = [
                    lang.text for lang in child.findall("language") if lang.text
                ]
            elif child.tag == "record_history":
                config["record_history"] = (child.text or "true").lower() == "true"
            else:
                config[child.tag] = child.text

        # 填充缺失的默认值
        for key, default_value in DEFAULT_CONFIG.items():
            if key not in config:
                config[key] = default_value

        return config
    except (ET.ParseError, OSError):
        return DEFAULT_CONFIG.copy()


def save_config(config: dict) -> None:
    """
    将配置保存到 XML 文件。

    :param config: 配置字典
    :raises IOError: 无法写入文件时抛出
    """
    root = ET.Element("config")

    for key, value in config.items():
        if key == "available_languages":
            langs_elem = ET.SubElement(root, "available_languages")
            for lang in value:
                lang_elem = ET.SubElement(langs_elem, "language")
                lang_elem.text = lang
        elif isinstance(value, bool):
            child = ET.SubElement(root, key)
            child.text = str(value).lower()
        else:
            child = ET.SubElement(root, key)
            child.text = str(value)

    tree = ET.ElementTree(root)
    os.makedirs(os.path.dirname(CONFIG_FILE), exist_ok=True)
    tree.write(CONFIG_FILE, encoding="utf-8", xml_declaration=True)
