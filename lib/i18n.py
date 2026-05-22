"""国际化（多语言）管理模块"""

import os
import xml.etree.ElementTree as ET

LANG_DIR = "lib/lang"


def load_language(lang_code: str) -> dict:
    """
    加载指定语言代码的语言文件。

    :param lang_code: 语言代码，如 'zh_cn', 'en_us'
    :return: 键值对字典，key 为标签名，value 为翻译文本
    """
    lang_file = os.path.join(LANG_DIR, f"{lang_code}.xml")
    if not os.path.exists(lang_file):
        return {}

    try:
        tree = ET.parse(lang_file)
        root = tree.getroot()
        return {child.tag: (child.text or "").strip() for child in root}
    except (ET.ParseError, OSError):
        return {}


def get_text(lang_dict: dict, key: str, default: str = "") -> str:
    """
    从语言字典中获取翻译文本，若不存在则返回默认值。

    :param lang_dict: 语言字典
    :param key: 翻译键
    :param default: 默认值
    :return: 翻译文本
    """
    return lang_dict.get(key, default) or default
