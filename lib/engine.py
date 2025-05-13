import os
import time
import xml.etree.ElementTree as ET

def load_config():
    """
    加载配置文件
    """
    try:
        tree = ET.parse('lib/config/config.xml')
        root = tree.getroot()
        config = {child.tag: child.text for child in root if child.tag != 'available_languages'}
        config['available_languages'] = [lang.text for lang in root.find('available_languages')]
        config['record_history'] = config.get('record_history', 'true').lower() == 'true'
        return config
    except Exception:
        return {'record_history': True}

def save_config(config):
    """
    保存配置文件
    """
    root = ET.Element("config")
    for key, value in config.items():
        if key == 'available_languages':
            available_languages = ET.SubElement(root, 'available_languages')
            for lang in value:
                lang_element = ET.SubElement(available_languages, 'language')
                lang_element.text = lang
        else:
            child = ET.SubElement(root, key)
            child.text = str(value).lower() if isinstance(value, bool) else value
    tree = ET.ElementTree(root)
    tree.write('lib/config/config.xml')

def load_language(lang_file):
    """
    加载语言文件
    """
    try:
        tree = ET.parse(lang_file)
        root = tree.getroot()
        return {child.tag: child.text for child in root}
    except Exception:
        return {}

def load_history():
    """
    加载历史记录
    """
    history_file = "lib/history.xml"
    if not os.path.exists(history_file) or os.path.getsize(history_file) == 0:
        return []
    try:
        tree = ET.parse(history_file)
        root = tree.getroot()
        return [{"smiles": entry.find("smiles").text, "timestamp": entry.find("timestamp").text, "input_text": entry.find("input_text").text} for entry in root.findall("entry")]
    except ET.ParseError:
        return []

def save_history(history):
    """
    保存历史记录
    """
    root = ET.Element("history")
    for entry in history:
        entry_element = ET.SubElement(root, "entry")
        ET.SubElement(entry_element, "smiles").text = entry["smiles"]
        ET.SubElement(entry_element, "timestamp").text = entry["timestamp"]
        ET.SubElement(entry_element, "input_text").text = entry["input_text"]
    tree = ET.ElementTree(root)
    tree.write("lib/history.xml")

def add_to_history(smiles, input_text, history, update_menu_callback):
    """
    添加到历史记录
    """
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    if smiles not in [entry["smiles"] for entry in history]:
        history.append({"smiles": smiles, "timestamp": timestamp, "input_text": input_text})
        update_menu_callback(history)
        save_history(history)

def clear_history(history, update_menu_callback):
    """
    清空历史记录
    """
    history.clear()
    save_history(history)
    update_menu_callback(history)

def update_history_menu(history_menu, history):
    """
    更新历史记录菜单
    """
    history_menu.delete(0, "end")
    for entry in history[-5:]:
        history_menu.add_command(label=entry["smiles"], command=lambda s=entry["smiles"]: print(f"Insert {s} into input field"))

def search_database(query, smiles_dict):
    """
    查询数据库
    """
    # 示例逻辑，具体实现根据需求调整
    return []

def advanced_search_gui():
    """
    高级查询 GUI
    """
    # 示例逻辑，具体实现根据需求调整
    pass

def sdf_converter_gui():
    """
    SDF 转换工具 GUI
    """
    # 示例逻辑，具体实现根据需求调整
    pass
