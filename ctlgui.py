#ctlgui.py
import tkinter as tk
from tkinter import filedialog, messagebox, Scrollbar, Canvas, Menu, ttk, Toplevel, StringVar, OptionMenu
from PIL import Image, ImageTk
from lib.ctlcore import *
from lib.engine import *
import time
import threading
import os
import sys
from lib.debug import enable_debug_mode  # 导入调试模块
from rdkit import Chem  # 添加此行以导入 Chem
from rdkit.Chem import Draw

database_path = "lib/db/default_database.xml"
VERSION = "1.4.1"
DEVELOPER = "Ziyang-Bai"
DATE = "2025-01-01"
CORE_VERSION = core_version()

history = []
history_file = "lib/history.xml"

def load_history():
    """
    从 XML 文件加载历史记录
    """
    if not os.path.exists(history_file) or os.path.getsize(history_file) == 0:
        # 如果文件不存在或为空，初始化空历史记录
        return []

    try:
        tree = ET.parse(history_file)
        root = tree.getroot()
        history = []
        for entry in root.findall("entry"):
            smiles = entry.find("smiles").text
            timestamp = entry.find("timestamp").text
            input_text = entry.find("input_text").text
            history.append({"smiles": smiles, "timestamp": timestamp, "input_text": input_text})
        return history
    except ET.ParseError:
        # 如果文件内容无效，返回空历史记录
        return []

def save_history():
    """
    将历史记录保存到 XML 文件
    """
    root = ET.Element("history")
    for entry in history:
        entry_element = ET.SubElement(root, "entry")
        smiles_element = ET.SubElement(entry_element, "smiles")
        smiles_element.text = entry["smiles"]
        timestamp_element = ET.SubElement(entry_element, "timestamp")
        timestamp_element.text = entry["timestamp"]
        input_text_element = ET.SubElement(entry_element, "input_text")
        input_text_element.text = entry["input_text"]
    tree = ET.ElementTree(root)
    tree.write(history_file)

def add_to_history(smiles):
    """
    添加SMILES到历史记录
    """
    if not config.get('record_history', True):
        return  # 如果禁用了历史记录，则不执行任何操作

    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    input_text = formula_entry.get().strip()
    if smiles not in [entry["smiles"] for entry in history]:
        history.append({"smiles": smiles, "timestamp": timestamp, "input_text": input_text})
        update_history_menu(history_menu, history)
        save_history()

def update_history_menu(history_menu, history):
    """
    更新历史记录菜单
    :param history_menu: 历史记录菜单对象
    :param history: 历史记录列表
    """
    history_menu.delete(0, tk.END)
    for entry in history[-5:]:
        history_menu.add_command(label=entry["smiles"], command=lambda s=entry["smiles"]: formula_entry.insert(0, s))

def show_long_history():
    """
    显示长历史记录窗口，使用 Windows 原生的窗口方式，并允许排序和即时删除
    """
    history_window = Toplevel(root)
    history_window.title(lang_dict.get("history_title", "历史记录"))
    history_window.geometry("600x400")
    history_window.iconbitmap("lib/media/nctl.ico")

    columns = ("timestamp", "input_text", "smiles")
    tree = ttk.Treeview(history_window, columns=columns, show="headings")
    tree.heading("timestamp", text=lang_dict.get("timestamp", "时间戳"))
    tree.heading("input_text", text=lang_dict.get("input_text", "输入文本"))
    tree.heading("smiles", text="SMILES")

    for entry in history:
        tree.insert("", "end", values=(entry["timestamp"], entry["input_text"], entry["smiles"]))

    tree.pack(fill="both", expand=True)

    def delete_selected():
        selected_items = tree.selection()
        for item in selected_items:
            item_values = tree.item(item, "values")
            entry_to_delete = next((entry for entry in history if entry["timestamp"] == item_values[0] and entry["input_text"] == item_values[1] and entry["smiles"] == item_values[2]), None)
            if entry_to_delete:
                history.remove(entry_to_delete)
                tree.delete(item)
        save_history()
        update_history_menu(history_menu, history)

    delete_button = tk.Button(history_window, text=(lang_dict.get("delete") or "删除"), command=delete_selected)
    delete_button.pack(pady=10)

    clear_button = tk.Button(history_window, text=(lang_dict.get("clear_history") or "清空历史记录"), command=clear_history)
    clear_button.pack(pady=10)

def delete_history_entry(entry):
    """
    删除历史记录条目
    """
    history.remove(entry)
    save_history()
    update_history_menu(history_menu, history)

def clear_history():
    """
    清空历史记录
    """
    history.clear()
    save_history()
    update_history_menu(history_menu, history)

def show_smiles_selection(smiles_list):
    """
    弹出窗口显示可供选择的 SMILES
    :param smiles_list: 包含多个 SMILES 的列表
    :return: 用户选择的 SMILES
    """
    selection_window = tk.Toplevel()
    selection_window.title(f"Chem2Line - {lang_dict.get('select_smiles_title', '选择 SMILES 表示')}")
    selection_window.geometry("600x400")
    selection_window.iconbitmap("lib/media/nctl.ico")

    selected_smiles = tk.StringVar()

    def select_smiles(smiles):
        selected_smiles.set(smiles)
        selection_window.destroy()

    # 添加滚动条
    canvas = Canvas(selection_window)
    scrollbar = Scrollbar(selection_window, orient="vertical", command=canvas.yview)
    scrollable_frame = tk.Frame(canvas)

    scrollable_frame.bind(
        "<Configure>",
        lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
    )

    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar.set)

    canvas.pack(side="left", fill="both", expand=True)
    scrollbar.pack(side="right", fill="y")

    for smiles in smiles_list:
        try:
            img = formula_to_bondline(smiles)
            img = img.convert("RGBA")  # 确保图像以 RGBA 模式加载
            img = ImageTk.PhotoImage(img)

            frame = tk.Frame(scrollable_frame, pady=10)
            frame.pack(fill="x")

            btn = tk.Button(frame, image=img, text=smiles, compound=tk.TOP,
                             command=lambda s=smiles: select_smiles(s))
            btn.image = img  # 保存引用，防止被垃圾回收
            btn.pack()

            label = tk.Label(frame, text=smiles, font=("Arial", 10))
            label.pack()

        except Exception as e:
            tk.Label(scrollable_frame, text=f"{lang_dict.get('error_loading_smiles', '无法加载')} {smiles}: {e}").pack()

    # 等待用户选择
    selection_window.wait_window()
    return selected_smiles.get()
def on_overlay_force_field():
    input_text = formula_entry.get().strip()
    if not input_text:
        messagebox.showwarning(f"Chem2Line - {lang_dict.get('input_empty_title', '输入为空')}", lang_dict.get("input_empty_message", "请输入化学式或 SMILES"))
        return

    try:
        smiles_list = get_smiles_options(input_text, smiles_dict)

        if not smiles_list:
            raise ValueError(f"{lang_dict.get('error_not_found', '找不到')} {input_text} {lang_dict.get('smiles_representation', '的 SMILES 表示，请检查输入或更换数据库')}")

        if len(smiles_list) == 1:
            selected_smiles = smiles_list[0]
        else:
            selected_smiles = show_smiles_selection(smiles_list)

        if selected_smiles:
            mol = Chem.MolFromSmiles(selected_smiles)
            mol = Chem.AddHs(mol)  # 添加显式氢原子
            img = overlay_force_field(mol)
            img = ImageTk.PhotoImage(img)

            result_label.config(image=img)
            result_label.image = img

    except ValueError as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_not_found_title', '未找到结果')}", f"{lang_dict.get('error_code', '错误代码')}: 1001\n{str(e)}")
    except RuntimeError as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_generation_failed_title', '生成失败')}", f"{lang_dict.get('error_code', '错误代码')}: 1002\n{str(e)}")
    except Exception as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_unknown_title', '未知错误')}", f"{lang_dict.get('error_code', '错误代码')}: 1000\n{lang_dict.get('error_unknown_message', '发生未知错误')}: {e}")

def load_database_with_progress(file_path=None):
    """
    在后台线程中加载数据库，并显示进度条。
    """
    def update_progress():
        for _ in range(100):  # 模拟加载进度
            time.sleep(0.05)
            progress_bar.step(1)
            progress_window.update_idletasks()

    def load_database():
        try:
            nonlocal file_path
            if not file_path:
                file_path = filedialog.askopenfilename(filetypes=[("XML files", "*.xml")])
            if file_path:
                progress_window.deiconify()  # 显示进度条窗口
                progress_bar.start()  # 启动进度条
                threading.Thread(target=update_progress).start()  # 启动进度条的后台更新
                new_smiles_dict = load_smiles_database(file_path)
                global smiles_dict, database_path
                smiles_dict = new_smiles_dict  # 更新全局数据库
                database_path = file_path  # 更新数据库路径
                database_info = get_database_info(file_path)
                info_str = "\n".join([f"{lang_dict.get(key, key)}: {value}" for key, value in database_info.items()])
                messagebox.showinfo(f"Chem2Line - {lang_dict.get('database_changed_title', '数据库已更换')}", f"{lang_dict.get('current_database_info', '当前数据库信息')}:\n{info_str}")
                progress_bar.stop()  # 停止进度条
                progress_window.withdraw()  # 隐藏进度条窗口
        except Exception as e:
            messagebox.showerror(f"Chem2Line - {lang_dict.get('error_unknown_title', '未知错误')}", f"{lang_dict.get('error_code', '错误代码')}: 1000\n{lang_dict.get('error_loading_database', '无法加载数据库')}: {e}")
            progress_bar.stop()  # 停止进度条
            progress_window.withdraw()  # 隐藏进度条窗口

    # 创建一个新的进度条窗口
    progress_window = tk.Toplevel()
    progress_window.title(f"Chem2Line - {lang_dict.get('loading_database_title', '加载数据库中')}")
    progress_window.geometry("300x100")
    progress_window.iconbitmap("lib/media/nctl.ico")
    progress_window.withdraw()  # 初始时隐藏窗口
    progress_window.attributes("-toolwindow", 2)
    progress_label = tk.Label(progress_window, text=lang_dict.get("loading_database_message", "加载中，请稍候..."), font=("Arial", 12))
    progress_label.pack(pady=10)

    progress_bar = ttk.Progressbar(progress_window, orient="horizontal", length=250, mode="indeterminate")
    progress_bar.pack(pady=10)

    threading.Thread(target=load_database).start()  # 启动加载数据库的后台线程

def on_submit():
    input_text = formula_entry.get().strip()
    if not input_text:
        messagebox.showwarning(f"Chem2Line - {lang_dict.get('input_empty_title', '输入为空')}", lang_dict.get("input_empty_message", "请输入化学式或 SMILES"))
        return

    try:
        smiles_list = get_smiles_options(input_text, smiles_dict)  # 添加 smiles_dict 参数
        if not smiles_list:
            raise ValueError(f"找不到 {input_text} 的 SMILES 表示，请检查输入或更换数据库")

        if len(smiles_list) == 1:
            selected_smiles = smiles_list[0]
        else:
            selected_smiles = show_smiles_selection(smiles_list)

        if selected_smiles:
            img = formula_to_bondline(selected_smiles)
            img = ImageTk.PhotoImage(img)

            result_label.config(image=img)
            result_label.image = img

            # 显示3D视图按钮
            view_3d_button.config(state=tk.NORMAL, command=lambda: show_3d_viewer(selected_smiles))
            analyze_button.config(state=tk.NORMAL)
            add_to_history(selected_smiles)

    except ValueError as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_not_found_title', '未找到结果')}", f"{lang_dict.get('error_code', '错误代码')}: 1001\n{str(e)}")
    except Exception as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_unknown_title', '未知错误')}", f"{lang_dict.get('error_code', '错误代码')}: 1000\n{lang_dict.get('error_unknown_message', '发生未知错误')}: {e}")

def on_analyze():
    input_text = formula_entry.get().strip()
    if not input_text:
        messagebox.showwarning(f"Chem2Line - {lang_dict.get('input_empty_title', '输入为空')}", lang_dict.get("input_empty_message", "请输入化学式或 SMILES"))
        return

    try:
        smiles_list = get_smiles_options(input_text, smiles_dict)  # 添加 smiles_dict 参数

        if not smiles_list:
            raise ValueError(f"{lang_dict.get('error_not_found', '找不到')} {input_text} {lang_dict.get('smiles_representation', '的 SMILES 表示，请检查输入或更换数据库')}")

        if len(smiles_list) == 1:
            selected_smiles = smiles_list[0]
        else:
            selected_smiles = show_smiles_selection(smiles_list)

        if selected_smiles:
            properties = analyze_molecule(selected_smiles, lang_dict)
            properties_str = "\n".join([f"{key}: {value}" for key, value in properties.items()])
            messagebox.showinfo(f"Chem2Line - {lang_dict.get('molecule_analysis', '分子分析')}", properties_str)
            add_to_history(selected_smiles)

    except ValueError as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_not_found_title', '未找到结果')}", f"{lang_dict.get('error_code', '错误代码')}: 1001\n{str(e)}")
    except RuntimeError as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_generation_failed_title', '生成失败')}", f"{lang_dict.get('error_code', '错误代码')}: 1002\n{str(e)}")
    except Exception as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_unknown_title', '未知错误')}", f"{lang_dict.get('error_code', '错误代码')}: 1000\n{lang_dict.get('error_unknown_message', '发生未知错误')}: {e}")

def save_image_as_png():
    if result_label.image:
        file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")])
        if file_path:
            result_label.image._PhotoImage__photo.write(file_path)
            messagebox.showinfo(f"Chem2Line - {lang_dict.get('save_success_title', '保存成功')}", f"{lang_dict.get('save_success_message', '图像已保存到')} {file_path}")

def save_image_as_svg():
    input_text = formula_entry.get().strip()
    if not input_text:
        messagebox.showwarning(f"Chem2Line - {lang_dict.get('input_empty_title', '输入为空')}", lang_dict.get("input_empty_message", "请输入化学式或 SMILES"))
        return

    try:
        smiles_list = get_smiles_options(input_text, smiles_dict)  # 添加 smiles_dict 参数
        if not smiles_list:
            raise ValueError(f"找不到 {input_text} 的 SMILES 表示，请检查输入或更换数据库")

        if len(smiles_list) == 1:
            selected_smiles = smiles_list[0]
        else:
            selected_smiles = show_smiles_selection(smiles_list)

        if selected_smiles:
            mol = Chem.MolFromSmiles(selected_smiles)
            if mol is None:
                raise ValueError(f"无效的SMILES字符串: {selected_smiles}")

            file_path = filedialog.asksaveasfilename(defaultextension=".svg", filetypes=[("SVG files", "*.svg")])
            if file_path:
                Draw.MolToFile(mol, file_path, imageType="svg")
                messagebox.showinfo(f"Chem2Line - {lang_dict.get('save_success_title', '保存成功')}", f"{lang_dict.get('save_success_message', '图像已保存到')} {file_path}")

    except ValueError as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_not_found_title', '未找到结果')}", f"{lang_dict.get('error_code', '错误代码')}: 1001\n{str(e)}")
    except RuntimeError as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_generation_failed_title', '生成失败')}", f"{lang_dict.get('error_code', '错误代码')}: 1002\n{str(e)}")
    except Exception as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_unknown_title', '未知错误')}", f"{lang_dict.get('error_code', '错误代码')}: 1000\n{lang_dict.get('error_unknown_message', '发生未知错误')}: {e}")

def change_database():
    load_database_with_progress()  # 使用新的加载数据库方法

def show_database_info():
    info = get_database_info(database_path)
    info_str = "\n".join([f"{lang_dict.get(key, key)}: {value}" for key, value in info.items()])
    messagebox.showinfo(f"Chem2Line - {lang_dict.get('database_info_title', '数据库信息')}", info_str)

def show_about_developer():
    messagebox.showinfo(f"Chem2Line - {lang_dict.get('about_developer_title', '关于开发者')}", f"{lang_dict.get('developer', '开发者')}: {DEVELOPER}\n{lang_dict.get('version', '版本')}: {VERSION}\n{lang_dict.get('date', '日期')}: {DATE}\n{lang_dict.get('core_version', '内核版本')}: {CORE_VERSION}")

def show_repository():
    messagebox.showinfo(f"Chem2Line - {lang_dict.get('repository_title', '软件仓库')}", "GitHub: https://github.com/Ziyang-Bai/Chem2Line")
def show_about_developer_with_icon():
    # 创建关于窗口
    about_window = tk.Toplevel(root)
    about_window.title(f"Chem2Line - {lang_dict.get('about_developer_title', '关于开发者')}")
    about_window.geometry("300x400")
    about_window.iconbitmap("lib/media/nctl.ico")

    # 加载图片
    try:
        icon_image = tk.PhotoImage(file="lib/media/chem2line.png")  # 替换为图标的路径
        icon_label = tk.Label(about_window, image=icon_image)
        icon_label.image = icon_image  # 保存引用，防止被垃圾回收
        icon_label.pack(pady=10)
    except Exception as e:
        tk.Label(about_window, text=f"{lang_dict.get('error_loading_icon', '无法加载图标')}: {e}").pack(pady=10)

    # 添加文字信息
    info_text = f"""
    {lang_dict.get('developer', '开发者')}: {DEVELOPER}
    {lang_dict.get('version', '版本')}: {VERSION}
    {lang_dict.get('date', '日期')}: {DATE}
    {lang_dict.get('core_version', '内核版本')}: {CORE_VERSION}
    """
    info_label = tk.Label(about_window, text=info_text, font=("Arial", 12), justify="left")
    info_label.pack(pady=10)

    about_window.mainloop()

# 加载配置文件
def load_config():
    """
    加载配置文件
    """
    try:
        tree = ET.parse('lib/config/config.xml')
        root = tree.getroot()
        config = {child.tag: child.text for child in root if child.tag != 'available_languages'}
        config['available_languages'] = [lang.text for lang in root.find('available_languages')]
        # 默认启用记录历史记录
        config['record_history'] = config.get('record_history', 'true').lower() == 'true'
        return config
    except Exception as e:
        messagebox.showerror(lang_dict.get("config_error_title", "配置错误"), f"{lang_dict.get('error_code', '错误代码')}: 2000\n{lang_dict.get('config_error_message', '无法加载配置文件')}: {e}")
        return {'record_history': True}

# 保存配置文件
def save_config(config):
    """
    保存配置文件
    """
    try:
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
    except Exception as e:
        messagebox.showerror(lang_dict.get("config_error_title", "配置错误"), f"{lang_dict.get('error_code', '错误代码')}: 2001\n{lang_dict.get('config_save_error_message', '无法保存配置文件')}: {e}")

# 加载语言文件
def load_language(lang_file):
    try:
        tree = ET.parse(f'{lang_file}')
        root = tree.getroot()
        lang_dict = {child.tag: child.text for child in root}
        return lang_dict
    except Exception as e:
        # 返回一个默认的空字典
        return {
            "language_error_title": "语言错误",
            "error_code": "错误代码",
            "language_error_message": "无法加载语言文件"
        }

# 切换语言
def change_language(lang):
    config['language'] = lang
    save_config(config)
    lang_dict = load_language(f'lib/lang/{lang}.xml')
    messagebox.showinfo(lang_dict.get("restart_title", "重启应用"), lang_dict.get("restart_message", "语言已更改，请重启应用以应用更改。"))
    root.quit()

def show_config_window():
    """
    显示配置窗口
    """
    config_window = Toplevel(root)
    config_window.title(lang_dict.get("config_title", "配置"))
    config_window.geometry("400x300")
    config_window.iconbitmap("lib/media/nctl.ico")

    # 语言配置
    lang_label = tk.Label(config_window, text=lang_dict.get("select_language") or "选择语言：", font=("Arial", 12))
    lang_label.pack(pady=10)
    lang_var = StringVar(config_window)
    lang_val = config.get('language', 'en_us')
    if not isinstance(lang_val, str):
        lang_val = 'en_us'
    lang_var.set(lang_val)
    available_langs = config.get('available_languages', ['en_us'])
    if not isinstance(available_langs, list):
        available_langs = ['en_us']
    lang_options = [load_language(f'lib/lang/{lang}.xml').get("language_name", lang) for lang in available_langs]
    lang_menu = OptionMenu(config_window, lang_var, *lang_options)
    lang_menu.pack(pady=10)

    # 默认数据库配置
    db_label = tk.Label(config_window, text=lang_dict.get("select_database") or "选择默认数据库：", font=("Arial", 12))
    db_label.pack(pady=10)
    db_var = StringVar(config_window)
    db_val = config.get('default_database', database_path)
    if not isinstance(db_val, str):
        db_val = database_path
    db_var.set(os.path.basename(db_val))
    db_files = [f for f in os.listdir('lib/db') if f.endswith('.xml')]
    db_menu = OptionMenu(config_window, db_var, *db_files)
    db_menu.pack(pady=10)

    # 是否记录历史记录
    record_history_val = config.get('record_history', True)
    if not isinstance(record_history_val, bool):
        record_history_val = str(record_history_val).lower() == 'true'
    record_history_var = tk.BooleanVar(value=record_history_val)
    record_history_check = tk.Checkbutton(config_window, text=lang_dict.get("record_history") or "记录历史记录", variable=record_history_var)
    record_history_check.pack(pady=10)

    # 3D模型配置
    model3d_label = tk.Label(config_window, text=lang_dict.get("select_3d_model") or "选择3D模型：", font=("Arial", 12))
    model3d_label.pack(pady=10)
    model3d_var_cfg = StringVar(config_window)
    model3d_val = config.get('default_3d_model', 'ball_and_stick')
    if not isinstance(model3d_val, str):
        model3d_val = 'ball_and_stick'
    model3d_var_cfg.set(model3d_val)
    model3d_menu_cfg = OptionMenu(config_window, model3d_var_cfg, 'ball_and_stick', 'space_filling')
    model3d_menu_cfg.pack(pady=10)

    def save_config_changes():
        # 语言
        selected_lang = 'en_us'
        for idx, lang in enumerate(available_langs):
            if lang_options[idx] == lang_var.get():
                selected_lang = lang
                break
        config['language'] = selected_lang
        # 数据库
        config['default_database'] = f'lib/db/{db_var.get()}'
        # 历史
        config['record_history'] = record_history_var.get()
        # 3D模型
        config['default_3d_model'] = model3d_var_cfg.get()
        save_config(config)
        messagebox.showinfo(lang_dict.get("config_saved_title") or "配置已保存", lang_dict.get("config_saved_message") or "配置已保存，请重启应用以应用更改。")
        config_window.destroy()

    save_button = tk.Button(config_window, text=lang_dict.get("save_button") or "保存", command=save_config_changes)
    save_button.pack(pady=20)

# 初始化配置和语言
config = load_config()
language = config.get('language', 'en_us')
lang_dict = load_language(f'lib/lang/{language}.xml')
database_path = config.get('default_database', 'lib/db/default_database.xml')

# 初始化主窗口
root = tk.Tk()
root.title("Chem2Line")
root.geometry("800x600")
root.iconbitmap("lib/media/nctl.ico")

# 检查是否启用调试模式
if "--debug" in sys.argv:
    enable_debug_mode(root, globals())

output_type = StringVar(value="bondline")

# 创建菜单栏
menu_bar = Menu(root)
root.config(menu=menu_bar)

# 文件菜单
file_menu = Menu(menu_bar, tearoff=0)
save_image_menu = Menu(file_menu, tearoff=0)
save_image_menu.add_command(label=lang_dict.get("save_as_png", "保存为PNG"), command=save_image_as_png)
save_image_menu.add_command(label=lang_dict.get("save_as_svg", "保存为SVG"), command=save_image_as_svg)
file_menu.add_cascade(label=lang_dict.get("save_image", "保存图像"), menu=save_image_menu)
file_menu.add_command(label=lang_dict.get("config", "配置"), command=show_config_window)
file_menu.add_separator()
file_menu.add_command(label=lang_dict.get("view_long_history", "历史记录"), command=show_long_history)
file_menu.add_command(label=lang_dict.get("exit", "退出"), command=root.quit)
menu_bar.add_cascade(label=lang_dict.get("file", "文件"), menu=file_menu)

# 输出菜单
output_menu = Menu(menu_bar, tearoff=0)
output_menu.add_radiobutton(label=lang_dict.get("bondline_output", "键线式输出"), variable=output_type, value="bondline")
menu_bar.add_cascade(label=lang_dict.get("output", "输出"), menu=output_menu)

# 数据库菜单
database_menu = Menu(menu_bar, tearoff=0)
database_menu.add_command(label=lang_dict.get("change_database", "更换数据库"), command=change_database)
database_menu.add_command(label=lang_dict.get("database_info", "关于数据库"), command=show_database_info)

# 常用数据库子菜单
common_db_menu = Menu(database_menu, tearoff=0)
for db_file in os.listdir('lib/db'):
    if db_file.endswith('.xml'):
        common_db_menu.add_command(label=db_file, command=lambda f=db_file: load_database_with_progress(f'lib/db/{f}'))
database_menu.add_cascade(label=lang_dict.get("common_databases", "常用数据库"), menu=common_db_menu)

menu_bar.add_cascade(label=lang_dict.get("database", "数据库"), menu=database_menu)

# 3D模型菜单
model3d_menu = Menu(menu_bar, tearoff=0)
def get_default_3d_model():
    val = config.get('default_3d_model', 'ball_and_stick')
    if not isinstance(val, str):
        return 'ball_and_stick'
    return val
model3d_var = StringVar(value=get_default_3d_model())
model3d_menu.add_radiobutton(label=(lang_dict.get('ball_and_stick') or '球棍模型'), variable=model3d_var, value='ball_and_stick')
model3d_menu.add_radiobutton(label=(lang_dict.get('space_filling') or '比例模型'), variable=model3d_var, value='space_filling')
menu_bar.add_cascade(label=(lang_dict.get('model3d_menu') or '3D模型'), menu=model3d_menu)

# 创建工具栏菜单
overlay_menu = Menu(menu_bar, tearoff=0)
overlay_menu.add_radiobutton(label=lang_dict.get("no_overlay", "无"), variable=output_type, value="none")
overlay_menu.add_radiobutton(label=lang_dict.get("force_field", "力场"), variable=output_type, value="force_field", command=on_overlay_force_field)
menu_bar.add_cascade(label=lang_dict.get("overlay", "叠加显示"), menu=overlay_menu)

# 创建历史记录菜单
history_menu = Menu(menu_bar, tearoff=0)

menu_bar.add_cascade(label=lang_dict.get("history", "历史记录"), menu=history_menu)

# 在菜单栏中添加查询工具
#tools_menu = Menu(menu_bar, tearoff=0)
#tools_menu.add_command(label=lang_dict.get("search_database", "查询数据库"), command=search_database)
#menu_bar.add_cascade(label=lang_dict.get("tools", "工具"), menu=tools_menu)

# 加载 SMILES 数据库
smiles_dict = load_smiles_database(database_path)

# 加载历史记录
history = load_history()

# 创建输入框和标签
input_label = tk.Label(root, text=lang_dict.get("input_label", "请输入化学式或 SMILES："), font=("Arial", 14))
input_label.pack(pady=10)

formula_entry = tk.Entry(root, font=("Arial", 14), width=30)
formula_entry.pack(pady=10)
update_history_menu(history_menu, history)

# 创建按钮框架
button_frame = tk.Frame(root)
button_frame.pack(pady=10)

# 创建提交按钮
submit_button = tk.Button(button_frame, text=lang_dict.get("submit_button", "生成键线式"), font=("Arial", 14), command=on_submit)
submit_button.grid(row=0, column=0, padx=5)

# 创建显示3D视图按钮
view_3d_button = tk.Button(button_frame, text=lang_dict.get("view_3d_button", "显示3D视图"), font=("Arial", 14), state=tk.DISABLED)
view_3d_button.grid(row=0, column=1, padx=5)

def show_3d_viewer(smiles):
    from lib.ctlcore import MoleculeViewer
    win = tk.Toplevel(root)
    win.title("3D Viewer")
    viewer = MoleculeViewer(win, smiles, model_type=model3d_var.get(), lang_dict=lang_dict)
    win.mainloop()

# 创建分析按钮
analyze_button = tk.Button(button_frame, text=lang_dict.get("analyze_button", "分析分子"), font=("Arial", 14), command=on_analyze, state=tk.DISABLED)
analyze_button.grid(row=0, column=2, padx=5)

# 创建显示化学信息按钮
#info_button = tk.Button(button_frame, text=lang_dict.get("info_button", "显示化学信息"), font=("Arial", 14), command=show_chemical_info)
#info_button.grid(row=0, column=3, padx=5)
# 关于菜单
about_menu = Menu(menu_bar, tearoff=0)
about_menu.add_command(label=lang_dict.get("developer", "开发者"), command=show_about_developer)
about_menu.add_command(label=lang_dict.get("repository", "软件仓库"), command=show_repository)
menu_bar.add_cascade(label=lang_dict.get("about", "关于"), menu=about_menu)
def show_chemical_info():
    input_text = formula_entry.get().strip()
    if not input_text:
        messagebox.showwarning(f"Chem2Line - {lang_dict.get('input_empty_title', '输入为空')}", lang_dict.get("input_empty_message", "请输入化学式或 SMILES"))
        return

    try:
        smiles_list = get_smiles_options(input_text, smiles_dict)

        if not smiles_list:
            raise ValueError(f"{lang_dict.get('error_not_found', '找不到')} {input_text} {lang_dict.get('smiles_representation', '的 SMILES 表示，请检查输入或更换数据库')}")

        if len(smiles_list) == 1:
            selected_smiles = smiles_list[0]
        else:
            selected_smiles = show_smiles_selection(smiles_list)

        if selected_smiles:
            info_str = show_chemical_info(selected_smiles)
            messagebox.showinfo(f"Chem2Line - {lang_dict.get('chemical_info', '化学信息')}", info_str)
            add_to_history(selected_smiles)

    except ValueError as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_not_found_title', '未找到结果')}", f"{lang_dict.get('error_code', '错误代码')}: 1001\n{str(e)}")
    except RuntimeError as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_generation_failed_title', '生成失败')}", f"{lang_dict.get('error_code', '错误代码')}: 1002\n{str(e)}")
    except Exception as e:
        messagebox.showerror(f"Chem2Line - {lang_dict.get('error_unknown_title', '未知错误')}", f"{lang_dict.get('error_code', '错误代码')}: 1000\n{lang_dict.get('error_unknown_message', '发生未知错误')}: {e}")

# 显示生成结果的标签
result_label = tk.Label(root)
result_label.pack(pady=20)

def search_database():
    """
    在已加载的数据库中查询化学物质
    """
    search_window = Toplevel(root)
    search_window.title(lang_dict.get("search_database_title", "查询数据库"))
    search_window.geometry("600x400")
    search_window.iconbitmap("lib/media/nctl.ico")

    def perform_search():
        query = search_entry.get().strip()
        if not query:
            messagebox.showwarning(f"Chem2Line - {lang_dict.get('input_empty_title', '输入为空')}", lang_dict.get("input_empty_message", "请输入查询内容"))
            return

        results = search_database(query, smiles_dict)

        for row in tree.get_children():
            tree.delete(row)

        for result in results:
            tree.insert("", "end", values=result)

        if not results:
            messagebox.showinfo(f"Chem2Line - {lang_dict.get('no_results_title', '无结果')}", lang_dict.get("no_results_message", "未找到匹配的化学物质"))

    # 搜索框
    search_label = tk.Label(search_window, text=lang_dict.get("search_label", "请输入名称或化学式："), font=("Arial", 12))
    search_label.pack(pady=10)

    search_entry = tk.Entry(search_window, font=("Arial", 12), width=40)
    search_entry.pack(pady=10)

    search_button = tk.Button(search_window, text=lang_dict.get("search_button", "查询"), command=perform_search)
    search_button.pack(pady=10)

    # 结果表格
    columns = ("name", "formula", "smiles")
    tree = ttk.Treeview(search_window, columns=columns, show="headings")
    tree.heading("name", text=lang_dict.get("name", "名称"))
    tree.heading("formula", text=lang_dict.get("formula", "化学式"))
    tree.heading("smiles", text="SMILES")
    tree.pack(fill="both", expand=True)

def advanced_search_gui():
    """
    高级查询界面，允许用户输入查询条件并选择查询类型
    """
    search_window = Toplevel(root)
    search_window.title(lang_dict.get("advanced_search_title", "高级查询"))
    search_window.geometry("600x400")
    search_window.iconbitmap("lib/media/nctl.ico")

    # 查询类型选项
    search_type_label = tk.Label(search_window, text=lang_dict.get("search_type_label", "选择查询类型："), font=("Arial", 12))
    search_type_label.pack(pady=5)

    search_type_var = StringVar(value="name")  # 默认查询类型为名称
    search_type_menu = ttk.Combobox(search_window, textvariable=search_type_var, state="readonly", values=[
        lang_dict.get("search_by_name", "名称"),
        lang_dict.get("search_by_formula", "化学式")
    ])
    search_type_menu.pack(pady=5)

    # 查询输入框
    search_label = tk.Label(search_window, text=lang_dict.get("search_label", "请输入查询内容："), font=("Arial", 12))
    search_label.pack(pady=5)

    search_entry = tk.Entry(search_window, font=("Arial", 12), width=40)
    search_entry.pack(pady=5)

    # 结果表格
    columns = ("name", "formula", "smiles")
    tree = ttk.Treeview(search_window, columns=columns, show="headings")
    tree.heading("name", text=lang_dict.get("name", "名称"))
    tree.heading("formula", text=lang_dict.get("formula", "化学式"))
    tree.heading("smiles", text="SMILES")
    tree.pack(fill="both", expand=True, pady=10)

    def perform_search():
        """
        执行查询操作，根据用户输入和查询类型筛选结果
        """
        query = search_entry.get().strip().lower()
        if not query:
            messagebox.showwarning(f"Chem2Line - {lang_dict.get('input_empty_title', '输入为空')}", lang_dict.get("input_empty_message", "请输入查询内容"))
            return

        search_type = "name" if search_type_var.get() == lang_dict.get("search_by_name", "名称") else "formula"
        results = []

        for smiles, info in smiles_dict.items():
            # 确保 info 是字典类型
            if isinstance(info, dict) and query in info.get(search_type, "").lower():
                results.append((info.get("name", "N/A"), info.get("formula", "N/A"), smiles))

        # 清空表格并插入新结果
        for row in tree.get_children():
            tree.delete(row)

        for result in results:
            tree.insert("", "end", values=result)

        if not results:
            messagebox.showinfo(f"Chem2Line - {lang_dict.get('no_results_title', '无结果')}", lang_dict.get("no_results_message", "未找到匹配的化学物质"))

    def reset_search():
        """
        重置查询界面，清空输入框和表格
        """
        search_entry.delete(0, tk.END)
        for row in tree.get_children():
            tree.delete(row)

    # 查询按钮
    search_button = tk.Button(search_window, text=lang_dict.get("search_button", "查询"), command=perform_search)
    search_button.pack(side="left", padx=10, pady=10)

    # 重置按钮
    reset_button = tk.Button(search_window, text=lang_dict.get("reset_button", "重新开始"), command=reset_search)
    reset_button.pack(side="right", padx=10, pady=10)

# 在工具菜单中添加高级查询功能
#tools_menu.add_command(label=lang_dict.get("advanced_search", "高级查询"), command=advanced_search_gui)

def sdf_converter_gui():
    """
    SDF 转换工具的 GUI 界面
    """
    converter_window = Toplevel(root)
    converter_window.title(lang_dict.get("sdf_converter_title", "SDF 转换工具"))
    converter_window.geometry("400x300")
    converter_window.iconbitmap("lib/media/nctl.ico")

    def convert_sdf_to_smiles():
        """
        将 SDF 文件转换为 SMILES 格式
        """
        sdf_file = filedialog.askopenfilename(filetypes=[("SDF files", "*.sdf")])
        if not sdf_file:
            return

        try:
            output_file = filedialog.asksaveasfilename(defaultextension=".txt", filetypes=[("Text files", "*.txt")])
            if not output_file:
                return

            with open(sdf_file, "r") as sdf, open(output_file, "w") as output:
                suppl = Chem.SDMolSupplier(sdf_file)
                for mol in suppl:
                    if mol is not None:
                        smiles = Chem.MolToSmiles(mol)
                        output.write(smiles + "\n")

            messagebox.showinfo(lang_dict.get("conversion_success_title", "转换成功"), lang_dict.get("conversion_success_message", "SDF 文件已成功转换为 SMILES 格式"))
        except Exception as e:
            messagebox.showerror(lang_dict.get("conversion_error_title", "转换错误"), f"{lang_dict.get('error_code', '错误代码')}: 3000\n{lang_dict.get('conversion_error_message', '无法转换文件')}: {e}")

    # 添加按钮
    convert_button = tk.Button(converter_window, text=lang_dict.get("convert_button", "选择 SDF 文件并转换"), command=convert_sdf_to_smiles)
    convert_button.pack(pady=20)

# 在工具菜单中添加 SDF 转换工具
tools_menu = Menu(menu_bar, tearoff=0)
tools_menu.add_command(label=lang_dict.get("sdf_converter", "SDF 转换工具"), command=sdf_converter_gui)
menu_bar.add_cascade(label=lang_dict.get("tools", "工具"), menu=tools_menu)

# 定义 update_menu_callback 函数
def update_menu_callback(history):
    """
    更新历史记录菜单的回调函数
    """
    update_history_menu(history_menu, history, formula_entry)

# 修复 formula_entry 的引用
update_history_menu(history_menu, history)

# 运行主窗口
root.mainloop()
