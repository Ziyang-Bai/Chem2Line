# ctlgui.py
import tkinter as tk
from tkinter import filedialog, messagebox, Scrollbar, Canvas, Menu, ttk
from PIL import Image, ImageTk
from ctlcore import load_smiles_database, get_smiles_options, formula_to_bondline, get_database_info, core_version, load_smiles_database_with_progress
import time
VERSION = "1.0"
DEVELOPER = "Ziyang-Bai"
DATE = "2025-01-01"
CORE_VERSION = core_version()
def show_smiles_selection(smiles_list):
    """
    弹出窗口显示可供选择的 SMILES
    :param smiles_list: 包含多个 SMILES 的列表
    :return: 用户选择的 SMILES
    """
    selection_window = tk.Toplevel()
    selection_window.title("选择 SMILES 表示")
    selection_window.geometry("600x400")

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
            tk.Label(scrollable_frame, text=f"无法加载 {smiles}: {e}").pack()

    # 等待用户选择
    selection_window.wait_window()
    return selected_smiles.get()


def show_progress_window():
    """
    显示进度条窗口
    :return: 更新进度的回调函数
    """
    progress_window = tk.Toplevel()
    progress_window.title("加载数据库中")
    progress_window.geometry("400x100")
    progress_window.resizable(False, False)

    label = tk.Label(progress_window, text="加载中，请稍候...", font=("Arial", 12))
    label.pack(pady=10)

    progress_bar = ttk.Progressbar(progress_window, orient="horizontal", length=300, mode="determinate")
    progress_bar.pack(pady=10)

    def update_progress(current, total):
        progress_bar["value"] = (current / total) * 100
        progress_window.update()

        if current == total:  # 完成时关闭窗口
            progress_window.destroy()

    return update_progress

def change_database():
    global smiles_dict
    file_path = filedialog.askopenfilename(filetypes=[("XML files", "*.xml")])
    if file_path:
        try:
            update_progress = show_progress_window_with_estimation()  # 显示进度条
            new_smiles_dict = load_smiles_database_with_progress(file_path, update_progress)
            smiles_dict = new_smiles_dict  # 更新全局数据库
            database_info = get_database_info(file_path)
            messagebox.showinfo("数据库已更换", f"当前数据库信息:\n{database_info}")
        except FileNotFoundError as e:
            messagebox.showerror("错误", str(e))
        except ValueError as e:
            messagebox.showerror("数据库格式错误", str(e))
        except RuntimeError as e:
            messagebox.showerror("加载错误", str(e))
        except Exception as e:
            messagebox.showerror("未知错误", f"无法加载数据库: {e}")

def save_image():
    if result_label.image:
        file_path = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG files", "*.png")])
        if file_path:
            result_label.image._PhotoImage__photo.write(file_path)
            messagebox.showinfo("保存成功", f"键线式图像已保存到 {file_path}")


def on_submit():
    input_text = formula_entry.get().strip()
    if not input_text:
        messagebox.showwarning("输入为空", "请输入化学式或 SMILES")
        return

    try:
        smiles_list = get_smiles_options(input_text, smiles_dict)

        if not smiles_list:
            raise ValueError(f"找不到 {input_text} 的 SMILES 表示，请检查输入或更换数据库")

        if len(smiles_list) == 1:
            selected_smiles = smiles_list[0]
        else:
            selected_smiles = show_smiles_selection(smiles_list)

        if selected_smiles:
            try:
                img = formula_to_bondline(selected_smiles)
                img = ImageTk.PhotoImage(img)

                result_label.config(image=img)
                result_label.image = img
            except Exception as e:
                raise RuntimeError(f"无法生成键线式: {e}")

    except ValueError as e:
        messagebox.showerror("未找到结果", str(e))
    except RuntimeError as e:
        messagebox.showerror("生成失败", str(e))
    except Exception as e:
        messagebox.showerror("未知错误", f"发生未知错误: {e}")
def show_progress_window_with_estimation():
    """
    显示带时间估算的进度条窗口
    :return: 更新进度的回调函数
    """
    progress_window = tk.Toplevel()
    progress_window.title("加载数据库中")
    progress_window.geometry("400x150")
    progress_window.resizable(False, False)

    label = tk.Label(progress_window, text="加载中，请稍候...", font=("Arial", 12))
    label.pack(pady=10)

    progress_bar = ttk.Progressbar(progress_window, orient="horizontal", length=300, mode="determinate")
    progress_bar.pack(pady=10)

    time_label = tk.Label(progress_window, text="", font=("Arial", 10), fg="blue")
    time_label.pack(pady=5)

    start_time = time.time()

    def update_progress(current, total):
        progress_bar["value"] = (current / total) * 100
        elapsed_time = time.time() - start_time

        # 估算剩余时间
        if current > 0:
            estimated_total_time = elapsed_time / current * total
            remaining_time = estimated_total_time - elapsed_time
            time_label.config(text=f"预计剩余时间: {int(remaining_time)} 秒")
        else:
            time_label.config(text="预计剩余时间: 计算中...")

        progress_window.update()

        if current == total:  # 完成时关闭窗口
            progress_window.destroy()

    return update_progress


def show_database_info():
    info = get_database_info()
    info_str = "\n".join([f"{key}: {value}" for key, value in info.items()])
    messagebox.showinfo("数据库信息", info_str)

def show_about_developer():
    messagebox.showinfo("关于开发者", "开发者: " +  DEVELOPER + "\n版本: " + VERSION +"\n日期: " + DATE + "\n" + "内核版本: " + CORE_VERSION)

def show_repository():
    messagebox.showinfo("软件仓库", "GitHub: https://github.com/Ziyang-Bai/Chem2Line")

# 初始化主窗口
root = tk.Tk()
root.title("Chem2Line")
root.geometry("800x600")

# 创建菜单栏
menu_bar = Menu(root)
root.config(menu=menu_bar)

# 文件菜单
file_menu = Menu(menu_bar, tearoff=0)
file_menu.add_command(label="保存键线式图像", command=save_image)
file_menu.add_separator()
file_menu.add_command(label="退出", command=root.quit)
menu_bar.add_cascade(label="文件", menu=file_menu)

# 数据库菜单
database_menu = Menu(menu_bar, tearoff=0)
database_menu.add_command(label="更换数据库", command=change_database)
database_menu.add_command(label="关于数据库", command=show_database_info)
menu_bar.add_cascade(label="数据库", menu=database_menu)

# 关于菜单
about_menu = Menu(menu_bar, tearoff=0)
about_menu.add_command(label="开发者", command=show_about_developer)
about_menu.add_command(label="软件仓库", command=show_repository)
menu_bar.add_cascade(label="关于", menu=about_menu)

# 加载 SMILES 数据库
smiles_dict = load_smiles_database("default_database.xml")

# 创建输入框和标签
input_label = tk.Label(root, text="请输入化学式或 SMILES：", font=("Arial", 14))
input_label.pack(pady=10)

formula_entry = tk.Entry(root, font=("Arial", 14), width=30)
formula_entry.pack(pady=10)

# 创建提交按钮
submit_button = tk.Button(root, text="生成键线式", font=("Arial", 14), command=on_submit)
submit_button.pack(pady=20)

# 显示生成结果的标签
result_label = tk.Label(root)
result_label.pack(pady=20)

# 运行主窗口
root.mainloop()
