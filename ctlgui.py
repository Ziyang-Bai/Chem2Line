"""Chem2Line - 桌面 GUI 入口

将化学式或 SMILES 转换为键线式图像的图形界面工具。
"""

import os
import sys
import threading
import tkinter as tk
from tkinter import filedialog, messagebox, Menu, ttk, Toplevel, StringVar, OptionMenu, Canvas, Scrollbar

from PIL import Image, ImageTk
from rdkit import Chem
from rdkit.Chem import Draw

from lib.version import VERSION, CORE_VERSION, DEVELOPER, DATE, APP_NAME, REPOSITORY
from lib.config import load_config, save_config
from lib.i18n import load_language, get_text
from lib.database import load_smiles_database, get_smiles_options, get_database_info
from lib.chemistry import formula_to_bondline, overlay_force_field, analyze_molecule
from lib.history import (
    load_history, save_history, add_to_history as _add_to_history, clear_history,
)
from lib.viewer3d import MoleculeViewer
from lib.debug import enable_debug_mode


class Chem2LineApp:
    """Chem2Line 主应用程序类"""

    def __init__(self):
        # 加载配置和语言
        self.config = load_config()
        language = self.config.get("language", "zh_cn")
        self.lang = load_language(language)
        self.database_path = self.config.get("default_database", "lib/db/default_database.xml")

        # 加载数据
        self.smiles_dict = load_smiles_database(self.database_path)
        self.history = load_history()

        # 创建主窗口
        self.root = tk.Tk()
        self.root.title(APP_NAME)
        self.root.geometry("800x600")
        self.root.iconbitmap("lib/media/nctl.ico")

        # 调试模式
        if "--debug" in sys.argv:
            enable_debug_mode(self.root, globals())

        # UI 状态变量
        self.output_type = StringVar(value="bondline")
        self.model3d_var = StringVar(value=self.config.get("default_3d_model", "ball_and_stick"))

        # 构建界面
        self._build_menu()
        self._build_main_ui()

        # 更新历史菜单
        self._update_history_menu()

    def run(self):
        """启动主循环"""
        self.root.mainloop()

    # ─── 辅助方法 ─────────────────────────────────────────────

    def _t(self, key: str, default: str = "") -> str:
        """获取翻译文本的快捷方法"""
        return get_text(self.lang, key, default)

    def _add_history(self, smiles: str):
        """添加到历史记录"""
        if not self.config.get("record_history", True):
            return
        input_text = self.formula_entry.get().strip()
        self.history = _add_to_history(smiles, input_text, self.history)
        self._update_history_menu()

    # ─── 菜单构建 ─────────────────────────────────────────────

    def _build_menu(self):
        menu_bar = Menu(self.root)
        self.root.config(menu=menu_bar)

        # 文件菜单
        file_menu = Menu(menu_bar, tearoff=0)
        save_menu = Menu(file_menu, tearoff=0)
        save_menu.add_command(label=self._t("save_as_png", "保存为PNG"), command=self._save_as_png)
        save_menu.add_command(label=self._t("save_as_svg", "保存为SVG"), command=self._save_as_svg)
        file_menu.add_cascade(label=self._t("save_image", "保存图像"), menu=save_menu)
        file_menu.add_command(label=self._t("config_title", "配置"), command=self._show_config)
        file_menu.add_separator()
        file_menu.add_command(label=self._t("view_long_history", "历史记录"), command=self._show_long_history)
        file_menu.add_command(label=self._t("exit", "退出"), command=self.root.quit)
        menu_bar.add_cascade(label=self._t("file", "文件"), menu=file_menu)

        # 输出菜单
        output_menu = Menu(menu_bar, tearoff=0)
        output_menu.add_radiobutton(
            label=self._t("bondline_output", "键线式输出"),
            variable=self.output_type, value="bondline",
        )
        menu_bar.add_cascade(label=self._t("output", "输出"), menu=output_menu)

        # 数据库菜单
        db_menu = Menu(menu_bar, tearoff=0)
        db_menu.add_command(label=self._t("change_database", "更换数据库"), command=self._change_database)
        db_menu.add_command(label=self._t("database_info", "关于数据库"), command=self._show_database_info)
        common_db_menu = Menu(db_menu, tearoff=0)
        for db_file in os.listdir("lib/db"):
            if db_file.endswith(".xml"):
                common_db_menu.add_command(
                    label=db_file,
                    command=lambda f=db_file: self._load_database(f"lib/db/{f}"),
                )
        db_menu.add_cascade(label=self._t("common_databases", "常用数据库"), menu=common_db_menu)
        menu_bar.add_cascade(label=self._t("database", "数据库"), menu=db_menu)

        # 3D 模型菜单
        model_menu = Menu(menu_bar, tearoff=0)
        model_menu.add_radiobutton(
            label=self._t("ball_and_stick", "球棍模型"),
            variable=self.model3d_var, value="ball_and_stick",
        )
        model_menu.add_radiobutton(
            label=self._t("space_filling", "比例模型"),
            variable=self.model3d_var, value="space_filling",
        )
        menu_bar.add_cascade(label=self._t("model3d_menu", "3D模型"), menu=model_menu)

        # 叠加显示菜单
        overlay_menu = Menu(menu_bar, tearoff=0)
        overlay_menu.add_radiobutton(
            label=self._t("no_overlay", "无"), variable=self.output_type, value="none",
        )
        overlay_menu.add_radiobutton(
            label=self._t("force_field", "力场"),
            variable=self.output_type, value="force_field",
            command=self._on_overlay_force_field,
        )
        menu_bar.add_cascade(label=self._t("overlay", "叠加显示"), menu=overlay_menu)

        # 历史记录菜单
        self.history_menu = Menu(menu_bar, tearoff=0)
        menu_bar.add_cascade(label=self._t("history", "历史记录"), menu=self.history_menu)

        # 关于菜单
        about_menu = Menu(menu_bar, tearoff=0)
        about_menu.add_command(label=self._t("developer", "开发者"), command=self._show_about)
        about_menu.add_command(label=self._t("repository", "软件仓库"), command=self._show_repository)
        menu_bar.add_cascade(label=self._t("about", "关于"), menu=about_menu)

    # ─── 主界面构建 ───────────────────────────────────────────

    def _build_main_ui(self):
        # 输入标签
        input_label = tk.Label(
            self.root,
            text=self._t("input_label", "请输入化学式或 SMILES："),
            font=("Arial", 14),
        )
        input_label.pack(pady=10)

        # 输入框
        self.formula_entry = tk.Entry(self.root, font=("Arial", 14), width=30)
        self.formula_entry.pack(pady=10)

        # 按钮框架
        btn_frame = tk.Frame(self.root)
        btn_frame.pack(pady=10)

        self.submit_btn = tk.Button(
            btn_frame, text=self._t("submit_button", "生成键线式"),
            font=("Arial", 14), command=self._on_submit,
        )
        self.submit_btn.grid(row=0, column=0, padx=5)

        self.view_3d_btn = tk.Button(
            btn_frame, text=self._t("view_3d_button", "显示3D视图"),
            font=("Arial", 14), state=tk.DISABLED,
        )
        self.view_3d_btn.grid(row=0, column=1, padx=5)

        self.analyze_btn = tk.Button(
            btn_frame, text=self._t("analyze_button", "分析分子"),
            font=("Arial", 14), command=self._on_analyze, state=tk.DISABLED,
        )
        self.analyze_btn.grid(row=0, column=2, padx=5)

        # 结果显示区域
        self.result_label = tk.Label(self.root)
        self.result_label.pack(pady=20, expand=True)

    # ─── 核心操作 ─────────────────────────────────────────────

    def _get_selected_smiles(self) -> str:
        """获取用户输入并解析为 SMILES，必要时弹出选择窗口。"""
        input_text = self.formula_entry.get().strip()
        if not input_text:
            messagebox.showwarning(
                f"{APP_NAME} - {self._t('input_empty_title', '输入为空')}",
                self._t("input_empty_message", "请输入化学式或 SMILES"),
            )
            return ""

        smiles_list = get_smiles_options(input_text, self.smiles_dict)
        if not smiles_list:
            raise ValueError(
                f"{self._t('error_not_found', '找不到')} {input_text} "
                f"{self._t('smiles_representation', '的 SMILES 表示')}"
            )

        if len(smiles_list) == 1:
            return smiles_list[0]
        return self._show_smiles_selection(smiles_list)

    def _on_submit(self):
        """生成键线式按钮回调"""
        try:
            smiles = self._get_selected_smiles()
            if not smiles:
                return

            img = formula_to_bondline(smiles)
            tk_img = ImageTk.PhotoImage(img)
            self.result_label.config(image=tk_img)
            self.result_label.image = tk_img

            # 启用 3D 和分析按钮
            self.view_3d_btn.config(
                state=tk.NORMAL,
                command=lambda s=smiles: self._show_3d_viewer(s),
            )
            self.analyze_btn.config(state=tk.NORMAL)
            self._add_history(smiles)

        except ValueError as e:
            messagebox.showerror(
                f"{APP_NAME} - {self._t('error_not_found_title', '未找到结果')}",
                f"{self._t('error_code', '错误代码')}: 1001\n{e}",
            )
        except Exception as e:
            messagebox.showerror(
                f"{APP_NAME} - {self._t('error_unknown_title', '未知错误')}",
                f"{self._t('error_code', '错误代码')}: 1000\n{e}",
            )

    def _on_analyze(self):
        """分析分子按钮回调"""
        try:
            smiles = self._get_selected_smiles()
            if not smiles:
                return

            properties = analyze_molecule(smiles, self.lang)
            props_str = "\n".join(f"{k}: {v}" for k, v in properties.items())
            messagebox.showinfo(
                f"{APP_NAME} - {self._t('molecule_analysis', '分子分析')}",
                props_str,
            )
            self._add_history(smiles)

        except ValueError as e:
            messagebox.showerror(
                f"{APP_NAME} - {self._t('error_not_found_title', '未找到结果')}",
                f"{self._t('error_code', '错误代码')}: 1001\n{e}",
            )
        except Exception as e:
            messagebox.showerror(
                f"{APP_NAME} - {self._t('error_unknown_title', '未知错误')}",
                f"{self._t('error_code', '错误代码')}: 1000\n{e}",
            )

    def _on_overlay_force_field(self):
        """力场叠加显示回调"""
        try:
            smiles = self._get_selected_smiles()
            if not smiles:
                return

            img = overlay_force_field(smiles)
            tk_img = ImageTk.PhotoImage(img)
            self.result_label.config(image=tk_img)
            self.result_label.image = tk_img

        except ValueError as e:
            messagebox.showerror(
                f"{APP_NAME} - {self._t('error_not_found_title', '未找到结果')}",
                f"{self._t('error_code', '错误代码')}: 1001\n{e}",
            )
        except RuntimeError as e:
            messagebox.showerror(
                f"{APP_NAME} - {self._t('error_generation_failed_title', '生成失败')}",
                f"{self._t('error_code', '错误代码')}: 1002\n{e}",
            )
        except Exception as e:
            messagebox.showerror(
                f"{APP_NAME} - {self._t('error_unknown_title', '未知错误')}",
                f"{self._t('error_code', '错误代码')}: 1000\n{e}",
            )

    # ─── 保存功能 ─────────────────────────────────────────────

    def _save_as_png(self):
        if not self.result_label.image:
            return
        file_path = filedialog.asksaveasfilename(
            defaultextension=".png", filetypes=[("PNG files", "*.png")]
        )
        if file_path:
            self.result_label.image._PhotoImage__photo.write(file_path)
            messagebox.showinfo(
                f"{APP_NAME} - {self._t('save_success_title', '保存成功')}",
                f"{self._t('save_success_message', '图像已保存到')} {file_path}",
            )

    def _save_as_svg(self):
        try:
            smiles = self._get_selected_smiles()
            if not smiles:
                return

            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"无效的 SMILES: {smiles}")

            file_path = filedialog.asksaveasfilename(
                defaultextension=".svg", filetypes=[("SVG files", "*.svg")]
            )
            if file_path:
                Draw.MolToFile(mol, file_path, imageType="svg")
                messagebox.showinfo(
                    f"{APP_NAME} - {self._t('save_success_title', '保存成功')}",
                    f"{self._t('save_success_message', '图像已保存到')} {file_path}",
                )
        except (ValueError, RuntimeError) as e:
            messagebox.showerror(f"{APP_NAME} - Error", str(e))

    # ─── 数据库操作 ───────────────────────────────────────────

    def _change_database(self):
        file_path = filedialog.askopenfilename(filetypes=[("XML files", "*.xml")])
        if file_path:
            self._load_database(file_path)

    def _load_database(self, file_path: str):
        """在后台线程中加载数据库"""
        progress_win = Toplevel(self.root)
        progress_win.title(f"{APP_NAME} - {self._t('loading_database_title', '加载数据库中')}")
        progress_win.geometry("300x100")
        progress_win.iconbitmap("lib/media/nctl.ico")
        progress_win.attributes("-toolwindow", 2)

        tk.Label(
            progress_win,
            text=self._t("loading_database_message", "加载中，请稍候..."),
            font=("Arial", 12),
        ).pack(pady=10)
        progress_bar = ttk.Progressbar(progress_win, orient="horizontal", length=250, mode="indeterminate")
        progress_bar.pack(pady=10)
        progress_bar.start()

        def _do_load():
            try:
                new_dict = load_smiles_database(file_path)
                self.smiles_dict = new_dict
                self.database_path = file_path
                info = get_database_info(file_path)
                info_str = "\n".join(f"{k}: {v}" for k, v in info.items())
                self.root.after(0, lambda: messagebox.showinfo(
                    f"{APP_NAME} - {self._t('database_changed_title', '数据库已更换')}",
                    f"{self._t('current_database_info', '当前数据库信息')}:\n{info_str}",
                ))
            except Exception as e:
                self.root.after(0, lambda: messagebox.showerror(
                    f"{APP_NAME} - {self._t('error_unknown_title', '错误')}",
                    f"{self._t('error_loading_database', '无法加载数据库')}: {e}",
                ))
            finally:
                self.root.after(0, progress_win.destroy)

        threading.Thread(target=_do_load, daemon=True).start()

    def _show_database_info(self):
        info = get_database_info(self.database_path)
        info_str = "\n".join(f"{self._t(k, k)}: {v}" for k, v in info.items())
        messagebox.showinfo(
            f"{APP_NAME} - {self._t('database_info', '数据库信息')}", info_str
        )

    # ─── 历史记录 ─────────────────────────────────────────────

    def _update_history_menu(self):
        self.history_menu.delete(0, tk.END)
        for entry in self.history[-5:]:
            self.history_menu.add_command(
                label=entry["smiles"],
                command=lambda s=entry["smiles"]: self.formula_entry.insert(0, s),
            )

    def _show_long_history(self):
        win = Toplevel(self.root)
        win.title(self._t("history_title", "历史记录"))
        win.geometry("600x400")
        win.iconbitmap("lib/media/nctl.ico")

        columns = ("timestamp", "input_text", "smiles")
        tree = ttk.Treeview(win, columns=columns, show="headings")
        tree.heading("timestamp", text=self._t("timestamp", "时间戳"))
        tree.heading("input_text", text=self._t("input_text", "输入文本"))
        tree.heading("smiles", text="SMILES")

        for entry in self.history:
            tree.insert("", "end", values=(entry["timestamp"], entry["input_text"], entry["smiles"]))
        tree.pack(fill="both", expand=True)

        def _delete_selected():
            for item in tree.selection():
                vals = tree.item(item, "values")
                to_del = next(
                    (e for e in self.history
                     if e["timestamp"] == vals[0] and e["input_text"] == vals[1] and e["smiles"] == vals[2]),
                    None,
                )
                if to_del:
                    self.history.remove(to_del)
                    tree.delete(item)
            save_history(self.history)
            self._update_history_menu()

        def _clear_all():
            self.history = clear_history()
            for item in tree.get_children():
                tree.delete(item)
            self._update_history_menu()

        btn_frame = tk.Frame(win)
        btn_frame.pack(pady=10)
        tk.Button(btn_frame, text=self._t("delete", "删除"), command=_delete_selected).pack(side=tk.LEFT, padx=5)
        tk.Button(btn_frame, text=self._t("clear_history", "清空历史记录"), command=_clear_all).pack(side=tk.LEFT, padx=5)

    # ─── SMILES 选择弹窗 ──────────────────────────────────────

    def _show_smiles_selection(self, smiles_list: list) -> str:
        """当有多个 SMILES 匹配时，弹窗让用户选择"""
        sel_win = Toplevel(self.root)
        sel_win.title(f"{APP_NAME} - {self._t('select_smiles_title', '选择 SMILES')}")
        sel_win.geometry("600x400")
        sel_win.iconbitmap("lib/media/nctl.ico")

        selected = StringVar()

        canvas = Canvas(sel_win)
        scrollbar = Scrollbar(sel_win, orient="vertical", command=canvas.yview)
        scroll_frame = tk.Frame(canvas)
        scroll_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
        canvas.create_window((0, 0), window=scroll_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        def _select(s):
            selected.set(s)
            sel_win.destroy()

        for smiles in smiles_list:
            try:
                img = formula_to_bondline(smiles)
                img = img.convert("RGBA")
                tk_img = ImageTk.PhotoImage(img)

                frame = tk.Frame(scroll_frame, pady=10)
                frame.pack(fill="x")

                btn = tk.Button(frame, image=tk_img, text=smiles, compound=tk.TOP,
                                command=lambda s=smiles: _select(s))
                btn.image = tk_img
                btn.pack()
                tk.Label(frame, text=smiles, font=("Arial", 10)).pack()
            except Exception as e:
                tk.Label(scroll_frame, text=f"{self._t('error_loading_smiles', '无法加载')} {smiles}: {e}").pack()

        sel_win.wait_window()
        return selected.get()

    # ─── 3D 查看器 ────────────────────────────────────────────

    def _show_3d_viewer(self, smiles: str):
        win = Toplevel(self.root)
        win.title(f"{APP_NAME} - 3D Viewer")
        win.iconbitmap("lib/media/nctl.ico")
        MoleculeViewer(win, smiles, model_type=self.model3d_var.get(), lang_dict=self.lang)

    # ─── 配置窗口 ─────────────────────────────────────────────

    def _show_config(self):
        win = Toplevel(self.root)
        win.title(self._t("config_title", "配置"))
        win.geometry("400x350")
        win.iconbitmap("lib/media/nctl.ico")

        # 语言选择
        tk.Label(win, text=self._t("select_language", "选择语言："), font=("Arial", 12)).pack(pady=5)
        available_langs = self.config.get("available_languages", ["zh_cn", "en_us"])
        lang_names = [load_language(l).get("language_name", l) for l in available_langs]
        lang_var = StringVar(value=self.lang.get("language_name", self.config.get("language", "zh_cn")))
        OptionMenu(win, lang_var, *lang_names).pack(pady=5)

        # 默认数据库
        tk.Label(win, text=self._t("select_database", "选择默认数据库："), font=("Arial", 12)).pack(pady=5)
        db_files = [f for f in os.listdir("lib/db") if f.endswith(".xml")]
        db_var = StringVar(value=os.path.basename(self.database_path))
        OptionMenu(win, db_var, *db_files).pack(pady=5)

        # 历史记录开关
        record_var = tk.BooleanVar(value=self.config.get("record_history", True))
        tk.Checkbutton(win, text=self._t("record_history", "记录历史记录"), variable=record_var).pack(pady=5)

        # 3D 模型选择
        tk.Label(win, text=self._t("select_3d_model", "选择3D模型："), font=("Arial", 12)).pack(pady=5)
        model_var = StringVar(value=self.config.get("default_3d_model", "ball_and_stick"))
        OptionMenu(win, model_var, "ball_and_stick", "space_filling").pack(pady=5)

        def _save():
            # 找到选中的语言代码
            selected_lang = self.config.get("language", "zh_cn")
            for idx, name in enumerate(lang_names):
                if name == lang_var.get():
                    selected_lang = available_langs[idx]
                    break

            self.config["language"] = selected_lang
            self.config["default_database"] = f"lib/db/{db_var.get()}"
            self.config["record_history"] = record_var.get()
            self.config["default_3d_model"] = model_var.get()
            save_config(self.config)
            messagebox.showinfo(
                self._t("config_saved_title", "配置已保存"),
                self._t("config_saved_message", "配置已保存，请重启应用以应用更改。"),
            )
            win.destroy()

        tk.Button(win, text=self._t("save_button", "保存"), command=_save).pack(pady=15)

    # ─── 关于 ─────────────────────────────────────────────────

    def _show_about(self):
        win = Toplevel(self.root)
        win.title(f"{APP_NAME} - {self._t('about_developer_title', '关于开发者')}")
        win.geometry("300x400")
        win.iconbitmap("lib/media/nctl.ico")

        try:
            icon = tk.PhotoImage(file="lib/media/chem2line.png")
            icon_label = tk.Label(win, image=icon)
            icon_label.image = icon
            icon_label.pack(pady=10)
        except Exception:
            pass

        info = (
            f"{self._t('developer', '开发者')}: {DEVELOPER}\n"
            f"{self._t('version', '版本')}: {VERSION}\n"
            f"{self._t('date', '日期')}: {DATE}\n"
            f"{self._t('core_version', '内核版本')}: {CORE_VERSION}"
        )
        tk.Label(win, text=info, font=("Arial", 12), justify="left").pack(pady=10)

    def _show_repository(self):
        messagebox.showinfo(
            f"{APP_NAME} - {self._t('repository_title', '软件仓库')}",
            f"GitHub: {REPOSITORY}",
        )


# ─── 程序入口 ─────────────────────────────────────────────────

if __name__ == "__main__":
    app = Chem2LineApp()
    app.run()
