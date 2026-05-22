"""调试工具模块

通过 --debug 命令行参数启用，提供运行时变量查看和修改功能。
"""

import tkinter as tk
from tkinter import simpledialog, messagebox


class Debugger:
    """运行时调试器窗口"""

    def __init__(self, root: tk.Tk, globals_dict: dict):
        self.root = root
        self.globals_dict = globals_dict
        self._create_window()

    def _create_window(self):
        self.window = tk.Toplevel(self.root)
        self.window.title("Chem2Line Debug")
        self.window.geometry("500x350")

        # 变量列表
        self.var_listbox = tk.Listbox(self.window, font=("Consolas", 9))
        self.var_listbox.pack(fill=tk.BOTH, expand=True)
        self.var_listbox.bind("<Double-1>", self._edit_variable)

        self._refresh_vars()
        self._create_controls()

    def _create_controls(self):
        frame = tk.Frame(self.window)
        frame.pack(pady=10)

        buttons = [
            ("Refresh", self._refresh_vars),
            ("Show Vars", self._show_variables),
            ("Call Func", self._call_function),
            ("Simulate Error", self._simulate_error),
        ]
        for text, cmd in buttons:
            tk.Button(frame, text=text, command=cmd).pack(side=tk.LEFT, padx=5)

    def _refresh_vars(self):
        self.var_listbox.delete(0, tk.END)
        for name, value in sorted(self.globals_dict.items()):
            if not name.startswith("_"):
                self.var_listbox.insert(tk.END, f"{name}: {repr(value)[:80]}")

    def _edit_variable(self, event):
        selected = self.var_listbox.curselection()
        if not selected:
            return
        var_name = self.var_listbox.get(selected).split(":")[0].strip()
        new_value = simpledialog.askstring("Edit Variable", f"New value for {var_name}:")
        if new_value is not None:
            try:
                self.globals_dict[var_name] = eval(new_value)  # noqa: S307
                self._refresh_vars()
            except Exception as e:
                messagebox.showerror("Error", f"Failed to set variable: {e}")

    def _show_variables(self):
        info = "\n".join(
            f"{k}: {repr(v)[:60]}"
            for k, v in sorted(self.globals_dict.items())
            if not k.startswith("_")
        )
        messagebox.showinfo("Variables", info[:2000])

    def _call_function(self):
        func_name = simpledialog.askstring("Call Function", "Enter function name:")
        if func_name and func_name in self.globals_dict:
            try:
                result = self.globals_dict[func_name]()
                messagebox.showinfo("Result", f"Result: {result}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed: {e}")

    def _simulate_error(self):
        messagebox.showerror("Chem2Line - Error", "Error Code: 1999 (Simulated)")


def enable_debug_mode(root: tk.Tk, globals_dict: dict):
    """启用调试模式，打开调试窗口"""
    Debugger(root, globals_dict)
