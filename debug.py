import tkinter as tk
from tkinter import simpledialog, messagebox
import threading
import time

class Debugger:
    def __init__(self, root, globals_dict):
        self.root = root
        self.globals_dict = globals_dict
        self.create_debug_window()

    def create_debug_window(self):
        self.debug_window = tk.Toplevel(self.root)
        self.debug_window.title("Debug Window")
        self.debug_window.geometry("400x300")

        self.var_listbox = tk.Listbox(self.debug_window)
        self.var_listbox.pack(fill=tk.BOTH, expand=True)

        self.update_var_list()

        self.var_listbox.bind("<Double-1>", self.edit_variable)

        # 添加调试功能按钮
        self.create_debug_controls()

    def update_var_list(self):
        self.var_listbox.delete(0, tk.END)
        for var_name, var_value in self.globals_dict.items():
            self.var_listbox.insert(tk.END, f"{var_name}: {var_value}")

    def edit_variable(self, event):
        selected = self.var_listbox.curselection()
        if selected:
            var_name = self.var_listbox.get(selected).split(":")[0]
            new_value = simpledialog.askstring("Edit Variable", f"New value for {var_name}:")
            if new_value is not None:
                try:
                    self.globals_dict[var_name] = eval(new_value)
                    self.update_var_list()
                except Exception as e:
                    messagebox.showerror("Error", f"Failed to set variable: {e}")

    def create_debug_controls(self):
        control_frame = tk.Frame(self.debug_window)
        control_frame.pack(pady=10)

        tk.Button(control_frame, text="Dead", command=self.freeze_program).pack(side=tk.LEFT, padx=5)
        tk.Button(control_frame, text="Error", command=self.simulate_error).pack(side=tk.LEFT, padx=5)
        tk.Button(control_frame, text="Show variables", command=self.show_variables).pack(side=tk.LEFT, padx=5)
        tk.Button(control_frame, text="Call Function", command=self.call_function).pack(side=tk.LEFT, padx=5)

    def freeze_program(self):
        self.root.after(100, self.root.quit)
        while True:
            pass

    def simulate_error(self):
        try:
            raise RuntimeError("1999")
        except RuntimeError as e:
            messagebox.showerror("Chem2Line - Error", f"Error Code: {e}")

    def show_variables(self):
        var_info = "\n".join([f"{var_name}: {var_value}" for var_name, var_value in self.globals_dict.items()])
        messagebox.showinfo("Variables", var_info)

    def call_function(self):
        func_name = simpledialog.askstring("Call Function", "Enter function name:")
        if func_name and func_name in self.globals_dict:
            try:
                result = self.globals_dict[func_name]()
                messagebox.showinfo("Function Result", f"Result: {result}")
            except Exception as e:
                messagebox.showerror("Error", f"Failed to call function: {e}")

def enable_debug_mode(root, globals_dict):
    Debugger(root, globals_dict)
