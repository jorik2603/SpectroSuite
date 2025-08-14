import tkinter as tk
from tkinter import ttk

class ModeSelectionDialog:
    def __init__(self, parent):
        self.parent = parent
        self.top = tk.Toplevel(parent)
        self.top.title("Select Mode")
        self.mode = None

        ttk.Label(self.top, text="Choose an analysis mode:", font="-weight bold").pack(pady=10)
        
        ttk.Button(self.top, text="View Mode (1 Dataset)", command=lambda: self.select_mode("view")).pack(pady=5, padx=20, fill=tk.X)
        ttk.Button(self.top, text="Compare Two Datasets", command=lambda: self.select_mode("compare_two")).pack(pady=5, padx=20, fill=tk.X)
        ttk.Button(self.top, text="Compare Multiple Datasets (Up to 4)", command=lambda: self.select_mode("compare_multi")).pack(pady=5, padx=20, fill=tk.X)
        
        self.top.update_idletasks()
        screen_width = self.top.winfo_screenwidth()
        screen_height = self.top.winfo_screenheight()
        dialog_width = self.top.winfo_width()
        dialog_height = self.top.winfo_height()
        x = (screen_width - dialog_width) // 2
        y = (screen_height - dialog_height) // 2
        self.top.geometry(f"+{x}+{y}")
        
        self.top.transient(parent)
        self.top.grab_set()
        self.parent.wait_window(self.top)

    def select_mode(self, mode):
        self.mode = mode
        self.top.destroy()