import tkinter as tk
from tkinter import ttk
# from spectro_suite.base_gui import BaseSpectroscopyGUI
from base_gui import BaseSpectroscopyGUI

class ViewModeGUI(BaseSpectroscopyGUI):
    def __init__(self, root, mode):
        super().__init__(root, mode)

    def _create_file_loading_controls(self, parent_frame):
        ttk.Label(parent_frame, text="Data Loading", font="-weight bold").pack(anchor='w')
        ttk.Button(parent_frame, text="Load Dataset", command=lambda: self.load_data(1)).pack(pady=2, fill=tk.X)
        self._add_common_norm_controls(parent_frame)

    def initial_plot_setup(self):
        self.fig_2d.clear()
        self.subplot_axes = [self.fig_2d.add_subplot(1, 1, 1)]
        self.subplot_axes[0].set_title("Load a Dataset")
        self.subplot_axes[0].set_xlabel("Wavelength (nm)"); self.subplot_axes[0].set_ylabel("Time (ps)")
        self.fig_2d.tight_layout()
        self.canvas_2d.draw()

    def _add_common_norm_controls(self, parent_frame):
        ttk.Label(parent_frame, text="Normalization:", font="-weight bold").pack(anchor='w', pady=(5,0))
        self.norm_type_var = tk.StringVar(value="None")
        norm_options = ["None", "Normalize to Max", "Full Range [-1, 1]", "Full Range [0, 1]"]
        norm_combo = ttk.Combobox(parent_frame, textvariable=self.norm_type_var, values=norm_options, state="readonly")
        norm_combo.pack(pady=2, fill=tk.X)
        ttk.Label(parent_frame, text="(Set before loading data)", style="TLabel.small.TLabel").pack(anchor='w')
        self.root.style = ttk.Style()
        self.root.style.configure("TLabel.small.TLabel", font=("TkDefaultFont", 7))