import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

# Import the new data loading factory function
from data_loaders import get_loader

try:
    from scipy.signal import savgol_filter
except ImportError:
    savgol_filter = None

try:
    from cmcrameri import cm
    CRAMERI_CMAPS = ['cmc.cork', 'cmc.batlow', 'cmc.vik']
except ImportError:
    CRAMERI_CMAPS = []

class BaseSpectroscopyGUI:
    """Base class containing all shared functionality for the spectroscopy GUI."""
    def __init__(self, root, mode):
        self.root = root
        self.root.title(f"SpectroSuite - {mode.replace('_', ' ').title()} Mode")
        self.mode = mode
        try:
            self.root.state('zoomed')
        except tk.TclError:
            self.root.attributes('-fullscreen', True)

        self.datasets = {} 
        self.subplot_axes = []
        self._block_callbacks = False

        self._create_controls_panel()
        self._create_plot_canvases()
        self.initial_plot_setup()
        
        # --- NEW: Attributes to cache the state of data processing settings ---
        self._last_norm_settings = None
        self._last_bg_settings = None

    def _create_controls_panel(self):
        main_controls_frame = ttk.LabelFrame(self.root, text="Plotting Controls", padding=10)
        main_controls_frame.pack(side=tk.TOP, fill=tk.X, padx=10, pady=5)

        # --- RESTRUCTURED: Use a consistent grid layout to prevent manager conflicts ---
        # Configure the grid columns to space out the control groups
        main_controls_frame.grid_columnconfigure(0, weight=1)
        main_controls_frame.grid_columnconfigure(1, weight=1)
        main_controls_frame.grid_columnconfigure(2, weight=1)
        main_controls_frame.grid_columnconfigure(3, weight=1)
        main_controls_frame.grid_columnconfigure(4, weight=1)
        main_controls_frame.grid_columnconfigure(5, weight=1)
        main_controls_frame.grid_columnconfigure(6, weight=1)
        main_controls_frame.grid_columnconfigure(7, weight=1)
        
        # --- Group 1: Loading & Processing ---
        loading_frame = ttk.Frame(main_controls_frame)
        loading_frame.grid(row=0, column=0, sticky='ns', padx=(0, 10))
        self._create_file_loading_controls(loading_frame) # Subclasses use .pack() here, which is fine inside their own frame
              
        # --- Group 2: Processing Controls ---
        processing_frame = ttk.Frame(main_controls_frame)
        processing_frame.grid(row=0, column=1, sticky='ns', padx=10)
        bg_frame = ttk.LabelFrame(processing_frame, text="Background Subtraction")
        bg_frame.pack(side=tk.TOP, fill=tk.Y) # .pack is okay here
        self.bg_subtract_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(bg_frame, text="Enable", variable=self.bg_subtract_var).grid(row=0, column=0, columnspan=2, sticky='w', padx=5)
        ttk.Label(bg_frame, text="BG Min:").grid(row=1, column=0, sticky='w', padx=5, pady=2)
        self.bg_time_min_var = tk.StringVar(value='-1.0')
        ttk.Entry(bg_frame, textvariable=self.bg_time_min_var, width=8).grid(row=1, column=1, sticky='ew', padx=5, pady=2)
        ttk.Label(bg_frame, text="BG Max:").grid(row=2, column=0, sticky='w', padx=5, pady=2)
        self.bg_time_max_var = tk.StringVar(value='0.0')
        ttk.Entry(bg_frame, textvariable=self.bg_time_max_var, width=8).grid(row=2, column=1, sticky='ew', padx=5, pady=2)
        
        # --- Group 3: Slicing Controls ---
        slice_frame = ttk.Frame(main_controls_frame)
        slice_frame.grid(row=0, column=2, sticky='ns', padx=10)
        self._create_slice_controls(slice_frame) # Uses grid

        # --- Group 4: Window controls ---
        win_frame = ttk.Frame(main_controls_frame)
        win_frame.grid(row=0, column=3, sticky='ns', padx=10)
        self._create_window_controls(win_frame) # Uses grid

        # --- Group 5: Display controls ---
        display_frame = ttk.Frame(main_controls_frame)
        display_frame.grid(row=0, column=4, sticky='ns', padx=10)
        self._create_display_controls(display_frame) # Uses grid        
      
        # --- Group 6: Axis controls
        axis_frame = ttk.Frame(main_controls_frame)
        axis_frame.grid(row=0, column=5, sticky='ns', padx=10)
        self._create_axis_controls(axis_frame) # Uses grid        

        # --- Group 7: Time Offsets ---
        offset_frame = ttk.LabelFrame(main_controls_frame, text="Time Offsets (ps)")
        offset_frame.grid(row=0, column=6, sticky='ns', padx=10)
        self.time_offset_vars = {}
        for i in range(1, 5):
            self.time_offset_vars[i] = tk.StringVar(value='0.0')
            ttk.Label(offset_frame, text=f"D{i}:").grid(row=i-1, column=0, sticky='w')
            ttk.Entry(offset_frame, textvariable=self.time_offset_vars[i], width=6).grid(row=i-1, column=1, sticky='ew')

        # --- Group 8: Apply button ---
        ttk.Button(main_controls_frame, text="Apply\nSettings", command=self.update_plots).grid(row=0, column=7, sticky='nsew', padx=(10, 0))

    def _create_file_loading_controls(self, parent_frame):
        raise NotImplementedError("This method must be implemented by a subclass.")

    def _create_slice_controls(self, parent_frame):
        # This frame will be placed with grid, so its children can use grid
        sl_frame = ttk.Frame(parent_frame)
        sl_frame.grid(row=0, column=0, sticky='ew') # Use grid
        ttk.Label(sl_frame, text="Slice Selection", font="-weight bold").grid(row=0, column=0, columnspan=3, sticky='w')
        self.y_slice_idx, self.x_slice_idx = tk.IntVar(), tk.IntVar()
        self.y_slice_val, self.x_slice_val = tk.StringVar(value='2.00'), tk.StringVar(value='450.00')
        self.y_slice_idx.trace_add("write", self._on_slider_change)
        self.x_slice_idx.trace_add("write", self._on_slider_change)
        ttk.Label(sl_frame, text="Time (Y):").grid(row=1, column=0, sticky='w')
        y_entry = ttk.Entry(sl_frame, textvariable=self.y_slice_val, width=8)
        y_entry.grid(row=1, column=1)
        self.slider_y = ttk.Scale(sl_frame, from_=0, to=100, orient=tk.HORIZONTAL, variable=self.y_slice_idx, length=120)
        self.slider_y.grid(row=1, column=2, sticky='ew')
        ttk.Label(sl_frame, text="Wave (X):").grid(row=2, column=0, sticky='w')
        x_entry = ttk.Entry(sl_frame, textvariable=self.x_slice_val, width=8)
        x_entry.grid(row=2, column=1)
        self.slider_x = ttk.Scale(sl_frame, from_=0, to=100, orient=tk.HORIZONTAL, variable=self.x_slice_idx, length=120)
        self.slider_x.grid(row=2, column=2, sticky='ew')
        y_entry.bind("<Return>", lambda e: self._on_entry_change('y'))
        x_entry.bind("<Return>", lambda e: self._on_entry_change('x'))
    
    def _create_display_controls(self, parent_frame):
        # This frame will be placed with grid, so its children can use grid
        disp_frame = ttk.Frame(parent_frame)
        disp_frame.grid(row=0, column=0, sticky='ew') # Use grid
        ttk.Label(disp_frame, text="Display Settings", font="-weight bold").grid(row=0, column=0, columnspan=2, sticky='w')
        ttk.Label(disp_frame, text="Colormap:").grid(row=1, column=0, sticky='w')
        self.cmap_var = tk.StringVar(value='cmc.cork')
        self.cmap_combo = ttk.Combobox(disp_frame, textvariable=self.cmap_var, values=['viridis'] + CRAMERI_CMAPS, width=12)
        self.cmap_combo.grid(row=1, column=1, pady=2, sticky='w')
        self.smoothing_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(disp_frame, text="Enable Smoothing", variable=self.smoothing_var).grid(row=2, column=0, columnspan=2, sticky='w')
        ttk.Label(disp_frame, text="Window:").grid(row=3, column=0, sticky='w')
        self.smooth_window = tk.StringVar(value='5')
        ttk.Entry(disp_frame, textvariable=self.smooth_window, width=5).grid(row=3, column=1, sticky='w')
        ttk.Label(disp_frame, text="Polyorder:").grid(row=4, column=0, sticky='w')
        self.smooth_poly = tk.StringVar(value='2')
        ttk.Entry(disp_frame, textvariable=self.smooth_poly, width=5).grid(row=4, column=1, sticky='w')


    def _create_axis_controls(self, parent_frame):
        # This frame will be placed with grid, so its children can use grid
        ax_frame = ttk.Frame(parent_frame)
        ax_frame.grid(row=2, column=0, sticky='ew') # Use grid
        ttk.Label(ax_frame, text="Axis Scaling", font="-weight bold").grid(row=0, column=0, columnspan=2, sticky='w')
        self.symlog_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(ax_frame, text="Symlog Time Axis", variable=self.symlog_var).grid(row=1, column=0, columnspan=2, sticky='w')
        ttk.Label(ax_frame, text="Linthresh:").grid(row=2, column=0, sticky='w')
        self.linthresh_var = tk.StringVar(value='1.0')
        ttk.Entry(ax_frame, textvariable=self.linthresh_var, width=10).grid(row=2, column=1, sticky='w')
        ttk.Label(ax_frame, text="Intensity Min:").grid(row=3, column=0, sticky='w')
        self.slice_ymin_var = tk.StringVar(value='-15')
        ttk.Entry(ax_frame, textvariable=self.slice_ymin_var, width=10).grid(row=3, column=1, sticky='w')
        ttk.Label(ax_frame, text="Intensity Max:").grid(row=4, column=0, sticky='w')
        self.slice_ymax_var = tk.StringVar(value='15')
        ttk.Entry(ax_frame, textvariable=self.slice_ymax_var, width=10).grid(row=4, column=1, sticky='w')
        self.sync_color_var = tk.BooleanVar(value=True)
        ttk.Checkbutton(ax_frame, text="Sync Color & Intensity", variable=self.sync_color_var).grid(row=5, column=0, columnspan=2, sticky='w')


    def _create_window_controls(self, parent_frame):
        # This frame will be placed with grid, so its children can use grid
        win_frame = ttk.Frame(parent_frame)
        win_frame.grid(row=2, column=0, sticky='ew') # Use grid
        ttk.Label(win_frame, text="Data Window", font="-weight bold").grid(row=0, column=0, columnspan=2, sticky='w')
        ttk.Label(win_frame, text="Time Min:").grid(row=1, column=0, sticky='w')
        self.time_min_var = tk.StringVar(value='0')
        ttk.Entry(win_frame, textvariable=self.time_min_var, width=8).grid(row=1, column=1)
        ttk.Label(win_frame, text="Time Max:").grid(row=2, column=0, sticky='w')
        self.time_max_var = tk.StringVar(value='20')
        ttk.Entry(win_frame, textvariable=self.time_max_var, width=8).grid(row=2, column=1)
        ttk.Label(win_frame, text="Wave Min:").grid(row=3, column=0, sticky='w')
        self.wave_min_var = tk.StringVar(value='365')
        ttk.Entry(win_frame, textvariable=self.wave_min_var, width=8).grid(row=3, column=1)
        ttk.Label(win_frame, text="Wave Max:").grid(row=4, column=0, sticky='w')
        self.wave_max_var = tk.StringVar(value='812')
        ttk.Entry(win_frame, textvariable=self.wave_max_var, width=8).grid(row=4, column=1)

    def _create_plot_canvases(self):
        canvas_frame = ttk.Frame(self.root, padding=(10,0))
        canvas_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        # --- MODIFIED: Adjusted figure sizes for better layout ---
        self.fig_2d = plt.Figure(figsize=(7, 4)) # Taller 2D plot
        self.canvas_2d = FigureCanvasTkAgg(self.fig_2d, master=canvas_frame)
        NavigationToolbar2Tk(self.canvas_2d, canvas_frame)
        self.canvas_2d.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        self.fig_slices = plt.Figure(figsize=(10, 5)) # Taller slice plots
        self.ax_slice_y, self.ax_slice_x = self.fig_slices.subplots(1, 2)
        self.canvas_slices = FigureCanvasTkAgg(self.fig_slices, master=canvas_frame)
        NavigationToolbar2Tk(self.canvas_slices, canvas_frame)
        self.canvas_slices.get_tk_widget().pack(fill=tk.BOTH, expand=True, pady=(10, 0))

    def _processing_settings_changed(self):
        """Check if data processing settings have changed since the last update."""
        current_norm_settings = self.norm_type_var.get()
        current_bg_settings = (
            self.bg_subtract_var.get(),
            self.bg_time_min_var.get(),
            self.bg_time_max_var.get()
        )

        if current_norm_settings != self._last_norm_settings or \
           current_bg_settings != self._last_bg_settings:
            return True
        return False
    
    def load_data(self, dataset_num):
        filepath = filedialog.askopenfilename(
            title="Select a data file",
            filetypes=(
                ("DAT files", "*.dat"),
                ("NPY files", "*.npy"),
                ("Text files", "*.txt"),
                ("CSV files", "*.csv"),
                ("DAC files", "*.dac"),
                ("All files", "*.*")
            )
        )
        if not filepath:
            return
        try:
            from data_loaders import get_loader
            loader = get_loader(filepath)
            x_axis, y_axis, data = loader.load(filepath)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load file: {e}"); return
        
        # Store original data and a working copy separately
        self.datasets[dataset_num] = {
            'data': np.copy(data), # The working copy that gets processed and plotted
            'original_data': data, # The pristine, untouched data
            'x_axis': x_axis,
            'y_axis': y_axis
        }
        
        if len(self.datasets) == 1:
            self.slider_y.config(to=data.shape[0] - 1)
            self.slider_x.config(to=data.shape[1] - 1)
            self._on_entry_change('y'); self._on_entry_change('x')
        self.update_plots()

    def _apply_background_subtraction(self, data, y_axis):
        """Applies background subtraction to the data based on UI settings."""
        if self.bg_subtract_var.get():
            try:
                bg_min = float(self.bg_time_min_var.get())
                bg_max = float(self.bg_time_max_var.get())
                if bg_min < bg_max:
                    bg_indices = np.where((y_axis >= bg_min) & (y_axis <= bg_max))[0]
                    if len(bg_indices) > 0:
                        background = data[bg_indices, :].mean(axis=0)
                        return data - background
            except (ValueError, KeyError):
                pass # Silently ignore errors from invalid UI input
        return data

    def _apply_normalization(self, data):
        """Applies normalization to the data based on UI settings."""
        norm_type = self.norm_type_var.get()
        if norm_type == "None":
            return data
        
        # Make a copy to avoid modifying the original data array in place
        data_to_norm = np.copy(data)

        if norm_type == "Normalize to Max":
            max_val = data_to_norm.max()
            if max_val != 0:
                return data_to_norm / max_val
        else:
            min_val, max_val = data_to_norm.min(), data_to_norm.max()
            range_val = max_val - min_val
            if range_val > 0:
                if norm_type == "Full Range [-1, 1]":
                    return 2 * (data_to_norm - min_val) / range_val - 1
                elif norm_type == "Full Range [0, 1]":
                    return (data_to_norm - min_val) / range_val
        return data # Return original data if normalization failed
    
    def initial_plot_setup(self):
        raise NotImplementedError("This method must be implemented by a subclass.")

    def _on_slider_change(self, *args):
        if 1 in self.datasets:
            dset1 = self.datasets[1]
            try:
                offset1 = float(self.time_offset_vars[1].get())
            except (ValueError, KeyError):
                offset1 = 0.0
            if self.y_slice_idx.get() < len(dset1['y_axis']): self.y_slice_val.set(f"{dset1['y_axis'][self.y_slice_idx.get()] + offset1:.2f}")
            if self.x_slice_idx.get() < len(dset1['x_axis']): self.x_slice_val.set(f"{dset1['x_axis'][self.x_slice_idx.get()]:.2f}")
        self._update_slice_plots(); self._draw_slice_lines()
        self.canvas_2d.draw_idle(); self.canvas_slices.draw_idle()

    def _on_entry_change(self, axis):
        if 1 not in self.datasets: return
        dset1 = self.datasets[1]
        try:
            offset1 = float(self.time_offset_vars[1].get())
        except (ValueError, KeyError):
            offset1 = 0.0
        
        if axis == 'y':
            try: self.y_slice_idx.set((np.abs(dset1['y_axis'] - (float(self.y_slice_val.get()) - offset1))).argmin())
            except (ValueError, IndexError): pass
        elif axis == 'x':
            try: self.x_slice_idx.set((np.abs(dset1['x_axis'] - float(self.x_slice_val.get()))).argmin())
            except (ValueError, IndexError): pass

    def update_plots(self, *args):
        if not self.datasets:
            self.initial_plot_setup(); return

        # --- MODIFIED: Only re-process data if settings have changed ---
        if self._processing_settings_changed():
            # Apply data processing dynamically on each update
            for key in self.datasets:
                dset = self.datasets[key]
                # Start with a fresh copy of the original data
                processed_data = np.copy(dset['original_data'])
                
                # Chain the processing functions
                processed_data = self._apply_background_subtraction(processed_data, dset['y_axis'])
                processed_data = self._apply_normalization(processed_data)
                
                # Store the newly processed data for the plotting functions to use
                dset['data'] = processed_data
            
            # --- NEW: Update the cached settings after processing ---
            self._last_norm_settings = self.norm_type_var.get()
            self._last_bg_settings = (
                self.bg_subtract_var.get(),
                self.bg_time_min_var.get(),
                self.bg_time_max_var.get()
            )
        
        # Clear the entire figure and recreate subplots on every update --- To prevent colorbar from duplicating
        self.fig_2d.clear()
        share_kwargs = {'sharex': True, 'sharey': True}
        if self.mode == 'view' or self.mode == 'fft': # Added fft mode here
            self.subplot_axes = [self.fig_2d.add_subplot(1, 1, 1)]
        elif self.mode == 'compare_two':
            self.subplot_axes = self.fig_2d.subplots(1, 2, **share_kwargs).flatten()
        elif self.mode == 'compare_multi':
            self.subplot_axes = self.fig_2d.subplots(2, 2, **share_kwargs).flatten()

        vmin, vmax = (float(self.slice_ymin_var.get()), float(self.slice_ymax_var.get())) if self.sync_color_var.get() and self.slice_ymin_var.get() and self.slice_ymax_var.get() else (None, None)
        cmap = self.cmap_var.get()
        plotted_keys = sorted(self.datasets.keys())
        
        for i, ax in enumerate(self.subplot_axes):
            if i < len(plotted_keys):
                key = plotted_keys[i]
                dset = self.datasets[key]
                self.plot_pcolormesh(ax, dset['x_axis'], dset['y_axis'], dset['data'], f"Dataset {key}", cmap, vmin, vmax)
            else:
                ax.clear(); ax.set_title(f"Load Dataset {i+1}")
                ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Time (ps)")
        self._draw_slice_lines()
        self.fig_2d.tight_layout(); self.canvas_2d.draw()
        self._update_slice_plots(); self.canvas_slices.draw()

    def _draw_slice_lines(self):
        if len(self.subplot_axes) == 0 or 1 not in self.datasets: return
        y_idx, x_idx = self.y_slice_idx.get(), self.x_slice_idx.get()
        dset1 = self.datasets[1]
        try:
            offset1 = float(self.time_offset_vars[1].get())
        except (ValueError, KeyError):
            offset1 = 0.0
        if y_idx >= len(dset1['y_axis']) or x_idx >= len(dset1['x_axis']): return
        y_val, x_val = dset1['y_axis'][y_idx] + offset1, dset1['x_axis'][x_idx]
        for ax in self.subplot_axes:
            for line in ax.lines: line.remove()
            ax.axhline(y_val, color='r', linestyle='--', lw=1.5)
            ax.axvline(x_val, color='cyan', linestyle='--', lw=1.5)

    def _update_slice_plots(self):
        self.ax_slice_y.clear(); self.ax_slice_x.clear()
        y_idx, x_idx = self.y_slice_idx.get(), self.x_slice_idx.get()
        
        data_slices = {}
        for key, dset in self.datasets.items():
            if y_idx < dset['data'].shape[0]: data_slices[f'{key}y_raw'] = dset['data'][y_idx, :]
            if x_idx < dset['data'].shape[1]: data_slices[f'{key}x_raw'] = dset['data'][:, x_idx]

        if self.smoothing_var.get() and savgol_filter:
            try:
                win, poly = int(self.smooth_window.get()), int(self.smooth_poly.get())
                if win % 2 == 0: win += 1
                for key in list(data_slices.keys()):
                    data = data_slices[key]
                    if data is not None and len(data) > win:
                        data_slices[key.replace('_raw', '_smooth')] = savgol_filter(data, win, poly)
            except Exception as e: print(f"Smoothing Error: {e}")

        colors = ['blue', 'orangered', 'green', 'purple']; styles = ['-', '--', '-.', ':']
        for i, (key, dset) in enumerate(sorted(self.datasets.items())):
            color, style = colors[i % len(colors)], styles[i % len(styles)]
            try:
                offset = float(self.time_offset_vars[key].get())
            except (ValueError, KeyError):
                offset = 0.0
            for axis_type, plot_axis_data in [('y', dset['x_axis']), ('x', dset['y_axis'] + offset)]:
                raw_key, smooth_key = f"{key}{axis_type}_raw", f"{key}{axis_type}_smooth"
                target_ax = self.ax_slice_y if axis_type == 'y' else self.ax_slice_x
                if raw_key in data_slices and data_slices[raw_key] is not None:
                    if smooth_key in data_slices:
                        target_ax.scatter(plot_axis_data, data_slices[raw_key], color=color, alpha=0.3, s=20)
                        target_ax.plot(plot_axis_data, data_slices[smooth_key], label=f"Data {key} (Smoothed)", c=color, ls=style)
                    else:
                        target_ax.plot(plot_axis_data, data_slices[raw_key], label=f"Data {key}", c=color, ls=style)

        self.ax_slice_y.set_title(f"Slice at Time = {self.y_slice_val.get()} ps"); self.ax_slice_y.set_xlabel("Wavelength (nm)")
        self.ax_slice_x.set_title(f"Slice at Wavelength = {self.x_slice_val.get()} nm"); self.ax_slice_x.set_xlabel("Time (ps)")
        for ax, ax_name in [(self.ax_slice_y, 'intensity'), (self.ax_slice_x, 'time')]:
            ax.set_ylabel("Intensity"); ax.legend(); ax.grid(True, linestyle=':')
            if self.symlog_var.get() and ax_name == 'time':
                try: ax.set_xscale('symlog', linthresh=float(self.linthresh_var.get()))
                except (ValueError, TypeError): ax.set_xscale('symlog')
        try:
            ymin = float(self.slice_ymin_var.get()) if self.slice_ymin_var.get() else None
            ymax = float(self.slice_ymax_var.get()) if self.slice_ymax_var.get() else None
            self.ax_slice_y.set_ylim(ymin, ymax); self.ax_slice_x.set_ylim(ymin, ymax)
            wave_min = float(self.wave_min_var.get()) if self.wave_min_var.get() else None
            wave_max = float(self.wave_max_var.get()) if self.wave_max_var.get() else None
            self.ax_slice_y.set_xlim(wave_min, wave_max)
            time_min = float(self.time_min_var.get()) if self.time_min_var.get() else None
            time_max = float(self.time_max_var.get()) if self.time_max_var.get() else None
            self.ax_slice_x.set_xlim(time_min, time_max)
        except ValueError: pass
        self.fig_slices.tight_layout()

    def plot_pcolormesh(self, ax, x, y, data, title, cmap, vmin, vmax):
        # ax.clear() is now handled by the main update_plots function
        try:
            dset_num = int(title.split()[-1])
            offset = float(self.time_offset_vars[dset_num].get())
            y_offsetted = y + offset
        except (ValueError, KeyError, IndexError):
            y_offsetted = y
        
        try:
            time_min = float(self.time_min_var.get()) if self.time_min_var.get() else y_offsetted.min()
            time_max = float(self.time_max_var.get()) if self.time_max_var.get() else y_offsetted.max()
            wave_min = float(self.wave_min_var.get()) if self.wave_min_var.get() else x.min()
            wave_max = float(self.wave_max_var.get()) if self.wave_max_var.get() else x.max()
            y_indices = np.where((y_offsetted >= time_min) & (y_offsetted <= time_max))[0]
            x_indices = np.where((x >= wave_min) & (x <= wave_max))[0]
            if len(y_indices) > 0 and len(x_indices) > 0:
                y_sliced, x_sliced = y_offsetted[y_indices], x[x_indices]
                data_sliced = data[y_indices[:, None], x_indices]
            else: x_sliced, y_sliced, data_sliced = x, y_offsetted, data
        except (ValueError, IndexError):
            x_sliced, y_sliced, data_sliced = x, y_offsetted, data

        im = ax.pcolormesh(x_sliced, y_sliced, data_sliced, cmap=cmap, shading='auto', vmin=vmin, vmax=vmax)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel("Wavelength (nm)"); ax.set_ylabel("Time (ps)")
        if len(x_sliced) > 1 and len(y_sliced) > 1:
            ax.set_xlim(x_sliced.min(), x_sliced.max())
            ax.set_ylim(y_sliced.min(), y_sliced.max())
        self.fig_2d.colorbar(im, ax=ax)