# src/spectro_suite/fft_mode.py
import tkinter as tk
from tkinter import ttk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from scipy.fft import fft, fftfreq

# from spectro_suite.base_gui import BaseSpectroscopyGUI
from base_gui import BaseSpectroscopyGUI

class FFTModeGUI(BaseSpectroscopyGUI):
    def __init__(self, root, mode):
        super().__init__(root, mode)
        # No need to define fft time range variables here, they are created in _create_window_controls

    def _create_window_controls(self, parent_frame):
        # First, call the base class method to create the standard window controls
        super()._create_window_controls(parent_frame)

        # Now, add the FFT specific controls to the parent_frame of window controls
        fft_frame = ttk.LabelFrame(parent_frame, text="FFT Time Range (ps)")
        # Add some vertical padding to separate from the offset controls
        fft_frame.grid(row=6, column=0, columnspan=2, pady=(10,0), sticky='ew')

        ttk.Label(fft_frame, text="FFT Min:").grid(row=0, column=0, sticky='w', padx=5, pady=2)
        self.fft_time_min_var = tk.StringVar(value='0.0')
        ttk.Entry(fft_frame, textvariable=self.fft_time_min_var, width=8).grid(row=0, column=1, sticky='ew', padx=5, pady=2)

        ttk.Label(fft_frame, text="FFT Max:").grid(row=1, column=0, sticky='w', padx=5, pady=2)
        self.fft_time_max_var = tk.StringVar(value='1.0')
        ttk.Entry(fft_frame, textvariable=self.fft_time_max_var, width=8).grid(row=1, column=1, sticky='ew', padx=5, pady=2)

    def _create_display_controls(self, parent_frame):
        # Add the FFT unit selection to the display controls
        super()._create_display_controls(parent_frame)
        self.fft_unit_var = tk.BooleanVar(value=False)
        self.blackman_var = tk.BooleanVar(value=False)
        ttk.Checkbutton(parent_frame, text="FFT in cm⁻¹", variable=self.fft_unit_var).grid(row=5, column=0, columnspan=2, sticky='w')
        ttk.Checkbutton(parent_frame, text="Enable Blackman Window", variable=self.blackman_var).grid(row=6, column=0, columnspan=2, sticky='w')



    def _create_file_loading_controls(self, parent_frame):
        # This is identical to ViewModeGUI
        ttk.Label(parent_frame, text="Data Loading", font="-weight bold").pack(anchor='w')
        ttk.Button(parent_frame, text="Load Dataset", command=lambda: self.load_data(1)).pack(pady=2, fill=tk.X)
        self._add_common_norm_controls(parent_frame)

    def _create_plot_canvases(self):
        # Override to create three subplots for slices and FFT
        canvas_frame = ttk.Frame(self.root, padding=(10,0))
        canvas_frame.pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        self.fig_2d = plt.Figure(figsize=(10, 5))
        self.canvas_2d = FigureCanvasTkAgg(self.fig_2d, master=canvas_frame)
        NavigationToolbar2Tk(self.canvas_2d, canvas_frame)
        self.canvas_2d.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        self.fig_slices = plt.Figure(figsize=(15, 4)) # Wider figure for 3 plots
        # Create three subplots. The third one is for the FFT.
        self.ax_slice_y, self.ax_slice_x, self.ax_fft = self.fig_slices.subplots(1, 3)
        self.canvas_slices = FigureCanvasTkAgg(self.fig_slices, master=canvas_frame)
        NavigationToolbar2Tk(self.canvas_slices, canvas_frame)
        self.canvas_slices.get_tk_widget().pack(fill=tk.BOTH, expand=True, pady=(10, 0))

    def initial_plot_setup(self):
        # Setup for the main 2D plot
        self.fig_2d.clear()
        self.subplot_axes = [self.fig_2d.add_subplot(1, 1, 1)]
        self.subplot_axes[0].set_title("Load a Dataset")
        self.subplot_axes[0].set_xlabel("Wavelength (nm)")
        self.subplot_axes[0].set_ylabel("Time (ps)")
        self.fig_2d.tight_layout()
        self.canvas_2d.draw()

        # Setup for the slice and FFT plots
        self.ax_slice_y.clear()
        self.ax_slice_x.clear()
        self.ax_fft.clear()
        self.ax_slice_y.set_title("Time Slice")
        self.ax_slice_y.set_xlabel("Wavelength (nm)")
        self.ax_slice_y.set_ylabel("Intensity")
        self.ax_slice_x.set_title("Wavelength Slice")
        self.ax_slice_x.set_xlabel("Time (ps)")
        self.ax_slice_x.set_ylabel("Intensity")
        self.ax_fft.set_title("FFT")
        self.ax_fft.set_xlabel("Frequency (1/ps)")
        self.ax_fft.set_ylabel("Amplitude")
        self.fig_slices.tight_layout()
        self.canvas_slices.draw()

    def _add_common_norm_controls(self, parent_frame):
        # This is identical to ViewModeGUI
        ttk.Label(parent_frame, text="Normalization:", font="-weight bold").pack(anchor='w', pady=(5,0))
        self.norm_type_var = tk.StringVar(value="None")
        norm_options = ["None", "Normalize to Max", "Full Range [-1, 1]", "Full Range [0, 1]"]
        norm_combo = ttk.Combobox(parent_frame, textvariable=self.norm_type_var, values=norm_options, state="readonly")
        norm_combo.pack(pady=2, fill=tk.X)
        ttk.Label(parent_frame, text="(Set before loading data)", style="TLabel.small.TLabel").pack(anchor='w')
        self.root.style = ttk.Style()
        self.root.style.configure("TLabel.small.TLabel", font=("TkDefaultFont", 7))

    def update_plots(self, *args):
        if not self.datasets:
            self.initial_plot_setup()
            return

        # --- FIX for 2D plot and colorbar ---
        # Clear the entire figure to remove old colorbars before redrawing
        self.fig_2d.clear()
        # Re-create the subplot, as clearing the figure removes it
        self.subplot_axes = [self.fig_2d.add_subplot(1, 1, 1)]

        vmin, vmax = (float(self.slice_ymin_var.get()), float(self.slice_ymax_var.get())) if self.sync_color_var.get() and self.slice_ymin_var.get() and self.slice_ymax_var.get() else (None, None)
        cmap = self.cmap_var.get()

        dset = self.datasets[1]
        # The plot_pcolormesh function will now add a colorbar to a clean figure
        self.plot_pcolormesh(self.subplot_axes[0], dset['x_axis'], dset['y_axis'], dset['data'], "Dataset 1", cmap, vmin, vmax)

        self._draw_slice_lines()
        self.fig_2d.tight_layout()
        self.canvas_2d.draw()
        self._update_slice_plots()
        self.canvas_slices.draw()

    def _update_slice_plots(self):
        # Call the base class method to handle the regular slice plots (y and x slices)
        # This will clear and redraw ax_slice_y and ax_slice_x
        super()._update_slice_plots()

        # Now, add the FFT plot logic
        self.ax_fft.clear()
        self.ax_fft.set_title("FFT")
        # Set initial axis labels, will be updated based on unit selection
        self.ax_fft.set_xlabel("Frequency (1/ps)")
        self.ax_fft.set_ylabel("Amplitude")
        self.ax_fft.grid(True, linestyle=':')

        # Check if a dataset is loaded
        if 1 not in self.datasets:
            self.ax_fft.text(0.5, 0.5, "No data loaded", ha='center', va='center')
            self.canvas_slices.draw()
            return

        try:
            # Get the currently selected wavelength slice
            x_idx = self.x_slice_idx.get()
            dset = self.datasets[1]
            time_axis = dset['y_axis']
            wavelength_slice = dset['data'][:, x_idx]

            # Get FFT time range from the UI and convert to float
            fft_min = float(self.fft_time_min_var.get())
            fft_max = float(self.fft_time_max_var.get())

            # --- Error Handling and Edge Cases ---
            if fft_min >= fft_max:
                self.ax_fft.text(0.5, 0.5, "FFT Min must be less than Max", ha='center', va='center')
                self.canvas_slices.draw()
                return

            # Find the indices of the time axis that fall within the FFT range
            time_indices = np.where((time_axis >= fft_min) & (time_axis <= fft_max))[0]

            if len(time_indices) < 2:
                self.ax_fft.text(0.5, 0.5, "Not enough data points for FFT\nin the selected time range.", ha='center', va='center')
                self.canvas_slices.draw()
                return

            # Select the relevant portions of the time axis and the data slice
            time_data = time_axis[time_indices]
            slice_data = wavelength_slice[time_indices]

            # --- Perform FFT ---
            N = len(slice_data)
            if N > 1:
                # Apply Blackman window if enabled.
                if self.blackman_var.get():
                    slice_data = slice_data * np.blackman(N)
                
                # Assuming the time data is evenly spaced, calculate the sampling interval.
                sampling_interval = np.mean(np.diff(time_data))
                if sampling_interval <= 0:
                     self.ax_fft.text(0.5, 0.5, "Time data is not increasing.", ha='center', va='center')
                     self.canvas_slices.draw()
                     return

                # Compute the FFT
                yf = fft(slice_data)
                # Compute the frequency axis.
                xf = fftfreq(N, sampling_interval)[:N//2]
                
                # Unit conversion for FFT x-axis.
                if self.fft_unit_var.get():
                    # Conversion from 1/ps to cm-1
                    # 1/ps = 1e12 Hz. c = 2.9979e10 cm/s.
                    # wavenumber (cm-1) = freq (Hz) / c (cm/s)
                    # wavenumber (cm-1) = (freq (1/ps) * 1e12) / 2.9979e10
                    # wavenumber (cm-1) = freq (1/ps) * 33.3564
                    xf = xf * 33.3564
                    self.ax_fft.set_xlabel("Wavenumber (cm⁻¹)")
                else:
                    self.ax_fft.set_xlabel("Frequency (1/ps)")


                # Plot the amplitude spectrum.
                self.ax_fft.plot(xf, 2.0/N * np.abs(yf[0:N//2]))

        except (ValueError, IndexError) as e:
            # Handle cases where data is not yet valid or UI values are not valid numbers
            self.ax_fft.text(0.5, 0.5, f"FFT Error:\nCheck input values.", ha='center', va='center', wrap=True)
        except Exception as e:
            # Catch any other unexpected errors
            self.ax_fft.text(0.5, 0.5, f"An unexpected error occurred:\n{e}", ha='center', va='center', wrap=True)

        # Redraw the canvas to show the new plots
        self.fig_slices.tight_layout()
        self.canvas_slices.draw_idle()