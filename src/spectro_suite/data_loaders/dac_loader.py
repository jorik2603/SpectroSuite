import numpy as np
from .base_loader import DataLoader

class DacLoader(DataLoader):
    """
    A data loader for .dac files.
    """
    def load(self, filepath):
        """
        Loads data from a .dac file.

        Args:
            filepath (str): The path to the .dac file.

        Returns:
            tuple: A tuple containing the x-axis (wavelengths),
                   y-axis (time), and z-data (intensity).
        """
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # First line is the wavelength axis. We need to average the pairs.
        raw_wavelengths = np.array(lines[0].strip().split('\t')[1:], dtype=float)
        wavelengths = (raw_wavelengths[::2] + raw_wavelengths[1::2]) / 2.0

        # Subsequent lines are the time and intensity data
        time_data = []
        intensity_data = []
        for line in lines[1:]:
            parts = line.strip().split('\t')
            if len(parts) > 1:
                time_data.append(float(parts[0]))
                # The rest of the values are intensity data, paired up.
                # We take the average of each pair.
                raw_intensity = np.array(parts[1:], dtype=float)
                avg_intensity = (raw_intensity[::2] + raw_intensity[1::2]) / 2.0
                intensity_data.append(avg_intensity)


        x_axis = wavelengths
        y_axis = np.array(time_data)
        z_data = np.array(intensity_data)

        # Ensure the dimensions of the final data are consistent
        if z_data.shape[1] != len(x_axis):
             raise ValueError("The number of wavelength points does not match the number of intensity data points.")
        if z_data.shape[0] != len(y_axis):
            raise ValueError("The number of time points does not match the number of intensity data rows.")

        return x_axis, y_axis, z_data