import numpy as np
from .base_loader import DataLoader # Changed import

class DatLoader(DataLoader):
    """Loads data from the specific .dat file format."""
    def load(self, filepath):
        with open(filepath, 'r') as f:
            x_axis = np.array([float(val) for val in f.readline().strip().split(',')[1:]])
            y_axis, data_rows = [], []
            for line in f:
                parts = line.strip().split(',')
                y_axis.append(float(parts[0]))
                data_rows.append([float(val) for val in parts[1:]])
        return x_axis, np.array(y_axis), np.array(data_rows)