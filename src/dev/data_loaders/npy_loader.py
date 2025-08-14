import numpy as np
from .base_loader import DataLoader # Changed import

class NpyLoader(DataLoader):
    """Loads data from a .npy file."""
    def load(self, filepath):
        data = np.load(filepath)
        x_axis = np.arange(data.shape[1])
        y_axis = np.arange(data.shape[0])
        return x_axis, y_axis, data