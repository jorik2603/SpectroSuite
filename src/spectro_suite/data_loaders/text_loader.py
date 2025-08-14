import numpy as np
from .base_loader import DataLoader # Changed import

class TextLoader(DataLoader):
    """Loads data from a generic text file (.txt, .csv, etc.)."""
    def load(self, filepath):
        data = np.loadtxt(filepath, delimiter=',')
        x_axis = np.arange(data.shape[1])
        y_axis = np.arange(data.shape[0])
        return x_axis, y_axis, data