#from .base_loader import DataLoader
from .dat_loader import DatLoader
from .npy_loader import NpyLoader
from .text_loader import TextLoader

def get_loader(filepath):
    """Factory function to get the appropriate loader for a file type."""
    if filepath.endswith('.dat'):
        return DatLoader()
    elif filepath.endswith('.npy'):
        return NpyLoader()
    elif filepath.endswith(('.txt', '.csv')):
        return TextLoader()
    else:
        # Fallback to text loader for unknown types
        return TextLoader()