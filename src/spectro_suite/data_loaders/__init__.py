from .base_loader import DataLoader
from .dat_loader import DatLoader
from .npy_loader import NpyLoader
from .text_loader import TextLoader
from .dac_loader import DacLoader

def get_loader(filepath):
    """
    Factory function to get the appropriate data loader based on the file extension.
    """
    if filepath.endswith('.dat'):
        return DatLoader()
    elif filepath.endswith('.npy'):
        return NpyLoader()
    elif filepath.endswith(('.txt', '.csv')):
        return TextLoader()
    elif filepath.endswith('.dac'):
        return DacLoader()
    else:
        raise ValueError(f"Unsupported file type: {filepath}")