class DataLoader:
    """Base class for all data loaders."""
    def load(self, filepath):
        """Loads data from a file. Must be implemented by subclasses."""
        raise NotImplementedError