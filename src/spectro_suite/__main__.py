import tkinter as tk
from spectro_suite.core import SpectroscopyGUI
from spectro_suite.dialog import ModeSelectionDialog

def main():
    """Main function to run the application."""
    root = tk.Tk()
    # Commented out because program will not load with this enabled. Fix this later
    # root.withdraw() 
    
    dialog = ModeSelectionDialog(root)
    mode = dialog.mode

    if mode:
        root.deiconify() 
        app = SpectroscopyGUI(root, mode)
        root.mainloop()
    else:
        root.destroy()

if __name__ == "__main__":
    main()