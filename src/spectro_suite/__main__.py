import tkinter as tk
from spectro_suite.dialog import ModeSelectionDialog
from spectro_suite.view_mode import ViewModeGUI
from spectro_suite.compare_two_mode import CompareTwoGUI
from spectro_suite.compare_multi_mode import CompareMultiGUI

def main():
    """Main function to run the application."""
    root = tk.Tk()
   
    dialog = ModeSelectionDialog(root)
    mode = dialog.mode

    if mode:
        root.deiconify()
        
        if mode == "view":
            app = ViewModeGUI(root, mode)
        elif mode == "compare_two":
            app = CompareTwoGUI(root, mode)
        elif mode == "compare_multi":
            app = CompareMultiGUI(root, mode)
        else:
            root.destroy()
            return
            
        root.mainloop()
    else:
        root.destroy()

if __name__ == "__main__":
    main()