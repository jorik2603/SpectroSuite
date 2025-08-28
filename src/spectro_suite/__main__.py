# src/spectro_suite/__main__.py
import tkinter as tk
# from spectro_suite.dialog import ModeSelectionDialog
# from spectro_suite.view_mode import ViewModeGUI
# from spectro_suite.compare_two_mode import CompareTwoGUI
# from spectro_suite.compare_multi_mode import CompareMultiGUI
# from spectro_suite.fft_mode import FFTModeGUI

from dialog import ModeSelectionDialog
from view_mode import ViewModeGUI
from compare_two_mode import CompareTwoGUI
from compare_multi_mode import CompareMultiGUI
from fft_mode import FFTModeGUI
def main():
    """Main function to run the application."""
    root = tk.Tk()

    dialog = ModeSelectionDialog(root)
    mode = dialog.mode

    if mode:
        root.deiconify()

        if mode == "view":
            app = ViewModeGUI(root, mode)
        elif mode == "fft":
            app = FFTModeGUI(root, mode)
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