from pathlib import Path
from platform import system
from typing import Union
from tkinter import PhotoImage, Tk, Toplevel


DIALOG_WINDOW_OFFSET_X = 25
DIALOG_WINDOW_OFFSET_Y = 25


def set_dialog_geometry(dialog: Toplevel, parent: Tk, x: int = DIALOG_WINDOW_OFFSET_X, y: int = DIALOG_WINDOW_OFFSET_Y):
    # A small worker function that will set a Toplevel dialog in front of a Tk window with a specified offset
    dialog.update_idletasks()  # do a quick pass to make sure the window has handled GUI events
    dialog.geometry(
        "%dx%d+%d+%d" % (
            dialog.winfo_width(),
            dialog.winfo_height(),
            parent.winfo_x() + x,
            parent.winfo_y() + y
        )
    )


def set_frame_or_top_level_icon(window: Union[Tk, Toplevel], icons_folder: Path) -> Path:
    """Sets the icon for the window and returns a Pathlib.Path object to the icon"""
    if system() == 'Darwin':
        icon_path = icons_folder / 'icon.icns'
        if icon_path.exists():
            window.iconbitmap(str(icon_path))
        else:
            print(f"Could not set icon for Mac, expecting to find it at {icon_path}")
    elif system() == 'Windows':
        icon_path = icons_folder / 'icon.png'
        img = PhotoImage(file=str(icon_path))
        if icon_path.exists():
            window.iconphoto(False, img)
        else:
            print(f"Could not set icon for Windows, expecting to find it at {icon_path}")
    else:  # Linux
        icon_path = icons_folder / 'icon.png'
        img = PhotoImage(file=str(icon_path))
        if icon_path.exists():
            window.iconphoto(False, img)
        else:
            print(f"Could not set icon for Windows, expecting to find it at {icon_path}")
    return icon_path
