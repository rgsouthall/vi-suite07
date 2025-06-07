from pathlib import Path
from typing import Callable, Dict
from tkinter import Tk, Label, Button, Toplevel, TOP, BOTH, LEFT

from eplaunch.interface import set_dialog_geometry, set_frame_or_top_level_icon


class TkGenericDialog:
    def __init__(self, default_padding: Dict[str, int]):
        self.pad = default_padding

    def display(self, parent_window: Tk, title: str, text: str):
        t = Toplevel(parent_window)
        set_frame_or_top_level_icon(t, Path(__file__).resolve().parent.parent / 'icons')
        t.title(title)
        Label(t, justify=LEFT, text=text).pack(side=TOP, expand=True, fill=BOTH, **self.pad)

        def dest():
            t.destroy()
        Button(t, text="OK", command=dest).pack(side=TOP, **self.pad)
        set_dialog_geometry(t, parent_window)
        t.grab_set()
        t.transient(parent_window)
        parent_window.wait_window(t)

    def display_with_alt_button(
            self, parent_window: Tk, title: str, text: str, button_text: str, button_action: Callable
    ):
        t = Toplevel(parent_window)
        t.title(title)
        Label(t, justify=LEFT, text=text).pack(side=TOP, expand=True, fill=BOTH, **self.pad)

        def dest():
            t.destroy()
        Button(t, text=button_text, command=button_action).pack(side=TOP, **self.pad)
        Button(t, text="OK", command=dest).pack(side=TOP, **self.pad)
        set_dialog_geometry(t, parent_window)
        t.grab_set()
        t.transient(parent_window)
        parent_window.wait_window(t)
