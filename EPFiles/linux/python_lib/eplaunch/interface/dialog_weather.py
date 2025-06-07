from pathlib import Path
from platform import system
from tkinter import Tk, Toplevel, StringVar, EW, Entry, ALL, DISABLED, filedialog, Label, OptionMenu, E, HORIZONTAL
if system() == 'Darwin':
    from tkmacosx import Button
else:
    from tkinter.ttk import Button
from tkinter.ttk import Separator
from typing import List, Optional

from eplaunch.interface import set_dialog_geometry, set_frame_or_top_level_icon


class TkWeatherDialog(Toplevel):
    CLOSE_SIGNAL_OK = 0
    CLOSE_SIGNAL_CANCEL = 1
    WEATHER_TYPE_DD = "DD"
    WEATHER_TYPE_EPW = "EPW"

    def __init__(self, parent_window, recent_files: List[Path], text: Optional[str] = None):
        super().__init__()
        set_frame_or_top_level_icon(self, Path(__file__).resolve().parent.parent / 'icons')
        self.alt_text = text
        self.title("Choose Weather Configuration")
        # assume cancel to allow for closing the dialog with the X
        self.exit_code = self.CLOSE_SIGNAL_CANCEL
        self.selected_weather_file: Optional[Path] = None
        self._epw_from_select = None
        # build the gui and call required modal methods
        self._define_tk_variables()
        self._build_gui(recent_files)
        set_dialog_geometry(self, parent_window)
        self.grab_set()
        self.transient(parent_window)

    def _define_tk_variables(self):
        self._tk_var_recent = StringVar(value="")
        self._tk_var_recent.trace('w', self._recent_changed)
        self._tk_var_select_epw_name = StringVar(value="<select_an_epw>")

    def _build_gui(self, recent_files: List[Path]):
        pad = {'padx': 3, 'pady': 3}
        button_text = "Choose This Configuration"
        if self.alt_text:
            Label(self, text=self.alt_text).grid(row=0, column=0, columnspan=3, sticky=EW, **pad)
            Separator(self, orient=HORIZONTAL).grid(row=1, column=0, columnspan=3, sticky=EW, **pad)
        Label(self, text="No Weather File (Design Run)").grid(row=2, column=0, columnspan=2, sticky=E, **pad)
        Button(self, text=button_text, command=self._go_dd).grid(row=2, column=2, sticky=EW, **pad)
        Separator(self, orient=HORIZONTAL).grid(row=3, column=0, columnspan=3, sticky=EW, **pad)
        if recent_files:
            self._epw_from_recent = None
            self._recent_file_paths = recent_files
            self._recent_file_names = [x.name for x in recent_files]
            Label(self, text="Use Recent").grid(row=4, column=0, sticky=E, **pad)
            self.recent_options = OptionMenu(self, self._tk_var_recent, *self._recent_file_names)
            self.recent_options.grid(row=4, column=1, sticky=EW, **pad)
            Button(self, text=button_text, command=self._go_recent).grid(row=4, column=2, sticky=EW, **pad)
            Separator(self, orient=HORIZONTAL).grid(row=5, column=0, columnspan=3, sticky=EW, **pad)
            self._tk_var_recent.set(str(self._recent_file_names[0]))
        Button(self, text="Select EPW", command=self._select_epw).grid(row=6, column=0, sticky=EW, **pad)
        Entry(self, textvariable=self._tk_var_select_epw_name, state=DISABLED).grid(row=6, column=1, sticky=EW, **pad)
        Button(self, text=button_text, command=self._go_select).grid(row=6, column=2, sticky=EW, **pad)
        Separator(self, orient=HORIZONTAL).grid(row=7, column=0, columnspan=3, sticky=EW, **pad)
        Button(self, text="Cancel", command=self._cancel).grid(row=8, column=0, columnspan=3, **pad)
        self.grid_rowconfigure(ALL, weight=1)
        self.grid_columnconfigure(ALL, weight=1)

    def _recent_changed(self, *_):
        selected_name = self._tk_var_recent.get()
        index_in_list = self._recent_file_names.index(selected_name)
        full_path = self._recent_file_paths[index_in_list]
        self._epw_from_recent = full_path

    def _select_epw(self):
        filetypes = (('EPW Weather Files', '*.epw'), ('All files', '*.*'))
        response = filedialog.askopenfilename(title="Specify Weather File", filetypes=filetypes)
        if response is not None:
            epw_path = Path(response)
            self._epw_from_select = Path(response)
            self._tk_var_select_epw_name.set(epw_path.name)

    def _go_dd(self):
        self.selected_weather_file = None
        self.exit_code = self.CLOSE_SIGNAL_OK
        self.grab_release()
        self.destroy()

    def _go_select(self):
        self.selected_weather_file = self._epw_from_select
        self.exit_code = self.CLOSE_SIGNAL_OK
        self.grab_release()
        self.destroy()

    def _go_recent(self):
        self.selected_weather_file = self._epw_from_recent
        self.exit_code = self.CLOSE_SIGNAL_OK
        self.grab_release()
        self.destroy()

    def _cancel(self):
        self.exit_code = self.CLOSE_SIGNAL_CANCEL
        self.grab_release()
        self.destroy()


if __name__ == "__main__":
    root = Tk()
    root.title('Root Window for Toplevel Demo')
    t = TkWeatherDialog(root, [Path('/hello/world.epw'), Path('/foo/bar.epw')], "hmm")
    root.wait_window(t)
    if t.exit_code == TkWeatherDialog.CLOSE_SIGNAL_CANCEL:
        print("Form was cancelled")
    else:
        if t.selected_weather_file is None:
            print("User selected DD ONLY")
        else:
            print(f"User selected this weather file: {t.selected_weather_file}")
