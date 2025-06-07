from platform import system
from tkinter import Tk, Toplevel, Label, Listbox, SINGLE, Variable, Entry, DISABLED, StringVar, Button, Frame, \
    filedialog, E, EW, W, LEFT, messagebox
from typing import Dict, List, Optional
from pathlib import Path

from eplaunch.interface import set_dialog_geometry, set_frame_or_top_level_icon


class TkViewerDialog(Toplevel):
    CLOSE_SIGNAL_OK = 0
    CLOSE_SIGNAL_CANCEL = 1

    def __init__(self, parent_window, list_of_suffixes: List[str], dict_of_viewer_overrides: Dict[str, Optional[Path]]):
        super().__init__()
        set_frame_or_top_level_icon(self, Path(__file__).resolve().parent.parent / 'icons')
        self.title("Viewers")
        self.exit_code = self.CLOSE_SIGNAL_CANCEL
        self.suffixes = ['txt'] + list_of_suffixes
        self.extension_to_viewer = dict_of_viewer_overrides
        self._define_tk_variables()
        self._build_gui()
        set_dialog_geometry(self, parent_window)
        self.grab_set()
        self.transient(parent_window)

    def _define_tk_variables(self):
        self._tk_var_extension_list = Variable()
        self._tk_var_application_path = StringVar()

    def _build_gui(self):
        extension_list = []
        # only keep the extensions for each suffix that are unique
        for suffix in self.suffixes:
            if '.' in suffix:
                s = suffix.split('.')[-1]
            else:
                s = suffix
            if s not in extension_list:
                extension_list.append(suffix)
        self._tk_var_extension_list.set(extension_list)
        for extension in extension_list:
            if extension not in self.extension_to_viewer:
                self.extension_to_viewer[extension] = None

        Label(self, text="Extensions").grid(row=0, column=0, padx=3, pady=3)
        # self._tk_var_extension_list.set(value=self.extension_to_viewer)
        self._lb_extensions = Listbox(self, height=10, listvariable=self._tk_var_extension_list, selectmode=SINGLE)
        self._lb_extensions.configure(exportselection=False)
        self._lb_extensions.grid(row=1, column=0, rowspan=3, padx=3, pady=3)
        self._lb_extensions.bind('<<ListboxSelect>>', self._update_for_selected_extension)

        Label(self, text="Application Path").grid(row=0, column=1, padx=3, pady=3)
        Entry(self, textvariable=self._tk_var_application_path, state=DISABLED).grid(
            row=1, column=1, padx=3, pady=3, sticky=EW
        )
        button_mid_frame = Frame(self)
        Button(button_mid_frame, text="Default", command=self._set_as_default).grid(row=0, column=0, padx=3, pady=3)
        Button(button_mid_frame, text="Select...", command=self._select_path).grid(row=0, column=1, padx=3, pady=3)
        button_mid_frame.grid(row=2, column=1, padx=3, pady=3, sticky=W)
        label_string_lines = [
            "Only set applications for extensions that don't open with",
            "the desired applications automatically. Typically,",
            "extension 'htm' opens using a web browser, extension",
            "'csv' opens using a spreadsheet program, extension 'txt'",
            "opens using a text editor, etc. By default, the 'txt'",
            "extension also opens all non-typical extensions such as",
            "'err' and 'eio'."
        ]
        Label(self, text='\n'.join(label_string_lines), justify=LEFT).grid(row=3, column=1, padx=3, pady=3, sticky=W)
        bottom_button_frame = Frame(self)
        Button(bottom_button_frame, text="OK", command=self._ok).grid(row=0, column=1, padx=3, pady=3)
        Button(bottom_button_frame, text="Cancel", command=self._cancel).grid(row=0, column=2, padx=3, pady=3)
        bottom_button_frame.grid(row=4, column=1, padx=3, pady=3, sticky=E)

    def _ok(self):
        self.exit_code = self.CLOSE_SIGNAL_OK
        current_mapping = self.extension_to_viewer
        self.extension_to_viewer = {k: v for k, v in current_mapping.items() if v is not None}
        self.grab_release()
        self.destroy()

    def _cancel(self):
        self.exit_code = self.CLOSE_SIGNAL_CANCEL
        self.grab_release()
        self.destroy()

    def _set_as_default(self):
        currently_selected = self._lb_extensions.curselection()
        if len(currently_selected) != 1:
            messagebox.showerror("Selection Error", "Must select one extension")
            return
        current_extension_index = currently_selected[0]
        current_extension = self.suffixes[current_extension_index]
        if current_extension in self.extension_to_viewer:
            self.extension_to_viewer[current_extension] = None
        self._update_for_selected_extension()

    def _select_path(self):
        currently_selected = self._lb_extensions.curselection()
        if len(currently_selected) != 1:
            messagebox.showerror("Selection Error", "Must select one extension")
            return
        exe = "*.exe" if system() == "Windows" else "*"
        response = filedialog.askopenfilename(title="Specify Weather File", filetypes=(('Programs', exe),))
        if response is None or response == '':
            return
        current_extension_index = currently_selected[0]
        current_extension = self.suffixes[current_extension_index]
        self._tk_var_application_path.set(response)
        self.extension_to_viewer[current_extension] = Path(response)

    def _update_for_selected_extension(self, _=None):
        currently_selected = self._lb_extensions.curselection()
        if len(currently_selected) != 1:
            messagebox.showerror("Selection Error", "Must select one extension")
            return
        current_extension_index = currently_selected[0]
        current_extension = self.suffixes[current_extension_index]
        viewer = self.extension_to_viewer.get(current_extension, None)
        if viewer is None:
            text = '<DEFAULT>'
        else:
            text = str(viewer)
        self._tk_var_application_path.set(text)


if __name__ == "__main__":
    root = Tk()
    root.title('Root Window for Toplevel Demo')
    file_listing = TkViewerDialog(root, ['doc', 'rdd', 'csv', 'err'], {'doc': Path('/word'), 'csv': Path('/excel')})
    root.mainloop()
