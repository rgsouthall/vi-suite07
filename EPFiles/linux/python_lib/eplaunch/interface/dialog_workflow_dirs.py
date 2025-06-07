from pathlib import Path
from tkinter import Tk, Toplevel, Frame, Button, EW, NSEW, Listbox, Variable, SINGLE, Scrollbar, LEFT, BOTH, RIGHT, ALL
from tkinter.filedialog import askdirectory
from tkinter.messagebox import showerror
from typing import List

from eplaunch.interface import set_dialog_geometry, set_frame_or_top_level_icon


class TkWorkflowsDialog(Toplevel):
    CLOSE_SIGNAL_OK = 0
    CLOSE_SIGNAL_CANCEL = 1

    def __init__(self, parent_window, current_workflow_dirs: List[Path], auto_found_workflow_dirs: List[Path]):
        super().__init__()
        set_frame_or_top_level_icon(self, Path(__file__).resolve().parent.parent / 'icons')
        self.title("Choose Workflow Directories")
        # assume cancel to allow for closing the dialog with the X
        self.exit_code = self.CLOSE_SIGNAL_CANCEL
        self.list_of_directories: List[Path] = []
        self._auto_find_dirs = auto_found_workflow_dirs
        self._tk_var_all_directories = Variable(value=current_workflow_dirs)
        self.selected_index = -1
        # build the gui and call required modal methods
        self._build_gui()
        set_dialog_geometry(self, parent_window)
        self.grab_set()
        self.transient(parent_window)

    def _update_from_traces(self, _=None):
        selected_indices = self.listbox.curselection()
        if len(selected_indices) != 1:
            self.selected_index = -1
            return
        self.selected_index = selected_indices[0]

    def _build_gui(self):
        f = Frame(self)
        self.listbox = Listbox(f, height=5, listvariable=self._tk_var_all_directories, selectmode=SINGLE)
        self.listbox.configure(exportselection=False)
        self.listbox.pack(side=LEFT, fill=BOTH, expand=True, padx=4, pady=4)
        self.listbox.bind('<<ListboxSelect>>', self._update_from_traces)
        columns_scroll = Scrollbar(f)
        columns_scroll.pack(side=RIGHT, fill=BOTH)
        self.listbox.config(yscrollcommand=columns_scroll.set)
        columns_scroll.config(command=self.listbox.yview)
        f.grid(row=0, column=0, sticky=NSEW, padx=4, pady=4)
        # f.grid_columnconfigure(ALL, weight=1)
        action_frame = Frame(self)
        Button(action_frame, text="Add Workflow Dir", command=self._add).grid(row=0, column=0, padx=4, pady=4)
        Button(action_frame, text="Remove Workflow Dir", command=self._remove).grid(row=0, column=1, padx=4, pady=4)
        Button(action_frame, text="Auto-find Workflows", command=self._auto).grid(row=0, column=2, padx=4, pady=4)
        action_frame.grid(row=1, column=0, sticky=EW, padx=4, pady=4)
        action_frame.grid_columnconfigure(ALL, weight=1)
        bottom_frame = Frame(self)
        Button(bottom_frame, text="OK", command=self._ok).grid(row=0, column=0, padx=4, pady=4)
        Button(bottom_frame, text="Cancel", command=self._cancel).grid(row=0, column=1, padx=4, pady=4)
        bottom_frame.grid_columnconfigure(ALL, weight=1)
        bottom_frame.grid(row=2, column=0, sticky=EW, padx=4, pady=4)
        self.grid_columnconfigure(ALL, weight=1)
        self.grid_rowconfigure(0, weight=5)
        self.grid_rowconfigure(1, weight=1)
        self.grid_rowconfigure(2, weight=1)

    def _add(self):
        response = askdirectory(mustexist=True)
        if response is None:
            return
        this_list = list(self._tk_var_all_directories.get())
        this_list.append(Path(response))
        self._tk_var_all_directories.set(this_list)
        self._update_from_traces()

    def _remove(self):
        if self.selected_index == -1:
            showerror('Selection Error', 'Could not remove directory, is a row selected?')
            return
        this_list = list(self._tk_var_all_directories.get())
        del this_list[self.selected_index]
        self._tk_var_all_directories.set(this_list)
        self._update_from_traces()

    def _auto(self):
        this_list = list(self._tk_var_all_directories.get())
        for auto_item in self._auto_find_dirs:
            if str(auto_item) not in this_list:
                this_list.append(auto_item)
        self._tk_var_all_directories.set(this_list)
        self._update_from_traces()

    def _ok(self):
        self.exit_code = self.CLOSE_SIGNAL_OK
        self.list_of_directories = [Path(p) for p in self._tk_var_all_directories.get()]
        self.grab_release()
        self.destroy()

    def _cancel(self):
        self.exit_code = self.CLOSE_SIGNAL_CANCEL
        self.grab_release()
        self.destroy()


if __name__ == "__main__":
    root = Tk()
    root.title('Root Window for Toplevel Demo')
    listing_set = [Path('/eplus/repos/1eplus'), Path('/home/edwin')]
    auto_find = [Path('/home/edwin/Documents')]
    file_listing = TkWorkflowsDialog(root, listing_set, auto_find)
    root.mainloop()
