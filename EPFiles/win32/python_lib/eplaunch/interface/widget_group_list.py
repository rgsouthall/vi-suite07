from tkinter import Listbox, Frame, Menu, END
from tkinter.ttk import LabelFrame
from typing import Callable, List, Union


class GroupListBox(Listbox):

    def __init__(self, parent_frame: Union[Frame, LabelFrame], delete_callback: Callable[[List[int]], None], *a, **kw):
        super().__init__(parent_frame, *a, **kw)
        self.cb_delete = delete_callback
        self.bind("<Button-3>", self._right_click)
        self.menu_already_open = False

    def _right_click(self, event):
        if self.menu_already_open:
            self.close_menu()
        self.selection_clear(0, END)
        self.selection_set(self.nearest(event.y))
        self.activate(self.nearest(event.y))
        self._popup_menu = Menu(self, tearoff=False)
        self._popup_menu.add_command(label="Delete", command=self.delete_selected)
        try:
            self._popup_menu.tk_popup(event.x_root, event.y_root, 0)
            self.menu_already_open = True
        finally:
            self._popup_menu.grab_release()
        self._popup_menu.bind("<FocusOut>", self.close_menu)

    def close_menu(self, *_):
        self._popup_menu.destroy()
        self.menu_already_open = False

    def delete_selected(self):
        self.cb_delete(self.curselection())


if __name__ == "__main__":
    from tkinter import Tk

    root = Tk()
    f = Frame(root)
    f.pack()
    flb = GroupListBox(f, lambda x: None)
    for n in range(10):
        flb.insert('end', n)
    flb.pack()
    root.mainloop()
