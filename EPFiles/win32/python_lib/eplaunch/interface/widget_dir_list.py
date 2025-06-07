from pathlib import Path
from platform import system
from string import ascii_uppercase
from tkinter import Tk, NSEW, VERTICAL, HORIZONTAL, Frame, END, NS, TOP, BOTH, EW, PhotoImage, BROWSE
from tkinter.ttk import Treeview, Scrollbar
from typing import Callable, List, Optional


class DirListWidget(Treeview):
    def __init__(
            self, parent_frame: Frame,
            on_select: Optional[Callable[[Path], None]] = None
    ):
        super().__init__(parent_frame, selectmode=BROWSE)
        self.heading('#0', text='Directory Selection Tree')
        # disable callbacks while widget is initializing
        self.callback_on_new_selection = None
        # set some paths to icons
        this_file_dir = Path(__file__).parent.resolve()
        folder_icon_path = this_file_dir / 'resources' / 'folder_opened.png'
        self.root_folder_image = PhotoImage(file=folder_icon_path)
        folder_closed_icon_path = this_file_dir / 'resources' / 'folder_closed.png'
        self.non_root_folder_image = PhotoImage(file=folder_closed_icon_path)
        # load the drive root(s), initially closed, with a single dummy item inside, so it allows expanding
        self.root_node_ids = []
        for r in DirListWidget.get_roots_by_platform():
            drive_root = self.insert('', END, text=r.anchor, open=False, tags=r.parts)
            self.root_node_ids.append(drive_root)
            self.insert(drive_root, END, text='{loading}', open=False, image=self.non_root_folder_image)
        # bind the click event
        self.bind('<<TreeviewSelect>>', self._item_selected)
        self.bind('<<TreeviewOpen>>', self._item_expanded)
        # connect bindings
        self.callback_on_new_selection = on_select

    def _item_expanded(self, *_):
        item_id = self.focus()
        self._rebuild_single_node(item_id)

    def _rebuild_single_node(self, item_id: str):
        item_details = self.item(item_id)
        item_path = Path(*item_details['tags'])
        # first delete any items from this node
        for item in self.get_children(item_id):
            self.delete(item)
        # now rebuild it
        raw_child_paths = []
        try:
            for f in item_path.iterdir():
                if f.is_dir():
                    raw_child_paths.append(f)
        except PermissionError:
            self.insert(item_id, END, text="{access_denied}", open=False, image=self.non_root_folder_image)
            return
        if len(raw_child_paths) == 0:
            self.insert(item_id, END, text="{empty}", open=False, image=self.non_root_folder_image)
        else:
            new_child_file_paths = sorted(raw_child_paths)
            for path in new_child_file_paths:
                new_child_id = self.insert(
                    item_id, END, text=path.name, open=False, image=self.non_root_folder_image, tags=path.parts
                )
                self.insert(new_child_id, END, text="{loading}", open=False, image=self.non_root_folder_image)

    @staticmethod
    def get_roots_by_platform() -> List[Path]:
        platform_name = system()
        if platform_name == 'Windows':
            from ctypes import windll
            bitmask = windll.kernel32.GetLogicalDrives()
            roots: List[Path] = []
            for letter in ascii_uppercase:
                if bitmask & 1:
                    roots.append(Path(f"{letter}:\\"))
                bitmask >>= 1
        else:  # Linux/Mac
            roots = [Path('/')]
        return roots

    def _item_selected(self, *_):
        if len(self.selection()) != 1:
            return  # must not have selected anything
        single_selected_item_id = self.selection()[0]
        item_contents = self.item(single_selected_item_id)
        if self.callback_on_new_selection:
            selected_path = Path(*item_contents['tags'])
            self.callback_on_new_selection(selected_path)

    def try_to_select_directory(self, target_path: Path):
        path_parts = target_path.parts
        # loop over each root, looking for the root node that matches the root of the target path
        # this is treated a little different because root paths could look weird on Windows, so
        # instead of a string match, I'm doing a pathlib.Path() equality check
        for found_root in self.root_node_ids:
            root_item = self.item(found_root)
            root_name = root_item['tags'][0]
            if Path(root_name) == Path(target_path.anchor):
                break
        else:
            # Could not find the root directory, just select the root
            raise Exception("Not an exception eventually")
        latest_parent_id = found_root
        items_to_expand = [latest_parent_id]
        self._rebuild_single_node(latest_parent_id)
        for ind, part in enumerate(path_parts):
            if ind == 0:
                continue  # we already got the root
            for child_node_id in self.get_children(latest_parent_id):
                this_child_node = self.item(child_node_id)
                if this_child_node['text'] == part:
                    latest_parent_id = child_node_id
                    items_to_expand.append(latest_parent_id)
                    self._rebuild_single_node(latest_parent_id)
                    break
        for i in items_to_expand:
            self.item(i, open=True)
        self.selection_set(latest_parent_id)
        self.see(latest_parent_id)


class DirListScrollableFrame(Frame):
    def __init__(
            self, parent,
            on_select: Optional[Callable[[Path], None]] = None
    ):
        super().__init__(parent)
        self.dir_list = DirListWidget(self, on_select)
        self.dir_list.grid(row=0, column=0, sticky=NSEW)
        ysb = Scrollbar(self, orient=VERTICAL, command=self.dir_list.yview)
        xsb = Scrollbar(self, orient=HORIZONTAL, command=self.dir_list.xview)
        self.dir_list.configure(yscrollcommand=ysb.set, xscrollcommand=xsb.set)
        ysb.grid(row=0, column=1, sticky=NS)
        xsb.grid(row=1, column=0, sticky=EW)
        self.grid_rowconfigure(0, weight=1)
        self.grid_columnconfigure(0, weight=1)


if __name__ == "__main__":
    root = Tk()
    root.title('Directory Listing Widget Demo')
    file_listing = DirListScrollableFrame(root)
    file_listing.pack(side=TOP, expand=True, fill=BOTH)
    file_listing.dir_list.try_to_select_directory(Path('/eplus/repos/1eplus'))
    root.mainloop()
