from datetime import datetime
from fnmatch import fnmatch

from json import dumps
from mimetypes import guess_type
from os import getenv
from pathlib import Path
from platform import system
from queue import Queue
from sys import version_info

from subprocess import Popen
from tkinter import Tk, PhotoImage, StringVar, Menu, DISABLED, Frame, Label, NSEW, E, VERTICAL, \
    SUNKEN, S, LEFT, BOTH, messagebox, END, BooleanVar, NORMAL, RIGHT, EW, NS, filedialog, \
    ALL, Scrollbar, SINGLE, Variable, HORIZONTAL
from tkinter.ttk import Combobox, PanedWindow as ttkPanedWindow, OptionMenu, LabelFrame
if system() == 'Darwin':
    from tkmacosx import Button
else:
    from tkinter.ttk import Button
from typing import Dict, List, Optional, Tuple
from uuid import uuid4
from webbrowser import open as open_web
from plan_tools.runtime import fixup_taskbar_icon_on_windows

from eplaunch import VERSION, DOCS_URL, NAME
from eplaunch.interface.dialog_generic import TkGenericDialog
from eplaunch.utilities.exceptions import EPLaunchFileException
from eplaunch.interface.config import ConfigManager
from eplaunch.interface.widget_dir_list import DirListScrollableFrame
from eplaunch.interface.widget_file_list import FileListScrollableFrame
from eplaunch.interface.widget_group_list import GroupListBox
from eplaunch.interface.dialog_external_viewers import TkViewerDialog
from eplaunch.interface.dialog_weather import TkWeatherDialog
from eplaunch.interface.dialog_workflow_dirs import TkWorkflowsDialog
from eplaunch.interface.dialog_output import TkOutputDialog
from eplaunch.workflows.workflow_thread import WorkflowThread
from eplaunch.utilities.cache import CacheFile
from eplaunch.workflows.base import BaseEPLaunchWorkflow1, EPLaunchWorkflowResponse1
from eplaunch.workflows.manager import WorkflowManager
from eplaunch.workflows.workflow import Workflow


class EPLaunchWindow(Tk):

    # region Construction and GUI building functions

    def __init__(self, called_from_ep_cli: bool):
        fixup_taskbar_icon_on_windows(NAME)
        super().__init__(className=NAME)
        # set the form title and icon, basic stuff
        self.title("EnergyPlus Launch")
        self.option_add('*Dialog.msg.font', 'Helvetica 12')
        if system() == 'Darwin':
            self.icon_path = Path(__file__).resolve().parent.parent / 'icons' / 'icon.icns'
            if self.icon_path.exists():
                self.iconbitmap(str(self.icon_path))
            else:
                print(f"Could not set icon for Mac, expecting to find it at {self.icon_path}")
        elif system() == 'Windows':
            self.icon_path = Path(__file__).resolve().parent.parent / 'icons' / 'icon.png'
            img = PhotoImage(file=str(self.icon_path))
            if self.icon_path.exists():
                self.iconphoto(False, img)
            else:
                print(f"Could not set icon for Windows, expecting to find it at {self.icon_path}")
        else:  # Linux
            self.icon_path = Path(__file__).resolve().parent.parent / 'icons' / 'icon.png'
            img = PhotoImage(file=str(self.icon_path))
            if self.icon_path.exists():
                self.iconphoto(False, img)
            else:
                print(f"Could not set icon for Windows, expecting to find it at {self.icon_path}")
        self.pad = {'padx': 3, 'pady': 3}
        self.dd_only_string = '<No_Weather_File>'  # Can we change to just using blank?  It's fine for now.
        self.generic_dialogs = TkGenericDialog(self.pad)

        # check the version of python for whether we can show extended UTF characters
        self.extended_utf8 = version_info >= (3, 7, 6)

        # initialize variables which will track output dialogs
        self.dialog_counter: int = 0
        self.output_dialogs: Dict[str, TkOutputDialog] = {}

        # create a config manager and load up the saved, or default, configuration
        self.conf = ConfigManager()
        self.conf.load(called_from_ep_cli)

        # initialize some dir/file selection variables
        self.previous_selected_directory: Optional[Path] = None

        # create a workflow manager, it will initialize workflows in predetermined locations
        self.workflow_manager = WorkflowManager()

        # but now, if the saved configuration exists, use that as the list of directories to use moving forward
        if self.conf.workflow_directories:
            self.workflow_manager.workflow_directories = self.conf.workflow_directories

        # now that we have a list of workflows, instantiate any/all of them
        self.workflow_manager.instantiate_all_workflows()
        if len(self.workflow_manager.warnings) > 0:
            self.generic_dialogs.display(
                self, 'Workflow Processing Errors', '\n'.join(self.workflow_manager.warnings)
            )

        # set up the GUI update queue
        self._gui_queue = Queue()
        self._check_queue()

        # set up some tk tracking variables and then build the GUI itself
        self._define_tk_variables()
        self._build_gui()

        # start filling in the UI with contents, optionally from the previous config
        self.available_workflows: List[Workflow] = []  # stores workflows under the current context
        self._init = True
        self._repopulate_workflow_context_menu()
        self._repopulate_workflow_instance_menu()
        self._repopulate_recent_weather_list()
        self._rebuild_recent_folder_menu()
        self._rebuild_favorite_folder_menu()
        self._rebuild_group_menu()
        self._init = False

        # finally set the initial directory and update the file listing
        self.dir_tree.dir_list.try_to_select_directory(self.conf.directory)
        self._update_file_list()

        # set the minimum size and redraw the app
        self.minsize(1050, 400)
        self.update()  # one quick redraw should be fine here to get updated geometry
        self.dir_files_pw.sashpos(0, self.conf.dir_file_paned_window_sash_position)
        self.conf.list_group_paned_window_sash_position = min(
            self.conf.list_group_paned_window_sash_position, self.winfo_height() - 250
        )
        self.list_group_pw.sashpos(0, self.conf.list_group_paned_window_sash_position)

        # one time update of the status bar
        self._update_status_bar("Program Initialized")

        # potentially show a welcome screen if we are on a new version
        self._open_welcome()

        # bind key presses for the app
        self.bind('<Key>', self._handle_keyboard_press)

        # and bind the focus event for the app
        self.bind("<FocusIn>", self._handle_focus_in)

    def _define_tk_variables(self):
        self._tk_var_workflow_context = StringVar(value='<context>')
        self._tk_var_workflow_instance = StringVar(value='<instance>')
        self._tk_var_output_suffix = StringVar(value='<output')
        self._tk_var_output_1 = StringVar(value='--')
        self._tk_var_output_2 = StringVar(value='--')
        self._tk_var_output_3 = StringVar(value='--')
        self._tk_var_output_4 = StringVar(value='--')
        self._tk_var_output_5 = StringVar(value='--')
        self._tk_var_status_dir = StringVar(value="Selected dir")
        self._tk_var_status_workflow = StringVar(value="Selected workflow")
        self._tk_var_status_message = StringVar(value="Status Message")
        self._tk_var_weather_recent = StringVar(value="<recent>")

    def _build_gui(self):
        self._build_top_menu()

        top_action_bar = Frame(self)
        self._build_top_icon_bar(top_action_bar)
        top_action_bar.grid(row=0, column=0, sticky=NSEW)

        listing_frame = Frame(self)
        self._build_listings(listing_frame)
        listing_frame.grid(row=1, column=0, sticky=NSEW)

        status_frame = Frame(self)
        self._build_status_bar(status_frame)
        status_frame.grid(row=2, column=0, sticky=NSEW)

        self.grid_rowconfigure(1, weight=1)
        self.grid_columnconfigure(0, weight=1)

        height = max(self.conf.height, 500)
        width = max(self.conf.width, 1000)
        x = max(self.conf.x, 128)
        y = max(self.conf.y, 128)
        self.wm_geometry(f"{width}x{height}+{x}+{y}")

    def _build_top_menu(self):
        menubar = Menu(self)

        menu_file = Menu(menubar, tearoff=False)
        menu_file.add_command(label="Run Current Workflow on Selection", command=self._run_workflow_on_selection)
        menu_file.add_command(label="Run Current Workflow on Current Group", command=self._run_workflow_on_group)
        menu_file.add_separator()
        menu_file.add_command(label="Quit", command=self._window_close)
        menubar.add_cascade(label="File", menu=menu_file)

        menu_nav = Menu(menubar, tearoff=False)
        self.menu_nav_recent = Menu(menu_nav, tearoff=False)
        menu_nav.add_cascade(label="Recent", menu=self.menu_nav_recent)
        menu_nav.add_command(label="Navigate to previous folder", command=self._navigate_to_previous_folder)
        menu_nav.add_separator()
        self.menu_nav_favorites = Menu(menu_nav, tearoff=False)
        menu_nav.add_cascade(label="Favorites", menu=self.menu_nav_favorites)
        menu_nav.add_command(label="Add Current Folder to Favorites", command=self._add_folder_to_favorites)
        menu_nav.add_command(label="Remove Current Folder from Favorites", command=self._remove_folder_from_favorites)
        menu_nav.add_separator()
        menu_nav.add_command(label="View Keyboard Shortcuts", command=self._open_shortcuts_dialog)
        menubar.add_cascade(label="Navigation", menu=menu_nav)

        menu_group = Menu(menubar, tearoff=False)
        menu_group.add_command(
            label="Add selected files to current group", command=self._add_current_files_to_group
        )
        menu_group.add_command(label="Clear Current Group", command=self._clear_group)
        menu_group.add_separator()
        menu_group.add_command(
            label="Set weather for current group", command=self._set_weather_for_current_group
        )
        menu_group.add_command(
            label="Cycle through current group entries", command=self._cycle_through_group
        )
        menu_group.add_separator()
        menu_group.add_command(label="Save Current Group to file...", command=self._save_group_file)
        menu_group.add_command(label="Load File to Current Group...", command=self._load_new_group_file)
        menubar.add_cascade(label="Group", menu=menu_group)

        menu_settings = Menu(menubar, tearoff=False)
        menu_settings.add_command(label="Workflow Directories", command=self._open_workflow_dir_dialog)
        self._tk_var_keep_dialogs_open = BooleanVar(value=True)
        menu_settings.add_checkbutton(
            label="Keep Output Dialog Open", onvalue=True, offvalue=False, variable=self._tk_var_keep_dialogs_open
        )

        def _update_keep_dialog_open(*_):
            """Called whenever the checkbox is checked, updates configuration value"""
            self.conf.keep_dialog_open = self._tk_var_keep_dialogs_open.get()

        self._tk_var_keep_dialogs_open.trace('w', _update_keep_dialog_open)
        menu_settings.add_command(label="Viewers...", command=self._open_viewers_dialog)
        menubar.add_cascade(label="Settings", menu=menu_settings)

        menu_help = Menu(menubar, tearoff=False)
        menu_help.add_command(label="EnergyPlus-Launch Documentation", command=self._open_documentation)
        menu_help.add_command(label="About...", command=self._open_about)
        menubar.add_cascade(label="Help", menu=menu_help)

        self.config(menu=menubar)

    # noinspection SqlNoDataSourceInspection
    def _build_top_icon_bar(self, container: Frame):
        lf = LabelFrame(container, text="Workflow Selection")
        Label(lf, text="Context:", justify=RIGHT).grid(row=0, column=0, sticky=E, **self.pad)
        self.option_workflow_context = OptionMenu(lf, self._tk_var_workflow_context, '<context>')
        self.option_workflow_context.grid(row=0, column=1, sticky=EW, **self.pad)
        Label(lf, text="Workflow:", justify=RIGHT).grid(row=1, column=0, sticky=E, **self.pad)
        self.option_workflow_instance = OptionMenu(lf, self._tk_var_workflow_instance, '<instance>')
        self.option_workflow_instance.grid(row=1, column=1, sticky=EW, **self.pad)
        lf.grid_rowconfigure(ALL, weight=1)
        lf.grid_columnconfigure(ALL, weight=1)
        lf.grid(row=0, column=0, sticky=NS, **self.pad)

        lf = LabelFrame(container, text="Run Workflow on...")
        if self.extended_utf8:
            b = Button(lf, text=u"\U000025B6 Selected File(s)", command=self._run_workflow_on_selection, state=NORMAL)
        else:
            b = Button(lf, text="Selected File(s)", command=self._run_workflow_on_selection, state=NORMAL)
        b.grid(row=0, column=0, sticky=EW, **self.pad)
        if self.extended_utf8:
            b = Button(lf, text=u"\U000025B6 Current Group", command=self._run_workflow_on_group, state=NORMAL)
        else:
            b = Button(lf, text="Current Group", command=self._run_workflow_on_group, state=NORMAL)
        b.grid(row=1, column=0, sticky=EW, **self.pad)
        lf.grid_rowconfigure(ALL, weight=1)
        lf.grid_columnconfigure(ALL, weight=1)
        lf.grid(row=0, column=1, sticky=NS, **self.pad)

        lf = LabelFrame(container, text="Weather Selection")
        Label(lf, text="Recent: ", justify=RIGHT).grid(row=0, column=0, **self.pad)
        self.option_weather_recent = OptionMenu(lf, self._tk_var_weather_recent, '<weather>')
        self.option_weather_recent.grid(row=0, column=1, sticky=EW, **self.pad)
        if self.extended_utf8:
            self.button_weather_select = Button(
                lf, text=u"\U0001f325 Select Weather File...", command=self._open_weather_dialog
            )
        else:
            self.button_weather_select = Button(
                lf, text="Select Weather File...", command=self._open_weather_dialog
            )
        self.button_weather_select.grid(row=1, column=0, columnspan=2, sticky=EW, **self.pad)
        lf.grid_rowconfigure(ALL, weight=1)
        lf.grid_columnconfigure(ALL, weight=1)
        lf.grid(row=0, column=2, sticky=NS, **self.pad)

        lf = LabelFrame(container, text="Open Selected File(s) or Directory")
        if self.extended_utf8:
            self.button_open_in_text = Button(
                lf, text=u"\U0001F5B9 Open File in Text Editor", command=self._open_text_editor, state=DISABLED
            )
        else:
            self.button_open_in_text = Button(
                lf, text="Open File in Text Editor", command=self._open_text_editor, state=DISABLED
            )
        self.button_open_in_text.grid(row=0, column=0, sticky=EW, **self.pad)
        if self.extended_utf8:
            dir_button = Button(
                lf, text=u"\U0001F5C0 Open Dir in File Browser", command=self._open_file_browser, state=NORMAL
            )
        else:
            dir_button = Button(
                lf, text="Open Dir in File Browser", command=self._open_file_browser, state=NORMAL
            )
        dir_button.grid(row=1, column=0, sticky=EW, **self.pad)
        idf_editor_icon_path = Path(__file__).resolve().parent / 'resources' / 'idf_editor_icon.png'
        self.idf_editor_image = PhotoImage(file=str(idf_editor_icon_path))
        self.button_idf_editor = Button(
            lf, image=self.idf_editor_image, command=self._open_idf_editor, state=DISABLED
        )
        self.button_idf_editor.grid(row=0, rowspan=2, column=1, sticky=EW, **self.pad)
        lf.grid_rowconfigure(ALL, weight=1)
        lf.grid_columnconfigure(ALL, weight=1)
        lf.grid(row=0, column=3, sticky=NS, **self.pad)

        lf = LabelFrame(container, text="Open Outputs")
        sub_frame = Frame(lf)
        self.button_open_output_1 = Button(
            sub_frame, textvariable=self._tk_var_output_1, command=self._open_output_1, state=DISABLED
        )
        self.button_open_output_1.grid(row=0, column=0, sticky=EW, **self.pad)
        self.button_open_output_2 = Button(
            sub_frame, textvariable=self._tk_var_output_2, command=self._open_output_2, state=DISABLED
        )
        self.button_open_output_2.grid(row=0, column=1, sticky=EW, **self.pad)
        self.button_open_output_3 = Button(
            sub_frame, textvariable=self._tk_var_output_3, command=self._open_output_3, state=DISABLED
        )
        self.button_open_output_3.grid(row=0, column=2, sticky=EW, **self.pad)
        self.button_open_output_4 = Button(
            sub_frame, textvariable=self._tk_var_output_4, command=self._open_output_4, state=DISABLED
        )
        self.button_open_output_4.grid(row=0, column=3, sticky=EW, **self.pad)
        self.button_open_output_5 = Button(
            sub_frame, textvariable=self._tk_var_output_5, command=self._open_output_5, state=DISABLED
        )
        self.button_open_output_5.grid(row=0, column=4, sticky=EW, **self.pad)
        sub_frame.grid_columnconfigure(ALL, weight=1)
        sub_frame.grid_rowconfigure(ALL, weight=1)
        sub_frame.grid(row=0, column=0, columnspan=3, sticky=EW, **self.pad)
        Label(lf, text="Full List: ", justify=RIGHT).grid(row=1, column=0, **self.pad)
        self.option_workflow_outputs = Combobox(lf, textvariable=self._tk_var_output_suffix)
        self.option_workflow_outputs.grid(row=1, column=1, **self.pad)
        self.option_workflow_outputs['state'] = 'readonly'
        if self.extended_utf8:
            self.button_open_output_file = Button(lf, text=u"\U0001f325 Open", command=self._open_output_file)
        else:
            self.button_open_output_file = Button(lf, text="Open", command=self._open_output_file)
        self.button_open_output_file.grid(row=1, column=2, **self.pad)
        lf.grid_columnconfigure(ALL, weight=1)
        lf.grid_rowconfigure(ALL, weight=1)
        lf.grid(row=0, column=4, sticky=NS, **self.pad)

    def _build_listings(self, container: Frame):
        # use vertical for the primary paned window
        self.list_group_pw = ttkPanedWindow(container, orient=VERTICAL)
        # create a sub paned window that is vertical to split the directory from the group box
        self.dir_files_pw = ttkPanedWindow(container, orient=HORIZONTAL)
        # the top part of this pane is simply the directory tree, add it with a heavier weight
        self.dir_tree = DirListScrollableFrame(self.dir_files_pw, on_select=self._new_dir_selected)
        self.dir_files_pw.add(self.dir_tree, weight=1)
        self.file_list = FileListScrollableFrame(
            self.dir_files_pw, on_selection_changed=self._callback_file_selection_changed
        )
        self.dir_files_pw.add(self.file_list, weight=1)
        self.list_group_pw.add(self.dir_files_pw, weight=5)
        # the bottom part of the primary pane will be a label frame containing a group listbox, add it with lower weight
        group_label_frame = LabelFrame(self.list_group_pw, text="Current Group")
        self.group_list_box = GroupListBox(group_label_frame, self._group_item_delete, height=5, selectmode=SINGLE)
        self.group_list_box.bind('<Double-1>', self._handle_group_selection_from_widget)
        from string import ascii_uppercase
        for i, e in enumerate(ascii_uppercase):
            self.group_list_box.insert(i + 1, e)
        self.group_list_box.pack(side=LEFT, fill=BOTH, expand=True, **self.pad)
        scroll = Scrollbar(group_label_frame)
        scroll.pack(side=RIGHT, fill=BOTH)
        self.group_list_box.config(yscrollcommand=scroll.set)
        scroll.config(command=self.group_list_box.yview)
        self.list_group_pw.add(group_label_frame, weight=1)
        # add the left paned window to the primary panes
        self.list_group_pw.pack(fill=BOTH, expand=True, **self.pad)

    def _build_status_bar(self, container: Frame):
        Label(container, relief=SUNKEN, anchor=S, textvariable=self._tk_var_status_dir).pack(
            side=LEFT, fill=BOTH, expand=True
        )
        Label(container, relief=SUNKEN, anchor=S, textvariable=self._tk_var_status_workflow).pack(
            side=LEFT, fill=BOTH, expand=True
        )
        Label(container, relief=SUNKEN, anchor=S, textvariable=self._tk_var_status_message).pack(
            side=LEFT, fill=BOTH, expand=True
        )

    def _check_queue(self):
        """Checks the GUI queue for actions and sets a timer to check again each time"""
        while True:
            # noinspection PyBroadException
            try:
                task = self._gui_queue.get(block=False)
                self.after_idle(task)
            except Exception:
                break
        self.after(100, self._check_queue)

    # endregion

    # region main window run/close/update functions

    def run(self):
        """Main entry point to register a close handler and execute the app main loop"""
        self.protocol('WM_DELETE_WINDOW', self._window_close)
        self.mainloop()

    def _window_close(self, *_):
        # block for running threads
        if len(self.workflow_manager.threads) > 0:
            msg = 'Program closing, but there are threads running; would you like to kill the threads and close?'
            if messagebox.askyesno(msg):
                for thread_id in self.workflow_manager.threads:
                    try:
                        self.workflow_manager.threads[thread_id].abort()
                        self.output_dialogs[thread_id].close()
                    except Exception as e:
                        print("Tried to abort thread, but something went awry: " + str(e))
            else:
                return  # will leave the window open
        # close completed windows
        keys = list(self.output_dialogs.keys())
        for dialog_id in keys:
            self.output_dialogs[dialog_id].close()
        # save the configuration and close
        self.conf.x = self.winfo_x()
        self.conf.y = self.winfo_y()
        self.conf.width = self.winfo_width()
        self.conf.height = self.winfo_height()
        self.conf.dir_file_paned_window_sash_position = self.dir_files_pw.sashpos(0)
        self.conf.list_group_paned_window_sash_position = self.list_group_pw.sashpos(0)
        self.conf.save()
        self.destroy()

    def _update_status_bar(self, message: str) -> None:
        """Updates all the status bar entries for current status and status by setting tk variables"""
        if self.conf.directory:
            self._tk_var_status_dir.set(str(self.conf.directory))
        else:
            self._tk_var_status_dir.set('<No Directory>')
        if self.workflow_manager.current_workflow:
            self._tk_var_status_workflow.set(self.workflow_manager.current_workflow.description)
        else:
            self._tk_var_status_workflow.set('<No Workflow>')
        self._tk_var_status_message.set(message)

    def _handle_keyboard_press(self, event) -> None:
        # Update docs/navigation.rst whenever this set of keybindings changes
        # Update self._list_keyboard_shortcuts() whenever this set of keybindings changes
        # relevant_modifiers
        # mod_shift = 0x1
        mod_control = 0x4
        # mod_alt = 0x20000
        if event.keysym == 'F5':
            self._update_file_list()
        elif event.keysym == 'w' and mod_control & event.state:
            self._open_weather_dialog()
        elif event.keysym == 'r' and mod_control & event.state:
            self._run_workflow_on_selection()
        elif event.keysym == 'g' and mod_control & event.state:
            self._run_workflow_on_group()
        elif event.keysym == 'm' and mod_control & event.state:
            self._cycle_through_group()
        elif event.keysym == 'z' and mod_control & event.state:
            self._navigate_to_previous_folder()

    def _handle_focus_in(self, _event) -> None:
        # This is firing twice, even on dummy hello world window examples.
        # The event argument is supposed to have a .detail member which you
        # can check to see if it equals NotifyAncestor and skip some callbacks.
        # Unfortunately that is missing, so there's not much we can do.
        # If we really needed to skip some, we could do a small timer that
        # won't let two callbacks happen within 0.5s of each other or something.
        # For now that's not necessary, but could be in the future.
        # As of right now, the only reason we need this is to try to refresh
        # the stale attribute, which should be accomplished with just the call here.
        # It's possible we also need to check the _event.widget to make sure we are only
        # calling this when the main window is focused.
        pass  # self._update_file_list()

    @staticmethod
    def _list_keyboard_shortcuts() -> List[Tuple[str, str]]:
        # Mimic the keyboard shortcuts in the function above
        # It helps the dialog text look nice if you make sure the widths are uniform
        return [
            ("Control + F5", "Update File List"),
            ("Control + w ", "Open Weather Dialog"),
            ("Control + r ", "Run Workflow on Selected Files"),
            ("Control + g ", "Run Workflow on Current Group"),
            ("Control + m ", "Cycle Through Group Files"),
            ("Control + z ", "Navigate to Previous Folder"),
        ]

    # endregion

    # region group operations

    def _clear_group(self):
        self.conf.group_locations.clear()
        self._rebuild_group_menu()

    def _load_new_group_file(self):
        group_file_path = filedialog.askopenfilename(
            title="Load a new EPLaunch Group File",
            initialdir=self.conf.directory,
            filetypes=(("EnergyPlus-Launch Group Files", "*.epg3"),)
        )
        if not group_file_path:
            return
        group_file_path_object = Path(group_file_path)
        entries = group_file_path_object.read_text().splitlines()
        self.conf.group_locations = []
        for e in sorted(entries):
            if e:
                self.conf.group_locations.append(Path(e))
        self._rebuild_group_menu()

    def _save_group_file(self):
        if len(self.conf.group_locations) == 0:
            messagebox.showerror("File Selection Issue", "Group is currently empty, add files to group before saving!")
            return
        group_file_path = filedialog.asksaveasfilename(
            title="Save EPLaunch Group File",
            initialdir=self.conf.directory,
            filetypes=(("EnergyPlus-Launch Group Files", "*.epg3"),)
        )
        if not group_file_path:
            return
        if not group_file_path.endswith('.epg3'):
            group_file_path += '.epg3'
        group_file_path_object = Path(group_file_path)
        try:
            group_file_path_object.write_text('\n'.join([str(p) for p in self.conf.group_locations]))
        except Exception as e:  # could be permission or just about anything
            messagebox.showerror("Error", f"Could not save group file, error message:\n{str(e)}")
            return

    def _add_current_files_to_group(self):
        if len(self.conf.file_selection) == 0:
            messagebox.showwarning("Warning", "No files selected, so none added to current group")
            return
        issue_warning = False
        for f in self.conf.file_selection:
            potential_path = self.conf.directory / f
            if potential_path in self.conf.group_locations:
                issue_warning = True
            else:
                self.conf.group_locations.append(potential_path)
        if issue_warning:
            messagebox.showwarning("Warning", "At least one file path was already in the group and skipped")
        self._rebuild_group_menu()

    def _group_item_delete(self, indices_to_delete):
        for i in indices_to_delete[::-1]:
            del self.conf.group_locations[i]
        self._rebuild_group_menu()

    def _cycle_through_group(self):
        if len(self.conf.group_locations) == 0:
            messagebox.showwarning("Invalid Group Cycle", "Could not cycle through empty group list")
            return
        self.group_list_box.selection_clear(0, END)
        self.group_list_box.selection_set(self.group_cycle_next_index)
        self._handle_group_selection(self.group_cycle_next_index)
        self.group_cycle_next_index += 1
        if self.group_cycle_next_index + 1 > len(self.conf.group_locations):
            self.group_cycle_next_index = 0

    def _get_weather_if_exists_in_cache(self, file_path: Path):
        if not self.workflow_manager.current_workflow.uses_weather:
            return ""
        cache = CacheFile(file_path.parent)
        files_in_current_workflow = cache.get_files_for_workflow(self.workflow_manager.current_workflow.name)
        if file_path.name in files_in_current_workflow:
            cached_file_info = files_in_current_workflow[file_path.name]
            if CacheFile.ParametersKey in cached_file_info:
                if CacheFile.WeatherFileKey in cached_file_info[CacheFile.ParametersKey]:
                    return cached_file_info[CacheFile.ParametersKey][CacheFile.WeatherFileKey]
        return ""

    def _rebuild_group_menu(self):
        self.group_cycle_next_index = 0
        self.group_list_box.delete(0, END)
        for index, entry in enumerate(self.conf.group_locations):
            known_weather_file_string = self._get_weather_if_exists_in_cache(entry)
            suffix = f"   --   {known_weather_file_string}" if known_weather_file_string else ""
            self.group_list_box.insert(index, f"{entry}{suffix}")

    def _handle_group_selection_from_widget(self, *_):
        self._handle_group_selection()

    def _handle_group_selection(self, specific_index: Optional[int] = None):
        if isinstance(specific_index, int):
            idx = specific_index
        else:
            idx = self.group_list_box.curselection()[0]
        try:
            entry = self.conf.group_locations[idx]
        except IndexError:
            messagebox.showerror("Error", f"Could not navigate to group entry at index {idx}")
            return
        if not entry.exists():
            messagebox.showerror("Error", f"Could not navigate to group entry, it doesn't exist: {entry}")
            return
        self.conf.directory = entry.parent
        self.dir_tree.dir_list.try_to_select_directory(self.conf.directory)
        self._update_file_list()
        self.file_list.tree.try_to_reselect([entry.name])

    # endregion

    # region handling file/folder navigation

    def _repopulate_control_list_columns(self):
        """Rebuilds the file list columns based on the current workflow status and column headers"""
        # add stale if the current workflow supports output suffixes
        column_list = []
        if len(self.workflow_manager.current_workflow.output_suffixes) > 0:
            column_list.append('Stale')
        # if the current workflow is set (should always be) extend the array dynamically
        if self.workflow_manager.current_workflow:
            if self.workflow_manager.current_workflow.uses_weather:
                column_list.append("Weather")
            column_list.extend(self.workflow_manager.current_workflow.columns)
        # ask the file listing widget to update columns
        self.file_list.tree.set_new_columns(column_list)

    def _add_folder_to_favorites(self):
        if self.conf.directory in self.conf.folders_favorite:
            messagebox.showwarning("Favorite Selection", "Folder already exists in favorites, skipping")
            return
        self.conf.folders_favorite.append(self.conf.directory)
        self._rebuild_favorite_folder_menu()

    def _remove_folder_from_favorites(self):
        if self.conf.directory not in self.conf.folders_favorite:
            messagebox.showwarning("Favorite Selection", "Folder is not in favorites, skipping")
            return
        self.conf.folders_favorite.remove(self.conf.directory)
        self._rebuild_favorite_folder_menu()

    def _rebuild_favorite_folder_menu(self):
        self.menu_nav_favorites.delete(0, END)
        for index, folder in enumerate(self.conf.folders_favorite):
            self.menu_nav_favorites.add_command(
                label=str(folder), command=lambda i=index: self._handle_favorite_folder_selection(i)
            )

    def _handle_favorite_folder_selection(self, folder_index: int):
        self.conf.directory = self.conf.folders_favorite[folder_index]
        self.dir_tree.dir_list.try_to_select_directory(self.conf.directory)
        self._update_file_list()

    def _rebuild_recent_folder_menu(self):
        self.menu_nav_recent.delete(0, END)
        for index, folder in enumerate(self.conf.folders_recent):
            self.menu_nav_recent.add_command(
                label=str(folder), command=lambda i=index: self._handle_recent_folder_selection(i)
            )

    def _navigate_to_previous_folder(self):
        if len(self.conf.folders_recent) > 1:
            self._handle_recent_folder_selection(1)

    def _handle_recent_folder_selection(self, folder_index: int):
        self.conf.directory = self.conf.folders_recent[folder_index]
        self.dir_tree.dir_list.try_to_select_directory(self.conf.directory)
        self._update_file_list()

    def _get_files_in_current_directory(self, workflow_file_patterns: List[str]) -> List[Tuple[str, str, str, str]]:
        """
        Returns a list of file information for each file in the currently selected directory
        :param workflow_file_patterns: List of file patterns that this workflow supports, for filtering
        :return: A list of tuples which each contain (file name, file type, file size, modified time) for each file
        """
        file_list = []
        # loop over all non-hidden files to retrieve data for each
        for iter_path in self.conf.directory.glob('*'):
            if iter_path.is_file() and not iter_path.name.startswith('.'):
                base_name = iter_path.name
                # check to see if this file type matches any of the workflow patterns
                for file_type in workflow_file_patterns:
                    if fnmatch(base_name, file_type):
                        break  # breaks out of for loop, allowing processing to continue 3 lines down
                else:  # if we never broke, we will hit this else block
                    continue  # in which case, we want to skip this file, so continue the outer for loop
                file_size = iter_path.stat().st_size
                file_modified_time = iter_path.lstat().st_mtime
                modified_time_string = str(datetime.fromtimestamp(file_modified_time).replace(microsecond=0))
                file_size_string = '{0:12,.0f} KB'.format(file_size / 1024)  # size
                guessed_type = guess_type(base_name)[0]
                file_type_string = "(unknown)" if guessed_type is None else guessed_type
                file_list.append((base_name, file_type_string, file_size_string, modified_time_string))
        # sort the list and return it
        file_list.sort(key=lambda x: x[0])
        return file_list

    def _new_dir_selected(self, selected_path: Path):
        self.previous_selected_directory = self.conf.directory
        self.conf.directory = selected_path
        if len(self.conf.folders_recent) == 0:
            self.conf.folders_recent.appendleft(selected_path)
        elif len(self.conf.folders_recent) > 0 and self.conf.folders_recent[0] != selected_path:
            self.conf.folders_recent.appendleft(selected_path)
        self._rebuild_recent_folder_menu()
        try:
            self._update_status_bar(f"Selected directory: {self.conf.directory}")
            self._update_file_list()
        except Exception as e:  # noqa -- status_bar and things may not exist during initialization, just ignore
            print(str(e))  # log it to the console for fun

    def _update_file_list(self):
        """Update the file listing widget by querying the directory and cache contents, try to reselect current files"""
        # If selected directory hasn't been set up yet then just carry on, this is only happening during app init
        if self._init or not self.conf.directory:
            return

        # If we aren't in a directory, just warn and abort, should not really be possible to get him after init
        if not self.conf.directory.exists() or not self.conf.directory.is_dir():
            self._update_status_bar(f"Bad directory selection: {self.conf.directory}")
            return

        # If we are staying in the same directory, try to select the files that were previously selected
        if self.previous_selected_directory == self.conf.directory:
            previous_selected_files = self.conf.file_selection
        else:
            previous_selected_files = []

        # there should be a cache file there, so get the cached data for the current workflow if it exists
        files_in_current_workflow = {}
        workflow_file_patterns = []
        workflow_columns = []
        if self.workflow_manager.current_workflow:
            cache = CacheFile(self.conf.directory)
            files_in_current_workflow = cache.get_files_for_workflow(
                self.workflow_manager.current_workflow.name
            )
            workflow_file_patterns = self.workflow_manager.current_workflow.file_types
            workflow_columns = self.workflow_manager.current_workflow.columns

        # then get the entire list of files in the current directory to build up the listview items
        # if they happen to match the filename in the workflow cache, then add that info to the row structure
        control_list_rows = []
        files_in_dir = self._get_files_in_current_directory(workflow_file_patterns)
        for file_structure in files_in_dir:
            # always add the columns to the raw list for all files
            file_name = file_structure[0]
            # listview row always includes the filename itself, so start the array with that
            row = [file_name]
            # if it in the cache then the listview row can include additional data
            if file_name in files_in_current_workflow:
                # potentially add a stale column token
                response = self._is_file_stale(file_name)
                if response is None:  # doesn't support stale, so ignore
                    pass
                elif response:  # does support stale and it's true
                    row.append('*')
                else:  # does support stale and it's not
                    row.append('')
                cached_file_info = files_in_current_workflow[file_name]
                if self.workflow_manager.current_workflow.uses_weather:
                    if CacheFile.ParametersKey in cached_file_info:
                        if CacheFile.WeatherFileKey in cached_file_info[CacheFile.ParametersKey]:
                            full_weather_path = cached_file_info[CacheFile.ParametersKey][CacheFile.WeatherFileKey]
                            weather_path_object = Path(full_weather_path)
                            row.append(weather_path_object.name)
                        else:
                            row.append('<no_weather_files>')
                    else:
                        row.append('<no_weather_file>')
                if CacheFile.ResultsKey in cached_file_info:
                    for column in workflow_columns:
                        if column in cached_file_info[CacheFile.ResultsKey]:
                            row.append(cached_file_info[CacheFile.ResultsKey][column])

            # always add the row to the main list
            control_list_rows.append(row)

        # clear all items from the listview and then add them back in
        self.file_list.tree.set_files(control_list_rows)
        if previous_selected_files:
            self.file_list.tree.try_to_reselect(previous_selected_files)
        # self._rebuild_group_menu()

    def _is_file_stale(self, input_file_name: str) -> Optional[bool]:
        """
        Returns a tri-state value trying to characterize whether the current file is "stale."
        It tries to determine this based on the output suffixes of the current workflow.
        If the current workflow does not include any output suffixes, "stale" cannot be determined.
        If '.err' is in the output suffixes, this is the only suffix checked.

        :returns: None if stale cannot be determined, True if it is stale, and False if not.
        """
        if len(self.workflow_manager.current_workflow.output_suffixes) == 0:
            return None  # can't support stale without output suffixes
        full_file_path = self.conf.directory / input_file_name
        if full_file_path.exists():
            input_file_date = full_file_path.lstat().st_mtime
            suffixes = self.workflow_manager.current_workflow.output_suffixes
            if '.err' in suffixes:  # for energyplus workflows just use the err file
                suffixes = ['.err']
            file_name_no_ext = full_file_path.with_suffix('').name
            for suffix in suffixes:
                output_sub_dir = f"EPLaunchRun_{file_name_no_ext}"
                output_file_name = file_name_no_ext + suffix
                tentative_output_file_path = self.conf.directory / output_sub_dir / output_file_name
                if tentative_output_file_path.exists():
                    output_file_date = tentative_output_file_path.lstat().st_mtime
                    if output_file_date < input_file_date:
                        return True
        return False

    def _callback_file_selection_changed(self, selected_file_names: List[str]) -> None:
        """This gets called back by the file listing widget when a selection changes"""
        self.conf.file_selection = selected_file_names
        status = NORMAL if len(self.conf.file_selection) > 0 else DISABLED
        self.button_open_in_text['state'] = status
        # disable the IDF editor button, then possibly enable it
        self.button_idf_editor['state'] = DISABLED
        if system() == 'Windows':
            if self.workflow_manager.current_workflow:
                if self.workflow_manager.current_workflow.is_energyplus:
                    if len(self.conf.file_selection) > 0:
                        self.button_idf_editor['state'] = NORMAL
        # update output buttons too
        self._refresh_output_suffix_buttons_based_on_selection()

    # endregion

    # region weather operations

    def _set_weather_widget_state(self, weather_enabled) -> None:
        self.option_weather_recent.configure(state=weather_enabled)
        self.button_weather_select.configure(state=weather_enabled)

    def _repopulate_recent_weather_list(self, try_to_select: Optional[Path] = None):
        self.option_weather_recent['menu'].delete(0, END)
        for x in self.conf.weathers_recent:
            weather_path = x
            self.option_weather_recent['menu'].add_command(
                label=weather_path.name,
                command=lambda w=weather_path: self._handler_weather_recent_option_changed(w)
            )
        if try_to_select in self.conf.weathers_recent:
            self._handler_weather_recent_option_changed(try_to_select)
        elif len(self.conf.weathers_recent) > 0:
            self._handler_weather_recent_option_changed(self.conf.weathers_recent[0])  # could persist in self.conf

    def _handler_weather_recent_option_changed(self, new_weather_path: Path):
        """This is called when the recent weather option menu changes value"""
        self._tk_var_weather_recent.set(str(new_weather_path.name))
        for selected_file_name in self.conf.file_selection:
            cache = CacheFile(self.conf.directory)
            cache.add_config(
                self.workflow_manager.current_workflow.name,
                selected_file_name,
                {'weather': str(new_weather_path)}
            )
        self._update_file_list()

    def _apply_weather_to_a_file_list(self, list_of_file_paths: List[Path]):
        dialog_weather = TkWeatherDialog(self, list(self.conf.weathers_recent))
        self.wait_window(dialog_weather)
        if dialog_weather.exit_code == TkWeatherDialog.CLOSE_SIGNAL_CANCEL:
            return
        if not dialog_weather.selected_weather_file:
            weather_file_to_use = self.dd_only_string
        else:
            weather_file_to_use = dialog_weather.selected_weather_file
            if weather_file_to_use not in self.conf.weathers_recent:
                self.conf.weathers_recent.appendleft(weather_file_to_use)
                self._repopulate_recent_weather_list(weather_file_to_use)
        for selected_path in list_of_file_paths:
            workflow_directory_cache = CacheFile(selected_path.parent)
            workflow_directory_cache.add_config(
                self.workflow_manager.current_workflow.name,
                selected_path.name,
                {'weather': str(weather_file_to_use)}
            )
        self._update_file_list()

    def _set_weather_for_current_group(self):
        self._apply_weather_to_a_file_list(self.conf.group_locations)

    def _open_weather_dialog(self) -> None:
        current_paths = [self.conf.directory / x for x in self.conf.file_selection]
        self._apply_weather_to_a_file_list(current_paths)

    # endregion

    # region workflow running, tracking, callbacks and handlers

    def _repopulate_workflow_context_menu(self):
        """Clears and repopulates the workflow context menu; tries to reselect the last context saved in config"""
        # save the currently selected workflow context for later
        desired_selected_workflow_context = self.conf.cur_workflow_context
        # get all known contexts from the workflow manager
        all_available_contexts = sorted(list(self.workflow_manager.workflow_contexts))
        # clear the context menu entirely
        self.option_workflow_context['menu'].delete(0, END)
        # add in all the contexts, registering a lambda that will set the tk var and request a workflow list update
        for opt in all_available_contexts:
            self.option_workflow_context['menu'].add_command(
                label=opt, command=lambda c=opt: self._handler_workflow_context_option_changed(c)
            )
        # call the handler method with either the saved context name, or the top one in the list
        if desired_selected_workflow_context in all_available_contexts:
            self._handler_workflow_context_option_changed(desired_selected_workflow_context)
        elif len(all_available_contexts) > 0:
            self._handler_workflow_context_option_changed(all_available_contexts[0])

    def _handler_workflow_context_option_changed(self, new_value: str):
        """This is called when the workflow context option menu changes value"""
        # set the tk var and the config var to the new value
        self._tk_var_workflow_context.set(new_value)
        self.conf.cur_workflow_context = self._tk_var_workflow_context.get()
        # then refresh the workflow list with the new context
        self._repopulate_workflow_instance_menu()

    def _repopulate_workflow_instance_menu(self):
        """Clears and repopulates the workflow instance menu; tries to reselect the last instance saved in config"""
        # save the currently selected workflow instance for later
        desired_selected_workflow_name = self.conf.cur_workflow_name
        # get all known workflows for the current context
        self.available_workflows = self.workflow_manager.workflow_instances(self.conf.cur_workflow_context)
        # clear the option menu entirely
        self.option_workflow_instance['menu'].delete(0, END)
        # add in all the instances, registering a lambda that will set the tk var and refresh the file list as needed
        just_names = []
        energyplus_ip_workflow_name = None
        for w in self.available_workflows:
            workflow_name = w.name
            workflow_name_upper = workflow_name.upper()
            if 'ENERGYPLUS' in workflow_name_upper and 'IP' in workflow_name_upper:
                energyplus_ip_workflow_name = workflow_name
            just_names.append(workflow_name)
            self.option_workflow_instance['menu'].add_command(
                label=workflow_name,
                command=lambda n=workflow_name: self._handler_workflow_instance_option_changed(n)
            )
        # call the handler method with either the saved instance name, or an E+ one, or the top one in the list
        if desired_selected_workflow_name in just_names:
            self._handler_workflow_instance_option_changed(desired_selected_workflow_name)
        elif energyplus_ip_workflow_name:
            self._handler_workflow_instance_option_changed(energyplus_ip_workflow_name)
        else:
            self._handler_workflow_instance_option_changed(just_names[0])

    def _handler_workflow_instance_option_changed(self, new_value: str):
        """This is called when the workflow instance option menu changes value"""
        # store the new instance name in the tk var and in the config
        self._tk_var_workflow_instance.set(new_value)
        self.conf.cur_workflow_name = new_value
        # find the new workflow in the list of currently available workflows
        new_workflow = None
        for w in self.available_workflows:
            if w.name == self._tk_var_workflow_instance.get():
                new_workflow = w
                break
        if new_workflow is None:
            messagebox.showerror("There was an unexpected error updating the workflow list, suggest restarting app.")
        # assign the current workflow instance in the workflow manager
        self.workflow_manager.current_workflow = new_workflow
        # update the weather buttons accordingly depending on if the workflow uses weather inputs
        self._set_weather_widget_state(NORMAL if new_workflow.uses_weather else DISABLED)
        # clear the output menu entirely, and set status conditionally
        self._repopulate_output_suffix_options()
        # now that the workflow has been set, repopulate the file list columns and the file list itself
        self._repopulate_control_list_columns()
        self._update_file_list()

    def _repopulate_output_suffix_options(self):
        sorted_suffixes = sorted(self.workflow_manager.current_workflow.output_suffixes)
        combobox_output_enabled = 'readonly' if len(sorted_suffixes) > 0 else 'disabled'
        output_enabled = NORMAL if len(sorted_suffixes) > 0 else DISABLED
        self.option_workflow_outputs.configure(state=combobox_output_enabled)
        self.button_open_output_file.configure(state=output_enabled)

        # rebuild the option menu if applicable
        current_selection = self._tk_var_output_suffix.get()
        if output_enabled == NORMAL:
            self.option_workflow_outputs['values'] = sorted_suffixes
            if current_selection not in sorted_suffixes:
                self._tk_var_output_suffix.set(sorted_suffixes[0])
        else:
            self.option_workflow_outputs['values'] = []
            self._tk_var_output_suffix.set('')
        self.option_workflow_outputs.selection_clear()
        suffixes = self.workflow_manager.current_workflow.output_suffixes
        self.button_open_output_1.configure(state=NORMAL if len(suffixes) > 0 else DISABLED)
        self._tk_var_output_1.set(suffixes[0] if len(suffixes) > 0 else '--')
        self.button_open_output_2.configure(state=NORMAL if len(suffixes) > 1 else DISABLED)
        self._tk_var_output_2.set(suffixes[1] if len(suffixes) > 1 else '--')
        self.button_open_output_3.configure(state=NORMAL if len(suffixes) > 2 else DISABLED)
        self._tk_var_output_3.set(suffixes[2] if len(suffixes) > 2 else '--')
        self.button_open_output_4.configure(state=NORMAL if len(suffixes) > 3 else DISABLED)
        self._tk_var_output_4.set(suffixes[3] if len(suffixes) > 3 else '--')
        self.button_open_output_5.configure(state=NORMAL if len(suffixes) > 4 else DISABLED)
        self._tk_var_output_5.set(suffixes[4] if len(suffixes) > 4 else '--')
        self._refresh_output_suffix_buttons_based_on_selection()

    def _suffixed_paths_exist(self, original_path: Path, new_suffix: str) -> Tuple[Optional[Path], Optional[Path]]:
        filename_no_ext = original_path.with_suffix('').name
        new_path = self.conf.directory / (filename_no_ext + new_suffix)
        eplus_sub_dir = f"EPLaunchRun_{filename_no_ext}"
        eplus_specific_output_path = self.conf.directory / eplus_sub_dir / f"{filename_no_ext}{new_suffix}"
        primary_path_to_return = new_path if new_path.exists() else None
        eplus_path_to_return = eplus_specific_output_path if eplus_specific_output_path.exists() else None
        return primary_path_to_return, eplus_path_to_return

    def _refresh_single_output_suffix_button(self, tk_var: Variable, button: Button):
        if tk_var.get() == '--':
            return
        else:
            suffix_to_open = tk_var.get()
            all_files_have_this_suffix = True
            for f in self.conf.file_selection:
                original_path = self.conf.directory / f
                new_path, eplus_path = self._suffixed_paths_exist(original_path, suffix_to_open)
                if (new_path and new_path.exists()) or (eplus_path and eplus_path.exists()):
                    pass  # good
                else:
                    all_files_have_this_suffix = False
                    break
            button.configure(state=NORMAL if all_files_have_this_suffix else DISABLED)

    def _refresh_output_suffix_buttons_based_on_selection(self):
        self._refresh_single_output_suffix_button(self._tk_var_output_1, self.button_open_output_1)
        self._refresh_single_output_suffix_button(self._tk_var_output_2, self.button_open_output_2)
        self._refresh_single_output_suffix_button(self._tk_var_output_3, self.button_open_output_3)
        self._refresh_single_output_suffix_button(self._tk_var_output_4, self.button_open_output_4)
        self._refresh_single_output_suffix_button(self._tk_var_output_5, self.button_open_output_5)

    def _open_workflow_dir_dialog(self):
        if len(self.workflow_manager.threads) > 0:
            messagebox.showerror("Workflow Selection Error", 'Cannot change workflows while threads are running')
            return  # will leave the window open
        # refresh the list of workflows auto-found on the machine
        self.workflow_manager.auto_find_workflow_directories()
        # pass the newly found folders in, as well as the current active list of workflow directories
        auto_find_workflows = self.workflow_manager.auto_found_workflow_dirs
        current_workflows = self.workflow_manager.workflow_directories
        wf_dialog = TkWorkflowsDialog(self, list(current_workflows), list(auto_find_workflows))
        self.wait_window(wf_dialog)
        if wf_dialog.exit_code == TkWorkflowsDialog.CLOSE_SIGNAL_CANCEL:
            return
        self.workflow_manager.workflow_directories = wf_dialog.list_of_directories
        self.workflow_manager.instantiate_all_workflows()
        if len(self.workflow_manager.warnings) > 0:
            self.generic_dialogs.display(
                self, 'Workflow Processing Errors', '\n'.join(self.workflow_manager.warnings)
            )
        self.conf.workflow_directories = self.workflow_manager.workflow_directories
        self._repopulate_workflow_context_menu()
        self._repopulate_workflow_instance_menu()

    def _get_or_set_weather_for_file(self,
                                     cur_workflow: Optional[BaseEPLaunchWorkflow1],
                                     directory: Path,
                                     selected_file: str,
                                     backup_weather_file_to_use: str,
                                     ) -> Tuple[bool, str, str]:
        """
        cur_workflow is "None" at init time, so it is marked technically optional
        """
        weather_file_to_use: Optional[str] = None
        # if we find this file in the current folder cache
        cache = CacheFile(directory)
        files_in_current_workflow = cache.get_files_for_workflow(cur_workflow.name)
        if selected_file in files_in_current_workflow:
            # and if we find a weather file key for that file, for that workflow, just use the weather file
            cached_file_info = files_in_current_workflow[selected_file]
            if CacheFile.ParametersKey in cached_file_info:
                if CacheFile.WeatherFileKey in cached_file_info[CacheFile.ParametersKey]:
                    weather_file_to_use = cached_file_info[CacheFile.ParametersKey][CacheFile.WeatherFileKey]
        # if we didn't find a weather file anywhere, we need to ask
        if not weather_file_to_use:
            # since we could be asking for many files, we don't want to ask for each
            # get a backup file to use only once, and apply that to any missing ones as needed
            if backup_weather_file_to_use:
                weather_file_to_use = backup_weather_file_to_use
            else:
                # if we need weather, didn't find one in the cache, and didn't have a backup, ask for one now
                recent_files = list(self.conf.weathers_recent)
                w = TkWeatherDialog(self, recent_files, "*At least one file is missing a weather configuration*")
                self.wait_window(w)
                if w.exit_code == TkWeatherDialog.CLOSE_SIGNAL_CANCEL:
                    return False, '', ''
                else:  # a valid response was encountered
                    if not w.selected_weather_file:
                        weather_file_to_use = self.dd_only_string
                    else:
                        weather_file_to_use = str(w.selected_weather_file)
                        if w.selected_weather_file not in self.conf.weathers_recent:
                            self.conf.weathers_recent.appendleft(w.selected_weather_file)
            # save the current selected one as the new backup for subsequent files
            backup_weather_file_to_use = weather_file_to_use
        # add the weather configuration to the cache regardless of how it was retrieved, and update file listing
        cache.add_config(cur_workflow.name, selected_file, {'weather': weather_file_to_use})
        return True, weather_file_to_use, backup_weather_file_to_use

    def _run_workflow_on_selection(self) -> None:
        cur_files = self.conf.file_selection
        file_paths = [self.conf.directory / f for f in cur_files]
        self._run_workflow_on_list_of_file_paths(file_paths)

    def _run_workflow_on_group(self) -> None:
        file_paths = self.conf.group_locations
        self._run_workflow_on_list_of_file_paths(file_paths)

    def _run_workflow_on_list_of_file_paths(self, file_paths: List[Path]):

        # error out early
        if self.conf.directory and file_paths and self.workflow_manager.current_workflow:
            pass
        else:
            messagebox.showerror(title='Selection', message='ERROR: Select a workflow, directory and a file')
            return

        already_running_instances = []
        cur_workflow = self.workflow_manager.current_workflow
        w_name = cur_workflow.name
        backup_weather_file_to_use: Optional[str] = None

        # loop over all the selected files and try to run the current workflow on each of them
        for path in file_paths:
            # if the current workflow uses weather, we need to determine what to pass into it
            weather_file_to_use: Optional[str] = None
            if cur_workflow.uses_weather:
                success, weather_file_to_use, backup_weather_file_to_use = self._get_or_set_weather_for_file(
                    cur_workflow, path.parent, path.name, backup_weather_file_to_use
                )
                if not success:
                    return  # TODO: probably shouldn't just return blindly here
            this_wea = '' if weather_file_to_use == self.dd_only_string else weather_file_to_use

            # now we actually need to run this current file
            this_file_good_to_go = True
            # first check to see if the current file/dir/workflow combination is already running
            for thread_id, t in self.workflow_manager.threads.items():
                t_path = Path(t.run_directory) / t.file_name
                if t_path == path and t.workflow_instance.name() == w_name:
                    # save it in a list and wait to issue error message until the end
                    already_running_instances.append(f"{w_name}: {path}")
                    this_file_good_to_go = False
                    break
            # if this file isn't good, just continue the file loop
            if not this_file_good_to_go:
                continue
            # otherwise, continue to instantiate a workflow thread instance to let it run
            new_uuid = str(uuid4())
            new_instance = cur_workflow.workflow_class()
            new_instance.register_standard_output_callback(new_uuid, self._callback_workflow_stdout)
            self.workflow_manager.threads[new_uuid] = WorkflowThread(
                new_uuid, new_instance, path.parent, path.name,
                {'weather': this_wea, 'workflow location': cur_workflow.workflow_directory},
                self._callback_workflow_done
            )
            self.output_dialogs[new_uuid] = self._create_output_dialog(new_uuid)
            self.output_dialogs[new_uuid].add_output("*** STARTING WORKFLOW ***")

        # emit an error dialog if needed
        if len(already_running_instances) > 0:
            out = '\n'.join(already_running_instances)
            messagebox.showerror("Run Error", f"Some configurations were already running, and were skipped: {out}")

        self._update_file_list()  # do one update once all threads are spawned to update for weather selection, etc.
        self._update_status_bar("Currently %s processes running" % len(self.workflow_manager.threads))

    def _create_output_dialog(self, workflow_id: str) -> TkOutputDialog:
        """Generates an output dialog with the specified ID, updating dialog counter and setting dialog position"""
        max_dialog_vertical_increments = 5.0
        self.dialog_counter += 1
        if self.dialog_counter == max_dialog_vertical_increments:
            self.dialog_counter = 1

        this_workflow = self.workflow_manager.threads[workflow_id]

        x_right_edge = self.winfo_x() + self.winfo_width()
        y_top = self.winfo_y()
        vertical_increment = int(self.winfo_height() / max_dialog_vertical_increments / 2.0)
        this_x = x_right_edge + 5
        this_y = y_top + vertical_increment * (self.dialog_counter - 1)
        thread_data = self.workflow_manager.threads[workflow_id]
        config_string = dumps(
            {
                'workflow_name:': thread_data.workflow_instance.name(),
                'workflow_dir': str(thread_data.workflow_directory),
                'file_name:': thread_data.file_name,
                'run_directory:': str(thread_data.run_directory),
            },
            indent=2
        )
        return TkOutputDialog(
            self,
            workflow_id,
            this_workflow.workflow_instance.name(),
            config_string,
            this_x,
            this_y
        )

    def _callback_workflow_done(self, workflow_response: EPLaunchWorkflowResponse1) -> None:
        self._gui_queue.put(lambda: self._handler_workflow_done(workflow_response))

    def _handler_workflow_done(self, workflow_response: EPLaunchWorkflowResponse1) -> None:
        try:
            if workflow_response.success:
                status_message = 'Successfully completed a workflow: ' + workflow_response.message
                try:
                    data_from_workflow = workflow_response.column_data
                    workflow_working_directory = self.workflow_manager.threads[workflow_response.id].run_directory
                    workflow_directory_cache = CacheFile(workflow_working_directory)
                    workflow_directory_cache.add_result(
                        self.workflow_manager.threads[workflow_response.id].workflow_instance.name(),
                        self.workflow_manager.threads[workflow_response.id].file_name,
                        data_from_workflow
                    )
                    if self.conf.directory == workflow_working_directory:
                        # only update file lists if we are still in that directory
                        self._update_file_list()
                except EPLaunchFileException:
                    pass
                if not self.conf.keep_dialog_open:
                    self.output_dialogs[workflow_response.id].close()
            else:
                status_message = 'Workflow failed: ' + workflow_response.message
                self.output_dialogs[workflow_response.id].add_output('Workflow FAILED: ' + workflow_response.message)
                self.output_dialogs[workflow_response.id].add_output("*** WORKFLOW FINISHED ***")
        except Exception as e:  # noqa -- there is *no* telling what all exceptions could occur inside a workflow
            print(e)
            status_message = 'Workflow response was invalid'
        self._update_status_bar(status_message)
        try:
            del self.workflow_manager.threads[workflow_response.id]
        except Exception as e:
            print(e)
        self._update_status_bar("Currently %s processes running" % len(self.workflow_manager.threads))

    def _callback_workflow_stdout(self, workflow_id: str, message: str) -> None:
        self._gui_queue.put(lambda: self._handler_workflow_stdout(workflow_id, message))

    def _handler_workflow_stdout(self, workflow_id: str, message: str) -> None:
        if workflow_id in self.output_dialogs:
            if self.output_dialogs[workflow_id].winfo_exists():
                self.output_dialogs[workflow_id].add_output(message)

    # endregion

    # region misc dialog and external tool launchers

    def _open_viewers_dialog(self) -> None:
        if self.workflow_manager.current_workflow:
            output_suffixes = self.workflow_manager.current_workflow.output_suffixes
        else:
            output_suffixes = []
        viewer_dialog = TkViewerDialog(self, output_suffixes, self.conf.viewer_overrides)
        self.wait_window(viewer_dialog)
        if viewer_dialog.exit_code == TkViewerDialog.CLOSE_SIGNAL_CANCEL:
            return
        self.conf.viewer_overrides = {
            **viewer_dialog.extension_to_viewer, **self.conf.viewer_overrides
        }

    def _open_welcome(self):
        if self.conf.welcome_shown and VERSION == self.conf.latest_welcome_shown:
            return
        message = """
EnergyPlus-Launch has been around for many years as a part of the EnergyPlus distribution.
Starting with the 3.0 release, it has changed drastically, completely redesigned and rewritten.
For full documentation or a quick start guide, click the "Open Docs" button below.
This dialog will only be shown once, but documentation is available in the Help menu."""
        self.generic_dialogs.display_with_alt_button(
            self, 'Welcome to EnergyPlus-Launch ' + VERSION, message, 'Open Documentation', self._open_documentation
        )
        self.conf.welcome_shown = True
        self.conf.latest_welcome_shown = VERSION

    def _open_shortcuts_dialog(self, *_) -> None:
        shortcuts = EPLaunchWindow._list_keyboard_shortcuts()
        text = '\n'.join([f"{r[0]}:\t{r[1]}" for r in shortcuts])
        self.generic_dialogs.display(self, "Available Keyboard Shortcuts", text)

    def _open_about(self) -> None:
        text = """
EnergyPlus-Launch

EnergyPlus-Launch is a graphical workflow manager for EnergyPlus.
Originally written in VB6 and released with essentially every version of EnergyPlus,
it is now a cross platform Python tool.
Users are encouraged to leverage EnergyPlus-Launch and write their own workflow scripts
to enable new capabilities.

Version %s
Copyright (c) 2023, Alliance for Sustainable Energy, LLC  and GARD Analytics, Inc

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
 list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

* The name of the copyright holder(s), any contributors, the United States
Government, the United States Department of Energy, or any of their employees
may not be used to endorse or promote products derived from this software
without specific prior written permission from the respective party.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDER(S) AND ANY CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER(S), ANY ,
CONTRIBUTORS THE UNITED STATES GOVERNMENT, OR THE UNITED STATES DEPARTMENT
OF ENERGY, NOR ANY OF THEIR EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
        """ % VERSION
        self.generic_dialogs.display(self, "About EnergyPlus-Launch", text)

    def _open_output_file_generic(self, suffix_to_open: str):
        if len(self.conf.file_selection) > 5:
            message = """More than 5 files were selected when choosing to open workflow outputs.
        This could generate many output windows and is usually not intentional.  Select up to 5 input files."""
            messagebox.showerror("File Selection Issue", message)
            return
        for f in self.conf.file_selection:
            original_path = self.conf.directory / f
            filename_no_ext = original_path.with_suffix('').name
            new_path = self.conf.directory / (filename_no_ext + suffix_to_open)
            new_path_str = str(new_path)
            eplus_sub_dir = f"EPLaunchRun_{filename_no_ext}"
            eplus_specific_output_path = self.conf.directory / eplus_sub_dir / f"{filename_no_ext}{suffix_to_open}"
            eplus_specific_output_file = str(eplus_specific_output_path)
            if new_path.exists():
                if suffix_to_open in self.conf.viewer_overrides:
                    if self.conf.viewer_overrides[suffix_to_open] is not None:
                        # then the viewer override was found, and not None, so use it
                        Popen([str(self.conf.viewer_overrides[suffix_to_open]), new_path_str])
                        continue
                # if we make it this far, we didn't open it with a custom viewer, try to use the default
                self._open_file_or_dir_with_default(new_path)
            elif eplus_specific_output_path.exists():
                if suffix_to_open in self.conf.viewer_overrides:
                    if self.conf.viewer_overrides[suffix_to_open] is not None:
                        # then the viewer override was found, and not None, so use it
                        Popen([str(self.conf.viewer_overrides[suffix_to_open]), eplus_specific_output_file])
                        continue
                # if we make it this far, we didn't open it with a custom viewer, try to use the default
                self._open_file_or_dir_with_default(eplus_specific_output_path)
            else:
                message = """At least some of the output files were not found.
Make sure that the workflow has run on the selected input file(s), and that the runs
actually generated the requested outputs.  Any found output files are being opened now."""
                messagebox.showwarning("Missing Output Files", message)

    def _open_output_file(self):
        self._open_output_file_generic(self._tk_var_output_suffix.get())

    def _open_output_1(self):
        self._open_output_file_generic(self._tk_var_output_1.get())

    def _open_output_2(self):
        self._open_output_file_generic(self._tk_var_output_2.get())

    def _open_output_3(self):
        self._open_output_file_generic(self._tk_var_output_3.get())

    def _open_output_4(self):
        self._open_output_file_generic(self._tk_var_output_4.get())

    def _open_output_5(self):
        self._open_output_file_generic(self._tk_var_output_5.get())

    @staticmethod
    def _open_file_or_dir_with_default(full_path_obj: Path, text_file_override: bool = False) -> None:
        full_path = str(full_path_obj)
        if system() == 'Windows':
            from os import startfile
            startfile(full_path)
        elif system() == 'Linux':
            Popen(['xdg-open', full_path])
            # could try to find the default text/plain editor and open using gtk-launch
            # if text_file_override:
            #     desktop_file = check_output(['xdg-mime', 'query', 'default', 'text/plain'], shell=False)
            #     if '.desktop' in desktop_file.decode('utf-8'):
            #         Popen(['gtk-launch', desktop_file, full_path])
            #     else:  # try falling back to the default application
            #         Popen(['xdg-open', full_path])
            # else:
            #     Popen(['xdg-open', full_path])
        else:  # assuming Mac
            if text_file_override:  # will open in the default text editor specifically
                Popen(['open', '-e', full_path])
            else:
                # alright, because we are trying to launch using `open`, we are severely limited
                # in how much information we can get out of the running child process
                # we can get info about the actual `open` process, but not the process that `open` spawns
                # it is launched as a new process whose parent is simply launchd (root)
                # I tried checking the running status of the open process, child subprocesses, and even
                # some hacky tricks to check the incremental PID after the open call, and nothing worked
                # reliably.  So I'm just going to pick a few known file extensions that I feel comfortable
                # dispatching to a default program (htm, html, csv), and then send everything else to text
                if full_path.endswith('htm') or full_path.endswith('html') or full_path.endswith('csv'):
                    Popen(['open', full_path])
                else:
                    Popen(['open', '-e', full_path])

    def _open_text_editor(self) -> None:
        for file_name in self.conf.file_selection:
            full_path_str = self.conf.directory / file_name
            if 'txt' in self.conf.viewer_overrides and self.conf.viewer_overrides['txt']:
                text_editor_binary = str(self.conf.viewer_overrides['txt'])
                Popen([text_editor_binary, full_path_str])
                return
            if system() == 'Windows':
                potential_program = self._find_default_text_editor_on_windows()
                if potential_program != '':
                    Popen([potential_program, full_path_str])
                    return
            # if neither of those work, just open with the default editor
            self._open_file_or_dir_with_default(full_path_str)

    def _open_idf_editor(self) -> None:
        # should only be able to get here if the current configuration/workflow/selection/platform allows it
        workflow_context_dir = self.workflow_manager.current_workflow.workflow_directory
        eplus_root_dir = workflow_context_dir.parent
        idf_editor_dir = eplus_root_dir / 'PreProcess' / 'IDFEditor'
        idf_editor_binary = idf_editor_dir / 'IDFEditor.exe'
        if len(self.conf.file_selection) > 3:
            messagebox.showerror("File Selection Issue", "Select up to 3 files to open with IDF Editor.")
            return
        any_skipped = False
        for f in self.conf.file_selection:
            if f.endswith('idf') or f.endswith('imf'):
                Popen([idf_editor_binary, self.conf.directory / f])
            else:
                any_skipped = True
        if any_skipped:
            messagebox.showerror("File Type Issue", "At least one file was skipped because it was not an IDF/IMF")

    @staticmethod
    def _find_default_text_editor_on_windows() -> str:
        from winreg import OpenKey, QueryValueEx, HKEY_CLASSES_ROOT
        try:
            key = OpenKey(HKEY_CLASSES_ROOT, '.txt')
            default_value, _ = QueryValueEx(key, '')
            key.Close()

            key = OpenKey(HKEY_CLASSES_ROOT, default_value + r'\shell\open\command')
            command, _ = QueryValueEx(key, '')
            key.Close()

            # resolve the SystemRoot item
            command = command.replace('%SystemRoot%', getenv('SystemRoot'))

            # then just take the program name by splitting any trailing % placeholders
            return command.split('%')[0].strip()
        except WindowsError:
            return ''

    def _open_file_browser(self) -> None:
        self._open_file_or_dir_with_default(self.conf.directory)

    @staticmethod
    def _open_documentation() -> None:
        open_web(DOCS_URL)

    # endregion


if __name__ == "__main__":
    window = EPLaunchWindow()
    window.run()
