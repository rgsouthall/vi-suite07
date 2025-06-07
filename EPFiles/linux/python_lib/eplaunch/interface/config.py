from collections import deque
from configparser import ConfigParser
from json import loads, dumps
from pathlib import Path
from re import search
from typing import Any, Dict, List, Optional


class ConfigManager:
    config_file_name = '.EP-Launch.json'

    def __init__(self):
        # specify configuration settings here
        self.keep_dialog_open: bool = True
        self.cur_workflow_name: str = ''
        self.cur_workflow_context: str = ''
        self.directory: Path = Path.home()  # default to the home directory
        self.file_selection: List[str] = []
        self.welcome_shown: bool = False
        self.latest_welcome_shown: str = ''
        self.height: int = -1
        self.width: int = -1
        self.x: int = -1
        self.y: int = -1
        self.workflow_directories: List[Path] = []
        self.folders_recent: deque[Optional[Path]] = deque(maxlen=5)
        self.folders_favorite: List[Path] = []
        self.weathers_recent: deque[Optional[Path]] = deque(maxlen=5)
        # self.weathers_favorite: List[Path] = []
        self.group_locations: List[Path] = []
        # self.groups_recent: List[str] = []
        # self.groups_favorite: List[str] = []
        self.viewer_overrides: Dict[str, Path] = {}
        self.dir_file_paned_window_sash_position: int = 400
        self.list_group_paned_window_sash_position: int = 500

    def load(self, called_from_ep_cli: bool):
        # load the config from the file on disk
        config_file_path = Path.home() / ConfigManager.config_file_name
        if not config_file_path.exists():
            pass  # just use the default values
        else:
            contents = config_file_path.read_text()
            if contents == '':
                pass  # just use the default values, something is wrong
            else:
                if contents[0] == '[':
                    # we have an INI style file, read a few things
                    parser = ConfigParser()
                    parser.read(config_file_path)

                    def parse(x1: str, x2: str, default: Any):
                        if x1 in parser and x2 in parser[x1]:
                            return parser[x1][x2]
                        return default
                    self.keep_dialog_open = parse('ActiveWindow', 'KeepDialogOpen', self.keep_dialog_open)
                    self.cur_workflow_name = parse('ActiveWindow', 'SelectedWorkflow', self.cur_workflow_name)
                    self.cur_workflow_context = parse('ActiveWindow', 'CurrentContext', self.cur_workflow_context)
                    self.directory = parse('ActiveWindow', 'CurrentDirectory', self.directory)
                    if isinstance(self.directory, str):
                        self.directory = Path(self.directory)
                    self.welcome_shown = parse('ActiveWindow', 'WelcomeAlreadyShown', self.welcome_shown)
                    self.latest_welcome_shown = parse(
                        'ActiveWindow', 'LatestWelcomeVersionShown', self.latest_welcome_shown
                    )
                    self.height = parse('ActiveWindow', 'height', self.height)
                    self.width = parse('ActiveWindow', 'width', self.width)
                    self.x = parse('ActiveWindow', 'x', self.x)
                    self.y = parse('ActiveWindow', 'y', self.y)
                    # These are more complex, I don't intend on reading these in
                    # self.workflow_directories: List[str] = []  /WorkflowDirectories/Path-##=/foo/bar
                    # self.folders_recent: List[str] = []  /FolderMenu/Recent/Path-##=/foo/bar
                    # self.folders_favorite: List[str] = []  /FolderMenu/Favorite/Path-##=/foo/bar
                    # self.weathers_recent: List[str] = []  /WeatherMenu/Recent/Path-##=/foo/bar
                    # self.weathers_favorite: List[str] = []  /WeatherMenu/Favorite/Path-##=/foo/bar
                    # self.groups_recent: List[str] = []  /GroupMenu/Recent/Path-##=/foo/bar
                    # self.groups_favorite: List[str] = []  /GroupMenu/Favorite/Path-##=/foo/bar
                    # this is different: /ViewerOverrides/Path-##=/foo/bar, /ViewerOverrides/Ext-##=.py # verify ext
                    # self.viewer_overrides: Dict[str, str] = {}

                elif contents[0] == '{':
                    # we have a JSON style file, read it
                    config = loads(contents)
                    self.keep_dialog_open = config.get('KeepDialogOpen', self.keep_dialog_open)
                    self.cur_workflow_name = config.get('SelectedWorkflow', self.cur_workflow_name)
                    self.cur_workflow_context = config.get('CurrentContext', self.cur_workflow_context)
                    self.directory = config.get('CurrentDirectory', self.directory)
                    if isinstance(self.directory, str):
                        self.directory = Path(self.directory)
                    # self.cur_filename = config.get('CurrentFileName', self.cur_filename)
                    self.file_selection = config.get('CurrentFileSelection', self.file_selection)
                    self.welcome_shown = config.get('WelcomeAlreadyShown', self.welcome_shown)
                    self.latest_welcome_shown = config.get('LatestWelcomeVersionShown', self.latest_welcome_shown)
                    self.height = config.get('height', self.height)
                    self.width = config.get('width', self.width)
                    self.x = config.get('x', self.x)
                    self.y = config.get('y', self.y)
                    self.workflow_directories = [
                        Path(p) for p in config.get('WorkflowDirectories', self.workflow_directories)
                    ]
                    if called_from_ep_cli:  # if we are called from E+ CLI, set the E+ dir directly from this file path
                        this_file = Path(__file__).resolve()
                        interface_dir = this_file.parent
                        ep_launch_package_dir = interface_dir.parent
                        python_lib_dir = ep_launch_package_dir.parent
                        eplus_install_dir = python_lib_dir.parent
                        workflow_path = eplus_install_dir / 'workflows'
                        eplus_workflow = workflow_path / 'energyplus.py'
                        context_pattern = r"EnergyPlus-(\d+\.\d+\.\d+)-([0-9a-f]+)"
                        context = search(context_pattern, eplus_workflow.read_text())
                        version = context.group(1)
                        sha = context.group(2)
                        self.cur_workflow_name = f"EnergyPlus-{version} SI"
                        self.cur_workflow_context = f"EnergyPlus-{version}-{sha}"
                        if workflow_path not in self.workflow_directories:
                            self.workflow_directories.append(workflow_path)
                    recent_folders = config.get('RecentFolders', self.folders_recent)
                    if recent_folders is not None:
                        for string_path in recent_folders:
                            self.folders_recent.appendleft(Path(string_path))
                    self.folders_favorite = [Path(p) for p in config.get('FavoriteFolders', self.folders_favorite)]
                    recent_weather = config.get('RecentWeather', self.weathers_recent)
                    if recent_weather is not None:
                        for string_path in recent_weather:
                            self.weathers_recent.appendleft(Path(string_path))
                    # self.weathers_favorite = [Path(p) for p in config.get('FavoriteWeather', self.weathers_favorite)]
                    # self.groups_recent = config.get('RecentGroup', self.groups_recent)
                    # self.groups_favorite = config.get('FavoriteGroup', self.groups_favorite)
                    temp_group_locations = config.get('GroupLocations', self.group_locations)
                    for t in temp_group_locations:
                        self.group_locations.append(Path(t))
                    for k, v in config.get('ViewerOverrides', self.viewer_overrides).items():
                        self.viewer_overrides[k] = Path(v) if v else None
                    self.dir_file_paned_window_sash_position = config.get(
                        'DirFilePanedWindowSash', self.dir_file_paned_window_sash_position
                    )
                    self.list_group_paned_window_sash_position = config.get(
                        'ListGroupPanedWindowSash', self.list_group_paned_window_sash_position
                    )
                else:
                    pass  # Bad config saved file format?  Indicates a crash?
                # fix up the current selected directory to initialize in case it doesn't exist (anymore)
                if not self.directory:
                    self.directory = Path.home()
                elif not self.directory.exists():
                    self.directory = Path.home()

    def save(self):
        # save the config file to disk
        output_dict = {
            'KeepDialogOpen': self.keep_dialog_open,
            'SelectedWorkflow': self.cur_workflow_name,
            'CurrentContext': self.cur_workflow_context,
            'CurrentDirectory': str(self.directory),
            'CurrentFileSelection': self.file_selection,
            'WelcomeAlreadyShown': self.welcome_shown,
            'LatestWelcomeVersionShown': self.latest_welcome_shown,
            'height': self.height,
            'width': self.width,
            'x': self.x,
            'y': self.y,
            'WorkflowDirectories': [str(p) for p in self.workflow_directories],
            'RecentFolders': [str(p) for p in self.folders_recent if p is not None],
            'FavoriteFolders': [str(p) for p in self.folders_favorite],
            'RecentWeather': [str(p) for p in self.weathers_recent if p is not None],
            # 'FavoriteWeather': [str(p) for p in self.weathers_favorite],
            # 'RecentGroup': self.groups_recent,
            # 'FavoriteGroup': self.groups_favorite,
            'GroupLocations': [str(t) for t in self.group_locations],
            'ViewerOverrides': {k: None if v is None else str(v) for k, v in self.viewer_overrides.items()},
            'DirFilePanedWindowSash': self.dir_file_paned_window_sash_position,
            'ListGroupPanedWindowSash': self.list_group_paned_window_sash_position,
        }
        config_file_path = Path.home() / ConfigManager.config_file_name
        config_file_path.write_text(dumps(output_dict, indent=2))
