from shutil import copyfile
from os import chmod, path, stat, symlink, fdopen, unlink
from pathlib import Path
from platform import system
from site import USER_BASE, getusersitepackages
from subprocess import check_call
from sys import argv
from sysconfig import get_path
from tempfile import mkstemp
from typing import List


class EntryPoint:
    """
    This class represents a single executable that is packaged up by setuptools using the entry_points method.
    The entry points may live in either the console_scripts or the gui_scripts subsections.
    Either way, when the package is pip installed, it results in a binary blob in a bin/ or Scripts/ directory.
    This class helps add user-friendly functionality to these binary blobs.
    """

    def __init__(self, package_source_dir: str, executable_name: str, nice_name: str, description: str, wm_class: str,
                 test_mode: bool = False):
        """
        Construct a single script instance, using arguments.

        :param package_source_dir: The name of the source dir containing the package code
        :param executable_name: The name of a resulting executable file (without the exe extension)
        :param nice_name: The "nice" name to use to reference the tool in links and menus
        :param description: The description to show up in Linux .desktop files
        :param wm_class: The wm-class for the Tk GUI (assigned in Tk(className='energyplus_regression_runner'))
        :param test_mode: If True, this will ensure that no files are actually written, and may respond to env vars
                          or some other test mocking later
        """
        # we are going to assume that this script file has been installed as part of a pip install plan-tools step
        self.package_root = Path(__file__).parent.parent / package_source_dir
        # this will ultimately become the name field in .desktop files, or the desktop link name,
        # so don't include a .lnk extension
        self.pretty_link_name = nice_name
        # these are used by the .desktop specification to help add a richer experience
        self.description = description
        self.wm_class = wm_class
        # this is the executable that will be linked; we can append the .exe to the binary name if needed
        self.installed_binary_name = executable_name
        if system() == 'Windows':
            self.installed_binary_name += '.exe'
        self.test_mode = test_mode
        self.desktop_file_data_check = ""  # used for unit test checks, will vary per platform

    def run(self) -> int:
        """
        This function will perform all relevant actions to set up this installed entry point on the system.
        Currently, this simply sets up a desktop icon on Windows or installs a .desktop file on Linux.
        Moving forward, this will create an .app bundle on Mac and create start menu entries on Windows.
        A future addition will also allow options to be set to only do certain actions.

        :return: zero for success, nonzero otherwise
        """
        return self.add_desktop_icon()  # eventually we'll accumulate any errors and return intelligently

    # def uninstall(self) -> int:
    #     """
    #     This function is not yet implemented, but could ultimately do the same set-up as the installation, and instead
    #     try to remove any artifacts that do exist.  Right now, since this is just dropping a desktop link, I'm not
    #     worrying about it
    #
    #     :return: zero for success, nonzero otherwise
    #     """
    #     return 0

    def add_desktop_icon(self) -> int:
        """
        Worker function to add the desktop icon for the package specified by the constructor arguments.
        This function is just a list of steps, which capture the required actions in a platform-agnostic way

        :return: Zero if successful, nonzero otherwise
        """
        scripts_dir = self.get_pip_entry_point_exe_dir()
        target_exe = scripts_dir / self.installed_binary_name
        icon_file = self.get_package_icon_for_link()
        desktop_file = self.get_path_to_link_file()
        self.write_desktop_file(desktop_file, target_exe, scripts_dir, icon_file)
        return 0

    def write_desktop_file(self, desktop_file: Path, target_exe: Path, scripts_dir: Path, icon_file: Path):
        if system() == 'Windows':
            def escape_path(path_: Path) -> str:
                return str(path_).replace('\\', '/')

            icon_file_string: str = escape_path(icon_file) if icon_file.exists() else ''
            self.desktop_file_data_check = {'exe': target_exe, 'cwd': scripts_dir, 'icon': icon_file_string}
            if self.test_mode:
                return
            else:  # pragma: no cover
                shortcut_path = escape_path(desktop_file)
                target = escape_path(target_exe)
                working_dir = escape_path(scripts_dir)
                js_content = f'''var sh = WScript.CreateObject("WScript.Shell");
var shortcut = sh.CreateShortcut("{shortcut_path}");
shortcut.TargetPath = "{target}";
shortcut.WorkingDirectory = "{working_dir}";
shortcut.IconLocation = "{icon_file_string}";
shortcut.Save();'''
                script_fd, script_path = mkstemp('.js')
                try:
                    with fdopen(script_fd, 'w') as f:
                        f.write(js_content)
                    check_call([R'wscript.exe', script_path])
                finally:
                    unlink(script_path)
        elif system() == 'Linux':
            desktop_file_contents = f"""[Desktop Entry]
Name={self.pretty_link_name}
Comment={self.description}
Exec={target_exe}
Icon={str(icon_file) if icon_file.exists() else ''}
Type=Application
Path={scripts_dir}
Terminal=false
StartupWMClass={self.wm_class}"""
            self.desktop_file_data_check = {'contents': desktop_file_contents}
            if self.test_mode:
                return
            else:  # pragma: no cover
                # maybe someday we could actually test this in CI
                with open(desktop_file, 'w') as f:
                    f.write(desktop_file_contents)
                mode = stat(desktop_file).st_mode
                mode |= (mode & 0o444) >> 2  # copy R bits to X
                chmod(desktop_file, mode)  # make it executable
        else:  # assuming 'Darwin'
            # the link is actually a folder
            contents_dir = desktop_file / 'Contents'
            mac_os_dir = contents_dir / 'MacOS'
            resource_dir = contents_dir / 'Resources'
            new_exe_path = mac_os_dir / target_exe.name
            new_icon_path = resource_dir / 'icon.icns'
            info_plist_file = contents_dir / 'Info.plist'
            plist_contents = f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple Computer//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
  <key>CFBundleGetInfoString</key>
  <string>{self.pretty_link_name}</string>
  <key>CFBundleExecutable</key>
  <string>{target_exe.name}</string>
  <key>CFBundleIdentifier</key>
  <string>com.plan-tools.www</string>
  <key>CFBundleName</key>
  <string>{target_exe.name}</string>
  <key>CFBundleIconFile</key>
  <string>icon</string>
  <key>CFBundleShortVersionString</key>
  <string>0.01</string>
  <key>CFBundleInfoDictionaryVersion</key>
  <string>6.0</string>
  <key>CFBundlePackageType</key>
  <string>APPL</string>
  <key>IFMajorVersion</key>
  <integer>0</integer>
  <key>IFMinorVersion</key>
  <integer>1</integer>
</dict>
</plist>"""
            self.desktop_file_data_check = {
                'contents_dir': contents_dir,
                'mac_os_dir': mac_os_dir,
                'resource_dir': resource_dir,
                'exe': new_exe_path,
                'icon': new_icon_path,
                'plist_path': info_plist_file,
                'plist_contents': plist_contents
            }
            if self.test_mode:
                return
            else:  # pragma: no cover
                # maybe someday we could actually test this in CI
                desktop_file.mkdir(parents=False)
                contents_dir.mkdir()
                mac_os_dir.mkdir()
                resource_dir.mkdir()
                with info_plist_file.open('w') as f:
                    f.write(plist_contents)
                symlink(target_exe, new_exe_path)
                if icon_file.exists():
                    copyfile(icon_file, new_icon_path)

    def get_pip_entry_point_exe_dir(self) -> Path:
        """
        Returns the path to the folder where the current executable will place entry_point binaries.
        On Linux, this will ultimately be a bin folder, while on Windows the folder name will be Scripts.
        For modern versions of Pip on Linux, a user-install will place binaries in USER_BASE/.local/bin,
        while global installs will install into /usr/bin or {venv}/bin.

        :return: Path instance where Pip is expected to have installed binaries
        """
        if system() == 'Windows':
            user_site_packages_dir = Path(getusersitepackages())
            user_scripts_dir = user_site_packages_dir.parent / 'Scripts'
            global_scripts_dir = Path(get_path('scripts'))
        elif system() == 'Linux':
            user_scripts_dir = Path(get_path('scripts'))
            global_scripts_dir = Path(USER_BASE) / 'bin'
        else:  # assuming 'Darwin'
            user_scripts_dir = Path(USER_BASE) / 'bin'
            global_scripts_dir = Path(get_path('scripts'))
        user_exe = user_scripts_dir / self.installed_binary_name
        global_exe = global_scripts_dir / self.installed_binary_name
        if user_exe.exists() and global_exe.exists():  # pragma: no cover
            # I'm going to allow this to be uncovered because it will be awkward to get this set up properly
            print(f"Detected the {self.installed_binary_name} binary in both user and global locations.")
            print("Due to this ambiguity, I cannot figure out to which one I should link.")
            print(f"User install location: {user_exe}")
            print(f"Global install location: {global_exe}")
            print("I am going to assume the user directory!!")
            return user_scripts_dir
        elif user_exe.exists():
            return user_scripts_dir
        elif global_exe.exists():
            return global_scripts_dir
        else:
            msg = f"Could not find {self.installed_binary_name} binary at either user or global location."
            msg += "This is weird since you are running this script...did you actually pip install this tool?"
            msg += "Make sure to pip install the tool and then retry"
            raise Exception(msg)

    def get_package_icon_for_link(self):
        if system() == 'Windows':
            return self.package_root / 'icons' / 'icon.ico'  # I think png might work here as well
        elif system() == 'Linux':
            return self.package_root / 'icons' / 'icon.png'
        else:  # assuming 'Darwin'
            return self.package_root / 'icons' / 'icon.icns'

    def get_path_to_link_file(self):
        if system() == 'Windows':
            from winreg import OpenKey, QueryValueEx, CloseKey, HKEY_CURRENT_USER as HKCU, KEY_READ as READ
            key = OpenKey(HKCU, r'Software\Microsoft\Windows\CurrentVersion\Explorer\User Shell Folders', 0, READ)
            desktop_value, _ = QueryValueEx(key, 'Desktop')
            CloseKey(key)
            desktop = Path(path.expandvars(desktop_value))
            link_name = f"{self.pretty_link_name}.lnk"
            return desktop / link_name
        elif system() == 'Linux':
            return Path.home() / '.local' / 'share' / 'applications' / f'{self.installed_binary_name}.desktop'
        else:  # assuming 'Darwin'
            return Path.home() / "Desktop" / f"{self.pretty_link_name}.app"  # actually a directory


def example_package_usage(_argv: List[str] = None) -> int:  # pragma: no cover
    """
    This is an example of the code that a packaged application should include to use this functionality.
    This function has an argument for options, so a wrapper function should be created to expose through an entry point.

    :param _argv: Currently, an unused argument, but eventually it could be how CLI args get passed in to
                  allow different run options
    :return: zero for success, nonzero otherwise
    """
    # Could handle arguments here
    source_dir = "cool_library"
    exe_name = "name_of_executable"
    nice_name = "Cool Program Name"
    s = EntryPoint(source_dir, exe_name, nice_name, "A really great tool description here", exe_name)
    return s.run()


def example_package_usage_wrapper() -> int:  # pragma: no cover
    """
    THis is an example of the wrapper function which is registered in the setup(entry_points...) section of setup.py.
    This function takes no Python arguments, but passes command line arguments to the worker function.

    :return: zero for success, nonzero otherwise
    """
    # Could handle argv here and pass them to example_package_usage(...)
    return example_package_usage(argv)


if __name__ == '__main__':  # pragma: no cover
    exit(example_package_usage_wrapper())
