from sys import platform


def fixup_taskbar_icon_on_windows(entry_point_id: str) -> None:
    """
    OK, so with setuptools, we can declare a GUI entry point as such:

    .. code-block:: python

      setup(
        ...
        entry_points={
          'gui_scripts': [
            'resulting_exe_name=path.to:function',
          ],
        },
        ...
      )

    When a user Pip installs the library, it will do a follow-up step to create
    an executable entry point to run the GUI, and place it in an appropriate location.
    Internally, this executable is a Python wrapper around the entry point function call.
    In the GUI source code, it's easy enough to set the icon shown on the GUI.
    For TkInter, this is usually with a call to root.iconbitmap() or root.iconphoto().
    When the program launches, this icon will display properly in the menu bar.
    But for Pip entry points, the taskbar icon will be the default Python icon.
    Windows detects that Python is executing the program, so that's the icon it uses.
    To date, the only identified override is to use a Windows function call to set the app ID.
    This is often used to manipulate the taskbar icon grouping for more complex GUI apps, but it works here.
    By overriding the ModelID to something unique for this program, the custom icon appears.

    :param entry_point_id: Unique for this program...just use the program name.
    :return: Nothing, if it fails, it just quietly fails.
    """
    # noinspection PyBroadException
    try:
        if platform.startswith('win'):
            import ctypes
            ctypes.windll.shell32.SetCurrentProcessExplicitAppUserModelID(f"Plan.Tools.Program.{entry_point_id}")
    except Exception:  # pragma: no cover
        pass
