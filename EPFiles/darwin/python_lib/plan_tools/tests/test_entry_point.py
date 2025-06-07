from pathlib import Path
from platform import system
from unittest import TestCase

from plan_tools.entry_point import EntryPoint


class TestEntryPoint(TestCase):

    def test_entry_point_success(self):
        # it will be difficult to test the resulting paths on CI
        # we can at least test that it is returning good types
        # for now I'm going to try to test it with pip... :)
        e = EntryPoint("pip", "pip", "Setup Tools", "Descriptive", "wm-class", test_mode=True)
        e.run()
        found_exe_dir = e.get_pip_entry_point_exe_dir()
        self.assertIsInstance(found_exe_dir, Path)
        # TODO: Extend these to actually test the values for these keys
        if system() == 'Windows':
            self.assertIn('exe', e.desktop_file_data_check)
            self.assertIn('cwd', e.desktop_file_data_check)
            self.assertIn('icon', e.desktop_file_data_check)
        elif system() == 'Linux':
            self.assertIn('contents', e.desktop_file_data_check)
        else:  # assuming Darwin
            self.assertIn('contents_dir', e.desktop_file_data_check)
            self.assertIn('mac_os_dir', e.desktop_file_data_check)
            self.assertIn('resource_dir', e.desktop_file_data_check)
            self.assertIn('exe', e.desktop_file_data_check)
            self.assertIn('icon', e.desktop_file_data_check)
            self.assertIn('plist_path', e.desktop_file_data_check)
            self.assertIn('plist_contents', e.desktop_file_data_check)

    def test_entry_point_cannot_find_binary(self):
        e = EntryPoint("does not matter", "does not exist!", "Nice Name", "Descriptive", "wm-class", test_mode=True)
        with self.assertRaises(Exception):
            e.get_pip_entry_point_exe_dir()
