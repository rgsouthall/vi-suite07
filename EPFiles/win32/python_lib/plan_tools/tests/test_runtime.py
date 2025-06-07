from unittest import TestCase

from plan_tools.runtime import fixup_taskbar_icon_on_windows


class TestRunTime(TestCase):
    def test_it_runs(self):
        # to test this function, we'd need to make a call into the Windows DLLs, and it feels unnecessary
        # that's basically just recreating the implementation.  It should run happy though, so let it continue
        fixup_taskbar_icon_on_windows('hello')
