from platform import system


class Platform:
    WINDOWS = 1
    LINUX = 2
    MAC = 3
    UNKNOWN = 4

    @staticmethod
    def get_current_platform(test_name: str = None) -> int:
        if test_name:
            platform_name = test_name
        else:  # pragma: no cover -- can't know ahead of time which system we will test on
            platform_name = system()
        if platform_name == 'Windows':
            return Platform.WINDOWS
        elif platform_name == 'Linux':
            return Platform.LINUX
        elif platform_name == 'Darwin':
            return Platform.MAC
        else:
            return Platform.UNKNOWN
