from eplaunch.interface.frame_main import EPLaunchWindow


def main_gui(called_from_ep_cli: bool = False):
    EPLaunchWindow(called_from_ep_cli).run()


if __name__ == "__main__":
    main_gui()
