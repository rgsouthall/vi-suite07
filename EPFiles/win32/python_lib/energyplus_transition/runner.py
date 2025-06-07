def main_cli():
    pass


def main_gui(called_from_ep_cli: bool = False):
    from energyplus_transition.gui import VersionUpdaterWindow

    # we will keep the form in a loop to handle requested restarts (language change, etc.)
    running = True
    while running:
        main_window = VersionUpdaterWindow(called_from_ep_cli)
        main_window.mainloop()
        running = main_window.doing_restart


if __name__ == "__main__":
    main_gui()
