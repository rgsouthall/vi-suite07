from pathlib import Path
from tkinter import Tk, Toplevel, Frame, TOP, Button, scrolledtext, BOTH, Scrollbar, BOTTOM, X, END

from eplaunch.interface import set_frame_or_top_level_icon


class TkOutputDialog(Toplevel):

    def __init__(self, parent_window, workflow_id, workflow_title: str, configuration: str, x: int, y: int):
        super().__init__(parent_window)
        set_frame_or_top_level_icon(self, Path(__file__).resolve().parent.parent / 'icons')
        self.x = x
        self.y = y
        self.desired_width = 500
        self.title(workflow_title)
        self.workflow_id = workflow_id
        self._build_gui(configuration)

    def _build_gui(self, configuration: str):
        config_frame = Frame(self)
        horizontal_scroller = Scrollbar(config_frame, orient='horizontal')
        horizontal_scroller.pack(side=BOTTOM, fill=X)
        self.text_config = scrolledtext.ScrolledText(
            config_frame, wrap='none', width=30, height=6, xscrollcommand=horizontal_scroller.set
        )
        self.text_config.pack(side=TOP, padx=3, pady=3, fill=BOTH, expand=True)
        self.text_config.insert(END, configuration)
        horizontal_scroller.config(command=self.text_config.xview)
        config_frame.pack(side=TOP, expand=True, fill=BOTH, padx=5, pady=5)

        output_frame = Frame(self)
        horizontal_scroller = Scrollbar(output_frame, orient='horizontal')
        horizontal_scroller.pack(side=BOTTOM, fill=X)
        self.text_output = scrolledtext.ScrolledText(
            output_frame, wrap='none', width=30, height=16, xscrollcommand=horizontal_scroller.set
        )
        self.text_output.pack(side=TOP, padx=3, pady=3, fill=BOTH, expand=True)
        horizontal_scroller.config(command=self.text_output.xview)
        output_frame.pack(side=TOP, expand=True, fill=BOTH, padx=5, pady=5)

        Button(self, text="Close", command=self.close).pack(side=TOP, expand=True, padx=5, pady=5)
        self.after(20, self._size_window)

    def _size_window(self):
        current_height = max(self.winfo_height(), 440)
        self.geometry('%dx%d+%d+%d' % (500, current_height, self.x, self.y))

    def add_output(self, message: str):
        self.text_output.insert(END, message)
        self.text_output.insert(END, '\n')
        self.text_output.see(END)

    def close(self):
        self.destroy()


if __name__ == "__main__":
    root = Tk()
    root.title('Root Window for Toplevel Demo')
    file_listing = TkOutputDialog(root, 'id', 'title', 'config', 20, 100)
    root.mainloop()
