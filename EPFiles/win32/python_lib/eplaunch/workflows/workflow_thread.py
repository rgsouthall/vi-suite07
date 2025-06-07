from pathlib import Path
from threading import Thread
from typing import Dict, Callable

from eplaunch.workflows.base import EPLaunchWorkflowResponse1


class WorkflowThread(Thread):
    """Worker Thread Class."""

    def __init__(self, identifier: str, workflow_instance,
                 run_directory: Path, file_name: str, main_args: Dict, done_callback: Callable):
        super().__init__()
        self._want_abort = 0
        self.id = identifier
        self.workflow_instance = workflow_instance
        self.workflow_directory = main_args['workflow location']
        self.run_directory = run_directory
        self.file_name = file_name
        self.workflow_main_args = main_args
        self.workflow_done_callback = done_callback
        self.start()

    def run(self):
        """Run Workflow Thread."""
        try:
            workflow_response = self.workflow_instance.main(self.run_directory, self.file_name, self.workflow_main_args)
            if type(workflow_response) is not EPLaunchWorkflowResponse1:
                workflow_response = EPLaunchWorkflowResponse1(
                    success=False,
                    message='Current workflow main function did not respond properly',
                    column_data=None
                )
        except Exception as e:
            # here is another location where we simply don't know what types of errors could occur in user defined files
            workflow_response = EPLaunchWorkflowResponse1(
                success=False,
                message='Current workflow main function failed unexpectedly:' + str(e),
                column_data=None
            )
        workflow_response.id = self.id
        try:
            self.workflow_done_callback(workflow_response)
        except RuntimeError:  # pragma: no cover  -- this is an exceedingly odd case
            pass
            # print("Could not post finished event to the GUI, did the GUI get force closed?")

    def abort(self):
        """abort worker thread."""
        # Method for use by main thread to signal an abort
        self.workflow_instance.abort()
