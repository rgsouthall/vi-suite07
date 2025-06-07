from pathlib import Path
from subprocess import Popen, PIPE, STDOUT, CalledProcessError
from typing import Callable, Dict, List, Optional


class EPLaunchWorkflowResponse1:

    def __init__(self, success: bool, message: str, column_data: Optional[Dict], **extra_data: Optional[Dict]):
        self.id: Optional[str] = None  # assigned by workflow thread manager
        self.success = success
        self.message = message
        self.column_data = column_data
        self.extra_data = extra_data


class BaseEPLaunchWorkflow1:

    def __init__(self):
        self.my_id: Optional[str] = None  # will be a UUID string
        self._callback: Optional[Callable] = None  # callback instance to pass messages back up to the GUI
        self._process: Optional[Popen] = None  # process instance if the workflow wants to use execute_for_callback

    def name(self) -> str:
        raise NotImplementedError("name function needs to be implemented in derived workflow class")

    def context(self) -> str:
        raise NotImplementedError("context function needs to be implemented in derived workflow class")

    def description(self) -> str:
        raise NotImplementedError("description function needs to be implemented in derived workflow class")

    # noinspection PyMethodMayBeStatic
    def uses_weather(self) -> bool:
        """
        If it returns True, this workflow accepts a "weather" key in the arguments to the workflow
        :return: Boolean
        """
        return False

    def get_file_types(self) -> List[str]:
        raise NotImplementedError("get_file_types needs to be implemented in derived workflow class")

    def get_output_suffixes(self) -> List[str]:
        raise NotImplementedError("get_output_suffixes needs to be implemented in derived workflow class")

    # noinspection PyMethodMayBeStatic
    def get_extra_data(self) -> Dict:
        """
        Allows a dictionary of extra data to be generated, defaults to empty, so it is not required
        :return: Dictionary of string, string
        """
        return {}

    def get_interface_columns(self) -> List[str]:
        """
        Returns an array of column names for the interface; defaults to empty, so it is not required
        :return: A list of interface column names
        """
        return []

    def register_standard_output_callback(self, workflow_id: str, callback: Callable) -> None:
        """
        Used to register the callback function from the UI for standard output from this workflow.
        This function is not to be inherited by derived workflows unless they are doing something really odd.
        Workflows should simply use self.callback(message) to send messages as necessary to the GUI during a workflow.

        :param workflow_id: A unique ID assigned by the program in order to track workflows
        :param callback: The GUI function to be called with message updates
        :return: None
        """
        self.my_id = workflow_id
        self._callback = callback

    def callback(self, message: str):
        """
        This is the actual callback function that workflows can call when they have an update.
        Internally here there are some other parameters passed up, but that is just to further isolate the users
        from having to pass extra data during their own calls

        :param message: A message to be sent to the GUI from the workflow
        :return: None
        """
        self._callback(self.my_id, message)

    def main(self, run_directory: Path, file_name: str, args: Dict) -> EPLaunchWorkflowResponse1:
        """
        The actual running operation for the workflow, should check `self.abort` periodically to allow exiting
        :return: Must return an EPLaunchWorkflowResponse1 instance
        """
        raise NotImplementedError("main function needs to be implemented in derived workflow class")

    def abort(self) -> None:
        if self._process:  # pragma: no cover  # not getting caught by coverage tools, not sure why
            self._process.kill()

    def execute_for_callback(self, command_line_tokens: List[str], working_directory: str):  # pragma: no cover
        # In Python 3.10+, Popen objects can be context managers, like such:
        # with Popen(["ifconfig"], stdout=PIPE) as proc:
        #     log.write(proc.stdout.read())
        # could be a useful change later
        self._process = Popen(
            command_line_tokens, cwd=working_directory, stdout=PIPE, stderr=STDOUT, universal_newlines=True
        )
        for stdout_line in iter(self._process.stdout.readline, ""):
            yield stdout_line.strip()
        self._process.stdout.close()
        return_code = self._process.wait()
        if return_code:
            raise CalledProcessError(return_code, command_line_tokens)
