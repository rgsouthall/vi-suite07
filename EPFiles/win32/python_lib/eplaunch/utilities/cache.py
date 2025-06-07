from json import dumps, loads
from os import path
from pathlib import Path
from time import sleep
from typing import Dict, List

from eplaunch.utilities.exceptions import EPLaunchFileException

# I'm trying to be good and catch specific exceptions, but the json library makes it difficult like this :-D
try:
    from json.decoder import JSONDecodeError
except ImportError:  # pragma: no cover
    JSONDecodeError = ValueError

#: This is used as the mutex queue, the list of unique directories being altered at a given time
cache_files_currently_updating_or_writing: List[Path] = []


class CacheFile:
    """
    Represents the file that is kept in each folder where workflows have been started
    Keeps track of the most recent state of the file, with some metadata that is workflow dependent

    Usage:

    To ensure thread-safety, this class employs a form of a mutex, where the unique id is the current directory
    Any worker function that wants to alter the queue should follow the following process:

    - The worker should call the ok_to_continue() function, which will check the mutex and then wait a predetermined
      amount of time for the mutex to clear, or fail.
    - The worker should check the return value of this function and if False, fail.  If True, it should set up a block
      on the directory by adding the current directory to the cache_files_currently_updating_or_writing array
    - The worker can then proceed to read the cache, modify ir, and write to disk
    - The worker must then release the mutex by removing the current directory from the list
    """

    FileName = '.eplaunch'
    RootKey = 'workflows'
    FilesKey = 'files'
    ParametersKey = 'config'
    ResultsKey = 'result'
    WeatherFileKey = 'weather'
    QueueCheckInterval = 0.1  # seconds
    QueueTotalCheckTime = 5  # seconds

    def _print(self, message) -> None:
        """
        Utility function for printing diagnostic messages -- useful for when debugging synchronous alterations

        :param message: The message to print
        :return: None
        """
        debug = False
        if debug:  # pragma: no cover
            print(f"{self.file_path}: {message}")

    def __init__(self, working_directory: Path):
        """
        Constructor for this class, stores the local file path and initializes the workflow_state

        :param working_directory:
        """
        self.file_path = working_directory / self.FileName
        self._print("Created cache file")
        self.workflow_state = None

    def _add_file_attribute(self, workflow_name, file_name, attribute, data, replace) -> None:
        """
        This function generically updates some attribute of a file within a given workflow context
        The hierarchy is:
         workflows
          - workflow_name
           - files
            - file_name
             - attribute
              - data
        The replace parameter states whether the content in attribute will be updated, or replaced, with data

        :param workflow_name: The name of the workflow to alter, as given by the workflow's name() method
        :param file_name: The file name of the file to alter
        :param attribute: The attribute to alter, currently the only two options are 'config' or 'result'
        :param data: A map of data to write to this attribute
        :param replace: A flag for whether this data should replace all prior data or just append to it
        :return: None
        """

        # if there is already a config for this workflow/file, update it
        # if something is missing from the structure, initialize it on each stage
        self._print(f"Adding file attribute for workflow: {workflow_name}, file: {file_name}")
        root = self.workflow_state[self.RootKey]
        if workflow_name in root:
            this_workflow = root[workflow_name]
            if self.FilesKey in this_workflow:
                these_files = this_workflow[self.FilesKey]
                if file_name in these_files:
                    this_file = these_files[file_name]
                    if replace:
                        this_file[attribute] = data
                    else:
                        if attribute in this_file:
                            this_config = this_file[attribute]
                            merged_config_data = {**this_config, **data}  # merge dicts, avail in Python 3.5+
                            this_file[attribute] = merged_config_data
                        else:
                            this_file[attribute] = data
                else:
                    these_files[file_name] = {attribute: data}
            else:  # pragma: no cover
                # There's really no way to get here...but I feel like I should leave this in
                this_workflow[self.FilesKey] = {file_name: {attribute: data}}
        else:
            root[workflow_name] = {self.FilesKey: {file_name: {attribute: data}}}

    def ok_to_continue(self) -> bool:
        """
        This function does the check-and-wait part of the mutex.  If the current directory is not blocked, it
        immediately returns.  If the current directory is blocked, it will attempt to check over a certain amount of
        time, at a tight interval, to wait on the mutex to be unlocked.  Ultimately if it can't pass, it returns False.

        :return: True or False, whether it is safe to write to this cache
        """
        self._print("Checking if its ok to continue")
        if self.file_path not in cache_files_currently_updating_or_writing:
            return True
        self._print("Found this file in the writing data, trying to sleep through it")
        for i in range(int(self.QueueTotalCheckTime / self.QueueCheckInterval)):
            sleep(self.QueueCheckInterval)
            if self.file_path not in cache_files_currently_updating_or_writing:
                self._print("Managed to sleep long enough, continuing!")
                return True
        else:
            self._print("Sleep didn't last long enough, aborting")
            return False
        # there is an **incredibly** small chance we could have a new file pop in between the check above and later code
        # I will have to noodle on whether we want to worry about that possibility

    def add_config(self, workflow_name, file_name, config_data) -> None:
        """
        This function is used to add a config data block for a workflow.  A config data block contains data that is
        generally thought of as "input data" for a workflow, such as a weather file for a simulation run.

        :param workflow_name: The name of the workflow to alter, as given by the workflow's name() method
        :param file_name: The file name of the file to alter
        :param config_data: A map of data to write to this config section
        :return: None
        """
        self._print(f"About to add a config attribute for workflow {workflow_name}; file {file_name}")
        if not self.ok_to_continue():  # pragma: no cover
            # I'm not sure if we should communicate this or not.
            # For now let's print to the term, so that we might catch it for debugging
            print("There was a problem adding a config to a cache file, something blocked for too long maybe?")
        cache_files_currently_updating_or_writing.append(self.file_path)
        self._print("Cache file locked")
        self.read()
        self._add_file_attribute(workflow_name, file_name, self.ParametersKey, config_data, False)
        self.write()
        cache_files_currently_updating_or_writing.remove(self.file_path)
        self._print("Cache file UN-locked")

    def add_result(self, workflow_name, file_name, column_data) -> None:
        """
        This function is used to add a result data block for a workflow.  A result data block contains data that is
        generally thought of as "output data" for a workflow, such as energy usage for a simulation run.

        :param workflow_name: The name of the workflow to alter, as given by the workflow's name() method
        :param file_name: The file name of the file to alter
        :param column_data: A map of data to write to this result section, the keys are expected to be defined by
                            the workflow itself as given by the get_interface_columns() method
        :return: None
        """
        self._print(f"About to add a result attribute for workflow {workflow_name}; file {file_name}")
        if not self.ok_to_continue():  # pragma: no cover
            # I'm not sure if we should communicate this or not.
            # For now let's print to the term, so that we might catch it for debugging
            print("There was a problem adding a config to a cache file, something blocked for too long maybe?")
        cache_files_currently_updating_or_writing.append(self.file_path)
        self._print("Cache file locked")
        self.read()
        self._add_file_attribute(workflow_name, file_name, self.ResultsKey, column_data, True)
        self.write()
        cache_files_currently_updating_or_writing.remove(self.file_path)
        self._print("Cache file UN-locked")

    def read(self) -> None:
        """
        Reads the existing cache file, if it exists, and stores the data in the workflow_state instance variable.
        If the cache file doesn't exist, this simply initializes the workflow_state instance variable.

        :return: None
        """
        if path.exists(self.file_path):
            try:
                with open(self.file_path, 'r') as f:
                    body_text = f.read()
            except IOError:  # pragma: no cover  -- would be difficult to mock up this weird case
                raise EPLaunchFileException(self.file_path, 'Could not open or read text from file')
            try:
                self.workflow_state = loads(body_text)
            except JSONDecodeError:
                raise EPLaunchFileException(self.file_path, 'Could not parse cache file JSON text')
        else:
            self.workflow_state = {self.RootKey: {}}

    def write(self) -> None:
        """
        Writes out the workflow state to the previously determined cache file location
        Note that this function does not protect for thread-safety!  It is expected that functions who are
        altering the state of the cache should call write() within their own blocking structure

        :return: None
        """
        body_text = dumps(self.workflow_state, indent=2)
        try:
            with open(self.file_path, 'w') as f:
                f.write(body_text)
        except IOError:  # pragma: no cover  -- would be difficult to mock up this weird case
            raise EPLaunchFileException(self.file_path, 'Could not write cache file')

    def get_files_for_workflow(self, current_workflow_name) -> Dict:
        """
        Gets a list of files that are found in this cache inside the given workflow name

        :param current_workflow_name: The name of a workflow (as determined by the name() function on the workflow)
        :return: A map with keys that are file names found in this workflow
        """
        self.read()
        workflows = self.workflow_state[CacheFile.RootKey]
        if current_workflow_name in workflows:
            return workflows[current_workflow_name][CacheFile.FilesKey]
        else:
            return {}
