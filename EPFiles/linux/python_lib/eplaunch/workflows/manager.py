from pathlib import Path
from string import ascii_uppercase
from typing import Dict, List, Optional, Set, Tuple, Type
from importlib import util as import_util
from inspect import getmembers, isclass

from eplaunch.utilities.crossplatform import Platform
from eplaunch.workflows.base import BaseEPLaunchWorkflow1
from eplaunch.workflows.workflow import Workflow
from eplaunch.workflows.workflow_thread import WorkflowThread


class WorkflowManager:
    def __init__(self):
        self.current_workflow: Optional[Workflow] = None
        self.threads: Dict[str, WorkflowThread] = dict()
        self.workflow_directories: List[Path] = []
        self.auto_found_workflow_dirs: List[Path] = []
        self.workflows: List[Workflow] = []
        self.workflow_contexts: Set[str] = set()
        self.warnings: List[str] = []
        self.auto_find_workflow_directories()
        self.workflow_directories = self.auto_found_workflow_dirs

    def workflow_instances(self, workflow_context: str) -> List[Workflow]:
        return [x for x in self.workflows if x.context == workflow_context]

    def auto_find_workflow_directories(self) -> None:
        """Locate the (EnergyPlus) workflow directories that are in predestined locations"""
        self.auto_found_workflow_dirs = []
        # then search for e+ workflows
        if Platform.get_current_platform() == Platform.WINDOWS:
            from ctypes import windll
            bitmask = windll.kernel32.GetLogicalDrives()
            roots: List[Path] = []
            for letter in ascii_uppercase:
                if bitmask & 1:
                    roots.append(Path(f"{letter}:\\"))
                bitmask >>= 1
        elif Platform.get_current_platform() == Platform.LINUX:
            roots = [Path('/usr/local/bin/'), Path('/tmp/'), Path('/opt')]
        else:  # Assuming Platform.get_current_platform() == Platform.MAC:
            roots = [Path('/Applications/'), Path('/tmp/')]
        search_names = ["EnergyPlus*", "EP*", "ep*", "E+*", "e+*"]
        if Platform.get_current_platform() == Platform.LINUX:
            search_names.append("energyplus*")  # add lower case check on case-sensitive file systems (typically Linux)
        for search_root in roots:
            try:
                for search_name in search_names:
                    eplus_folder_matches = search_root.glob(search_name)
                    for ep_folder in eplus_folder_matches:  # pragma: no cover, would have to install in system folders
                        ep_workflow_dir = ep_folder / 'workflows'
                        if ep_workflow_dir.exists():
                            self.auto_found_workflow_dirs.append(ep_workflow_dir)
            except PermissionError:  # pragma: no cover -- this could be quite hard to reproduce :)
                continue  # just skip it, it could be like an empty DVD drive

    def instantiate_all_workflows(self, disable_builtins=False, extra_workflow_dir: Optional[Path] = None,
                                  skip_ep_search: bool = False) -> None:
        this_file_directory_path = Path(__file__).parent.resolve()
        this_project_root_dir = this_file_directory_path.parent
        built_in_workflow_dir = this_project_root_dir / 'workflows' / 'default'
        if skip_ep_search:
            all_workflow_directories = []
        else:
            all_workflow_directories = self.workflow_directories
        if disable_builtins:
            # don't add built-in default workflows
            pass
        elif built_in_workflow_dir not in all_workflow_directories and built_in_workflow_dir.exists():
            # add the built-in directory if it exists
            all_workflow_directories.append(built_in_workflow_dir)
        if extra_workflow_dir is not None:
            all_workflow_directories.append(extra_workflow_dir)

        self.workflows = []
        self.workflow_contexts.clear()
        self.warnings = []
        for i, workflow_directory in enumerate(all_workflow_directories):
            sanitized_directory_upper_case = str(workflow_directory).upper().replace('-', '.').replace('\\', '/')
            version_id = None
            dir_is_eplus = False
            # I tried regexes, and they worked using online Python regex testers, but using the same patterns
            # and strings in here resulting in false responses...bogus.  So here I go, manually chopping up a string
            # re_dots = re.compile('(?P<version>(\d.\d.\d))')
            if Platform.get_current_platform() == Platform.WINDOWS:  # pragma: no cover, skipping platform specifics
                energyplus_uc_search_string = 'ENERGYPLUSV'
            else:  # pragma: no cover, skipping platform specifics
                energyplus_uc_search_string = 'ENERGYPLUS.'
            if energyplus_uc_search_string in sanitized_directory_upper_case:
                dir_is_eplus = True
                san = sanitized_directory_upper_case
                trailing_string = san[san.index(energyplus_uc_search_string) + 11:]
                if '/' in trailing_string:
                    version_id = trailing_string[:trailing_string.index('/')]

            modules = []
            for this_file_path in workflow_directory.glob('*.py'):
                if this_file_path.name == '__init__.py':
                    continue
                module_spec = import_util.spec_from_file_location(('workflow_module_%s' % i), this_file_path)
                this_module = import_util.module_from_spec(module_spec)
                try:
                    modules.append([this_file_path, this_module])
                    module_spec.loader.exec_module(this_module)
                except ImportError as ie:
                    # this error generally means they have a bad workflow class or something
                    self.warnings.append(f"Import error occurred on workflow file {str(this_file_path)}: {ie.msg}")
                    continue
                except SyntaxError as se:
                    # syntax errors are, well, syntax errors in the Python code itself
                    self.warnings.append(f"Syntax error in workflow {str(this_file_path)}, line {se.lineno}: {se.msg}")
                    continue
                except Exception as e:  # pragma: no cover
                    # there's always the potential of some other unforeseen thing going on when a workflow is executed
                    self.warnings.append(f"Unexpected error importing workflow: {str(this_file_path)}: {str(e)}")
                    continue

            for module_file_path, this_module in modules:
                class_members: List[Tuple[str, Type[BaseEPLaunchWorkflow1]]] = getmembers(this_module, isclass)
                for this_class in class_members:
                    this_class_name, this_class_type = this_class
                    # so right here, we could check issubclass, but this also matches the BaseEPLaunchWorkflow1, which
                    # is imported in each workflow class.  No need to do that.  For now, I'm going to check the direct
                    # parent class of this class to verify we only get direct descendants.  We can evaluate this later.
                    # if issubclass(this_class_type, BaseEPLaunchWorkflow1):
                    num_inheritance = len(this_class_type.__bases__)
                    base_class_name = this_class_type.__bases__[0].__name__
                    workflow_base_class_name = 'BaseEPLaunchWorkflow1'
                    if num_inheritance == 1 and workflow_base_class_name in base_class_name:
                        try:
                            # we've got a good match, grab more data and get ready to load this into the Detail class
                            workflow_instance: BaseEPLaunchWorkflow1 = this_class_type()
                            workflow_name: str = workflow_instance.name()
                            workflow_file_types = workflow_instance.get_file_types()
                            workflow_output_suffixes = workflow_instance.get_output_suffixes()
                            workflow_columns = workflow_instance.get_interface_columns()
                            workflow_context = workflow_instance.context()
                            workflow_weather = workflow_instance.uses_weather()

                            file_type_string = "("
                            first = True
                            for file_type in workflow_file_types:
                                if first:
                                    first = False
                                else:
                                    file_type_string += ", "
                                file_type_string += file_type
                            file_type_string += ")"

                            description = f"{workflow_name} {file_type_string}"
                            self.workflow_contexts.add(workflow_context)

                            self.workflows.append(
                                Workflow(
                                    this_class_type,
                                    workflow_name,
                                    workflow_context,
                                    workflow_output_suffixes,
                                    workflow_file_types,
                                    workflow_columns,
                                    workflow_directory,
                                    description,
                                    dir_is_eplus,
                                    workflow_weather,
                                    version_id
                                )
                            )
                        except NotImplementedError as nme:
                            self.warnings.append(
                                f"Import error for file \"{module_file_path}\"; class: \"{this_class_name}\"; error: "
                                f"\"{str(nme)}\" "
                            )
                        except Exception as e:
                            # there's always the potential of some other thing going on when a workflow is executed
                            self.warnings.append(
                                f"Unexpected error in file \"{module_file_path}\"; class: \"{this_class_name}\"; error:"
                                f" \"{str(e)}\" "
                            )
                            continue

        self.workflows.sort(key=lambda w: w.description)

    def reset_workflow_array(self, filter_context=None) -> None:
        self.instantiate_all_workflows()
        if filter_context:
            self.workflows = [w for w in self.workflows if w.context == filter_context]
