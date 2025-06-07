#!/usr/bin/env python

"""
This file is a standalone EnergyPlus-Launch workflow tester.
It is called using a single argument, the path to a workflow file
"""

from inspect import isclass, getmembers
from pathlib import Path
from sys import argv, exit
from importlib import util as import_util


class WorkflowTesting:

    def __init__(self, verbose: bool):
        self.verbose = verbose

    def printer(self, message: str):
        if self.verbose:  # pragma: no cover
            print(message)

    def workflow_file_tester(self, file_path: str):
        modules = []

        path = Path(file_path)

        if path.exists():
            self.printer(f"   OK: File path exists at: {path}")
        else:  # pragma: no cover
            print(f"ERROR: File path does not exist!  Path: {path}")
            return 1

        if file_path.endswith('.py'):
            self.printer("   OK: File ends with .py")
        else:  # pragma: no cover
            print("ERROR: File path does NOT end with .py")
            return 1

        module_spec = import_util.spec_from_file_location('workflow_module', file_path)
        this_module = import_util.module_from_spec(module_spec)
        try:
            modules.append(this_module)
            module_spec.loader.exec_module(this_module)
            self.printer("   OK: Python import process completed successfully!")
        except ImportError as ie:  # pragma: no cover
            # this error generally means they have a bad workflow class or something
            print(f"ERROR: Import error occurred on workflow file {path}: {ie.msg}")
            return 1
        except SyntaxError as se:  # pragma: no cover
            # syntax errors are, well, syntax errors in the Python code itself
            print(f"ERROR: Syntax error occurred on workflow file {path}, line {se.lineno}: {se.msg}")
            return 1
        except Exception as e:  # pragma: no cover
            # there's always the potential of some other unforeseen thing going on when a workflow is executed
            print(f"ERROR: Unexpected error occurred trying to import workflow: {path}: {e}")
            return 1

        successful_classes = []
        for this_module in modules:
            class_members = getmembers(this_module, isclass)
            for this_class in class_members:
                this_class_name, this_class_type = this_class
                self.printer(f" INFO: Encountered class: \"{this_class_name}\", testing now...")
                # so right here, we could check issubclass, but this would also match the BaseEPLaunchWorkflow1, which
                # is imported in each workflow class.  No need to do that.  For now, I'm going to check the direct
                # parent class of this class to verify we only get direct descendants.  We can evaluate this later.
                # if issubclass(this_class_type, BaseEPLaunchWorkflow1):
                num_inheritance = len(this_class_type.__bases__)
                base_class_name = this_class_type.__bases__[0].__name__
                workflow_base_class_name = 'BaseEPLaunchWorkflow1'
                if num_inheritance == 1 and workflow_base_class_name in base_class_name:
                    self.printer(f"   OK: Basic inheritance checks out OK for class: {this_class_name}")
                    successful_classes.append(this_class_name)

                    try:
                        workflow_instance = this_class_type()
                        self.printer("   OK: Instantiation of derived class works")
                    except Exception as e:  # pragma: no cover
                        print(f"ERROR: Instantiation of derived class malfunctioning; reason: {e}")
                        return 1

                    required_method_names = [
                        "name", "description", "get_file_types",
                        "get_output_suffixes", "get_interface_columns", "context"
                    ]

                    for method in required_method_names:
                        try:
                            func = getattr(workflow_instance, method)
                            func()
                            self.printer(f"   OK: Overridden {method}() function execution works")
                        except Exception as e:  # pragma: no cover
                            print(f"ERROR: {method}() function not overridden, or malfunctioning; reason: {e}")
                            return 1

                else:
                    self.printer(" INFO: Inheritance does not check out, will continue with other classes in this file")
                    continue

        if len(successful_classes) > 0:
            t = '\n        - '
            class_list = t.join(successful_classes)
            self.printer(f"   OK: Found {len(successful_classes)} successful workflow imports:{t}{class_list}")
            return 0
        else:  # pragma: no cover
            print("ERROR: Did not find ANY successful workflow imports in this file!")
            return 1


def cli() -> int:  # pragma: no cover
    if len(argv) != 2:
        print("Bad call to workflow_tester.cli, give one command line argument, the full path to a workflow file")
        return 2
    else:
        return WorkflowTesting(verbose=True).workflow_file_tester(argv[1])


if __name__ == "__main__":  # pragma: no cover
    exit(cli())
