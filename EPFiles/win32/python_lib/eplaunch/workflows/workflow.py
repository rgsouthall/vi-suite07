from pathlib import Path
from typing import List, Type

from eplaunch.workflows.base import BaseEPLaunchWorkflow1


class Workflow:
    def __init__(self,
                 workflow_class: Type[BaseEPLaunchWorkflow1],
                 name: str,
                 context: str,
                 output_suffixes: List[str],
                 file_types: List[str],
                 columns: List[str],
                 directory: Path,
                 description: str,
                 is_energyplus: bool,
                 uses_weather: bool,
                 version_id: str):
        self.workflow_class = workflow_class
        self.name = name
        self.context = context
        self.output_suffixes = output_suffixes
        self.file_types = file_types
        self.columns = columns
        self.workflow_directory = directory
        self.description = description
        self.is_energyplus = is_energyplus
        self.uses_weather = uses_weather
        self.version_id = version_id

    def __str__(self) -> str:
        return f"Workflow {self.context}:{self.name}"
