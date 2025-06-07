import os
from pathlib import Path
from typing import Dict, List

from eplaunch import NAME, VERSION
from eplaunch.workflows.base import BaseEPLaunchWorkflow1, EPLaunchWorkflowResponse1


class ColumnNames:
    FileType = 'File Type'
    FileSize = 'File Size [kB]'


class FileDetailsWorkflow1(BaseEPLaunchWorkflow1):

    def name(self) -> str:
        return "Get File Details"

    def context(self) -> str:
        return f"{NAME} {VERSION}"

    def description(self) -> str:
        return "Retrieves details about the file"

    def get_file_types(self) -> List[str]:
        return ["*"]

    def get_output_suffixes(self) -> List[str]:
        return []

    def get_interface_columns(self) -> List[str]:
        return [ColumnNames.FileType, ColumnNames.FileSize]

    def main(self, run_directory: Path, file_name: str, args: Dict) -> EPLaunchWorkflowResponse1:  # pragma: no cover
        self.callback(f"In {type(self).__name__}, about to process file: {file_name}")
        file_path = os.path.join(run_directory, file_name)
        info = os.stat(file_path)
        file_type = os.path.splitext(file_path)[1]
        file_size = round(info.st_size / 1024)
        self.callback(f"Completed {type(self).__name__}")
        return EPLaunchWorkflowResponse1(
            success=True,
            message='Parsed File data successfully',
            column_data={
                ColumnNames.FileSize: file_size,
                ColumnNames.FileType: file_type,
            }
        )
