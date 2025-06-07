from pathlib import Path
from typing import Dict, List

from eplaunch import NAME, VERSION
from eplaunch.workflows.base import BaseEPLaunchWorkflow1, EPLaunchWorkflowResponse1


class ColumnNames:
    Location = 'Site:Location []'


class SiteLocationWorkflow(BaseEPLaunchWorkflow1):

    def name(self) -> str:
        return "Get Site:Location"

    def context(self) -> str:
        return f"{NAME} {VERSION}"

    def description(self) -> str:
        return "Retrieves the Site:Location name"

    def get_file_types(self) -> List[str]:
        return ["*.idf"]

    def get_output_suffixes(self) -> List[str]:
        return []

    def get_interface_columns(self) -> List[str]:
        return [ColumnNames.Location]

    def main(self, run_directory: Path, file_name: str, args: Dict) -> EPLaunchWorkflowResponse1:  # pragma: no cover
        self.callback(f"In {type(self).__name__}, about to process file: {file_name}")
        content = (run_directory / file_name).read_text()
        new_lines = []
        for line in content.split('\n'):
            if line.strip() == '':
                continue
            if '!' not in line:
                new_lines.append(line.strip())
            else:
                line_without_comment = line[0:line.index('!')].strip()
                if line_without_comment != '':
                    new_lines.append(line_without_comment)
        one_long_line = ''.join(new_lines)
        objects = one_long_line.split(';')
        for obj in objects:
            if obj.upper().startswith('SITE:LOCATION'):
                location_fields = obj.split(',')
                location_name = location_fields[1]
                break
        else:
            return EPLaunchWorkflowResponse1(
                success=False,
                message='Could not parse location object!',
                column_data={ColumnNames.Location: '*unknown*'}
            )
        self.callback(f"Completed {type(self).__name__}")
        return EPLaunchWorkflowResponse1(
            success=True,
            message='Parsed Location object successfully',
            column_data={ColumnNames.Location: location_name}
        )
