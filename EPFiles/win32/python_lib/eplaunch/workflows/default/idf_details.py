import os
from pathlib import Path
from typing import Dict, List

from eplaunch import NAME, VERSION
from eplaunch.workflows.base import BaseEPLaunchWorkflow1, EPLaunchWorkflowResponse1


class ColumnNames:
    Version = 'Version'
    NumDesignDays = '# Design Days'
    NumRunPeriods = '# Run Periods'
    NumZones = '# Zones'


class IDFDetailsWorkflow1(BaseEPLaunchWorkflow1):

    def name(self) -> str:
        return "IDF Details"

    def context(self) -> str:
        return f"{NAME} {VERSION}"

    def description(self) -> str:
        return "Retrieves IDF Details"

    def get_file_types(self) -> List[str]:
        return ["*.idf"]

    def get_output_suffixes(self) -> List[str]:
        return []

    def get_interface_columns(self) -> List[str]:
        return [
            ColumnNames.Version,
            ColumnNames.NumDesignDays,
            ColumnNames.NumRunPeriods,
            ColumnNames.NumZones
        ]

    def main(self, run_directory: Path, file_name: str, args: Dict) -> EPLaunchWorkflowResponse1:  # pragma: no cover
        self.callback(f"In {type(self).__name__}, about to process file: {file_name}")
        file_path = os.path.join(run_directory, file_name)
        content = open(file_path).read()
        new_lines = []
        try:
            for line in content.split('\n'):
                if line.strip() == '':
                    continue
                if '!' not in line:
                    new_lines.append(line.strip())
                else:
                    line_without_comment = line[0:line.index('!')].strip()
                    if line_without_comment != '':
                        new_lines.append(line_without_comment)
        except Exception as e:
            self.callback("Could not process IDF; error: " + str(e))
            return EPLaunchWorkflowResponse1(
                success=False,
                message='Could not parse IDF object data!',
                column_data={
                    ColumnNames.Version: '*unknown*',
                    ColumnNames.NumDesignDays: '*unknown*',
                    ColumnNames.NumRunPeriods: '*unknown*',
                    ColumnNames.NumZones: '*unknown*'
                }
            )
        one_long_line = ''.join(new_lines)
        objects = one_long_line.split(';')
        num_dds = 0
        num_rps = 0
        num_zones = 0
        version_id = '*unknown*'
        for obj in objects:
            if obj.upper().startswith('SIZINGPERIOD:DESIGNDAY'):
                num_dds += 1
            elif obj.upper().startswith('RUNPERIOD'):
                num_rps += 1
            elif obj.upper().startswith('ZONE,'):
                num_zones += 1
            elif obj.upper().startswith('VERSION,'):
                version_fields = obj.split(',')
                version_id = version_fields[1]
        self.callback(f"Completed {type(self).__name__}")
        return EPLaunchWorkflowResponse1(
            success=True,
            message='Parsed IDF information successfully',
            column_data={
                ColumnNames.Version: version_id,
                ColumnNames.NumDesignDays: num_dds,
                ColumnNames.NumRunPeriods: num_rps,
                ColumnNames.NumZones: num_zones
            }
        )
