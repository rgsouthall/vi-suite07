import json
import os
from typing import Tuple

# TODO: this built-in eplaunch/utilities/version module should be removed, workflows should manage this themselves


class Version:

    @staticmethod
    def check_energyplus_version(file_path: str) -> Tuple[bool, str, int]:
        """Gets the version number information for a given EnergyPlus input file"""
        extension = os.path.splitext(file_path)[1].upper()
        if extension == '.IDF' or extension == '.IMF':
            return Version.check_idf_imf_energyplus_version(file_path)
        elif extension == '.EPJSON':
            return Version.check_json_energyplus_version(file_path)
        else:
            return False, '', 0

    @staticmethod
    def check_idf_imf_energyplus_version(file_path: str) -> Tuple[bool, str, int]:
        """Attempts to read a version number from an IDF syntax file"""
        # noinspection PyBroadException
        try:
            with open(file_path, "r") as f:
                cur_line = f.readline()
                while cur_line:
                    cur_line = cur_line.strip()
                    if len(cur_line) > 0 and cur_line[0] != "!" and "VERSION" in cur_line.upper():
                        cur_line = Version.line_with_no_comment(cur_line)
                        if ";" in cur_line:  # one liner version object ("Version, 8.4;")
                            poss_obj = cur_line
                        else:  # hoping for a two-liner version object ("Version,\n  8.4;")
                            next_line = Version.line_with_no_comment(f.readline().strip())
                            poss_obj = cur_line + next_line
                        hopeful_object = poss_obj[:-1] if poss_obj[-1] == ";" else poss_obj
                        fields = hopeful_object.split(',')
                        return True, fields[1], Version.numeric_version_from_string(fields[1])
                    cur_line = f.readline()  # get the next line
        except:  # noqa: E722
            return False, '', 0

    @staticmethod
    def line_with_no_comment(in_string: str) -> str:
        """Returns an IDF line with any comments removed"""
        return in_string[0:in_string.find("!")].strip() if in_string.find("!") >= 0 else in_string.strip()

    @staticmethod
    def numeric_version_from_string(string_version: str) -> int:
        """Gets the coded version string from a version string like 5.0.0-abcdef"""
        parts = [int(x) for x in string_version.split("-")[0].split('.')]
        return 10000 * parts[0] + 100 * parts[1]

    @staticmethod
    def numeric_version_from_dash_string(string_version: str) -> int:
        """Gets the coded version string from a version string like V5-0-0"""
        string_version = string_version[1:] if string_version[0] == 'V' else string_version
        # the rest of the version number should just be separated by periods
        parts = [int(x) for x in string_version.split("-")]
        return 10000 * parts[0] + 100 * parts[1]

    @staticmethod
    def string_version_from_number(version_number: int) -> str:
        """Converts a coded number like 50200 (fictional version 5.2) to string with leading zeros 'V050200'"""
        return 'V' + str(version_number).zfill(6)

    @staticmethod
    def check_json_energyplus_version(file_path: str) -> Tuple[bool, str, int]:
        """Reads a version number from an EpJSON file"""
        # noinspection PyBroadException
        try:
            with open(file_path, "r") as readfile:
                data = json.load(readfile)
            current_version = data['Version']['Version 1']['version_identifier']
            return True, current_version, Version.numeric_version_from_string(current_version)
        except:  # noqa: E722  # could be a file-not-found, a permission issue, a key error, a JSON exception.....
            return False, '', 0
