from pathlib import Path


class TransitionBinary(object):
    """
    This class describes a single transition binary in the installation
    The class is constructed from the full path to the binary, then source and target versions are parsed from this path
    and the filenames

    :param full_path: Full path object to the transition binary itself

    :ivar full_path_to_binary: Copy of the full path to binary passed into the constructor
    :ivar binary_name: This is just the filename portion of the binary executable
    :ivar source_version: This is the source version of this particular transition, for example,
                          in V8-5-0-to-8-6-0, this will be 8.5
    :ivar target_version: This is the target version of this particular transition, for example,
                          in V8-5-0-to-8-6-0, this will be 8.6
    """

    def __init__(self, full_path: Path):
        self.full_path_to_binary = full_path
        self.binary_name = self.full_path_to_binary.name
        split_by_v = self.binary_name.split('V')
        source_token = split_by_v[1].split('-')
        source_string = f"{source_token[0]}.{source_token[1]}"
        target_token = split_by_v[2].split('-')
        target_string = f"{target_token[0]}.{target_token[1]}"
        self.source_version = float(source_string)
        self.target_version = float(target_string)

    def __str__(self) -> str:
        return f"TransitionBinary ({self.source_version} -> {self.target_version}) - {self.full_path_to_binary}"
