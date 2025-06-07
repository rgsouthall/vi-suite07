from plan_tools.entry_point import EntryPoint
from eplaunch import NAME


def configure_cli() -> None:
    source_dir = "eplaunch"
    name = "energyplus_launch"
    description = "EnergyPlus Launch"
    s = EntryPoint(source_dir, name, "EnergyPlus Launch", description, NAME)
    s.run()
