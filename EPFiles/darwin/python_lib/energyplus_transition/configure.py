from plan_tools.entry_point import EntryPoint  # type: ignore
from energyplus_transition import NAME


def configure_cli() -> None:
    source_dir = "energyplus_transition"
    name = "energyplus_transition_gui"
    description = "EnergyPlus Transition, for Updating Input Files"
    s = EntryPoint(source_dir, name, "EnergyPlus Transition", description, NAME)
    s.run()
