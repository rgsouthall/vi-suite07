"""pybind meshvis module"""
import netgen.libngpy.meshvis
import typing
import netgen.libngpy._meshing
__all__  = [
"VisualSceneMesh",
"GetGlobalMesh",
"MouseMove",
"SelectFace",
"VS",
"_Redraw"
]
class VisualSceneMesh():
    def Draw(self) -> None: ...
    pass
def GetGlobalMesh() -> netgen.libngpy._meshing.Mesh:
    pass
def MouseMove(arg0: VisualSceneMesh, arg1: int, arg2: int, arg3: int, arg4: int, arg5: str) -> None:
    pass
def SelectFace(arg0: int) -> None:
    pass
def VS(arg0: netgen.libngpy._meshing.Mesh) -> VisualSceneMesh:
    pass
def _Redraw(blocking: bool = False, fr: float = 25) -> bool:
    """
    Redraw all

    Parameters:

    blocking : bool
      input blocking

    fr : double
      input framerate
    """
