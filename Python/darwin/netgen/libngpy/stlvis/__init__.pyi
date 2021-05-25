"""pybind stlvis module"""
import netgen.libngpy.stlvis
import typing
import netgen.libngpy._stl
__all__  = [
"VisualSceneSTLGeometry",
"SetBackGroundColor",
"VS"
]
class VisualSceneSTLGeometry():
    def Draw(self) -> None: ...
    pass
def SetBackGroundColor(arg0: float) -> None:
    pass
def VS(arg0: netgen.libngpy._stl.STLGeometry) -> VisualSceneSTLGeometry:
    pass
