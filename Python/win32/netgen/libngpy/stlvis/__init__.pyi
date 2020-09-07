"""pybind stlvis module"""
import netgen.libngpy.stlvis
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
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
