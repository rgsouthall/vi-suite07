"""pybind csgvis module"""
import netgen.libngpy.csgvis
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import netgen.libngpy._csg
__all__  = [
"VisualSceneGeometry",
"MouseMove",
"SetBackGroundColor",
"VS"
]
class VisualSceneGeometry():
    def Draw(self) -> None: ...
    pass
def MouseMove(arg0: VisualSceneGeometry, arg1: int, arg2: int, arg3: int, arg4: int, arg5: str) -> None:
    pass
def SetBackGroundColor(arg0: float) -> None:
    pass
def VS(arg0: netgen.libngpy._csg.CSGeometry) -> VisualSceneGeometry:
    pass
