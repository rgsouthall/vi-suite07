"""pybind stl module"""
import netgen.libngpy._stl
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import netgen.libngpy._meshing
__all__  = [
"STLGeometry",
"LoadSTLGeometry"
]
class STLGeometry(netgen.libngpy._meshing.NetgenGeometry):
    def GenerateMesh(self, mp: netgen.libngpy._meshing.MeshingParameters = None, **kwargs) -> netgen.libngpy._meshing.Mesh: 
        """
        Meshing Parameters
        -------------------

        maxh: float = 1e10
          Global upper bound for mesh size.

        grading: float = 0.3
          Mesh grading how fast the local mesh size can change.

        meshsizefilename: str = None
          Load meshsize from file. Can set local mesh size for points
          and along edges. File must have the format:

            nr_points
            x1, y1, z1, meshsize
            x2, y2, z2, meshsize
            ...
            xn, yn, zn, meshsize

            nr_edges
            x11, y11, z11, x12, y12, z12, meshsize
            ...
            xn1, yn1, zn1, xn2, yn2, zn2, meshsize

        segmentsperedge: float = 1.
          Minimal number of segments per edge.

        quad_dominated: bool = False
          Quad-dominated surface meshing.

        blockfill: bool = True
          Do fast blockfilling.

        filldist: float = 0.1
          Block fill up to distance

        delaunay: bool = True
          Use delaunay meshing.

        Optimization Parameters
        -----------------------

        optimize3d: str = "cmdmustm"
          3d optimization strategy:
            m .. move nodes
            M .. move nodes, cheap functional
            s .. swap faces
            c .. combine elements
            d .. divide elements
            p .. plot, no pause
            P .. plot, Pause
            h .. Histogramm, no pause
            H .. Histogramm, pause

        optsteps3d: int = 3
          Number of 3d optimization steps.

        optimize2d: str = "smcmSmcmSmcm"
          2d optimization strategy:
            s .. swap, opt 6 lines/node
            S .. swap, optimal elements
            m .. move nodes
            p .. plot, no pause
            P .. plot, pause
            c .. combine

        optsteps2d: int = 3
          Number of 2d optimization steps.

        elsizeweight: float = 0.2
          Weight of element size w.r.t. element shape in optimization.


        STL Specific Meshing Parameters
        -------------------------------

        yangle: float = 30.
          Angle for edge detection

        contyangle: float = 20.
          Edges continue if angle > contyangle

        edgecornerangle: float = 60.
          Angle of geometry edge at which the mesher should set a point.

        closeedgefac: Optional[float] = 1.
          Factor for meshing close edges, if None it is disabled.

        minedgelen: Optional[float] = 0.001
          Minimum edge length to be used for dividing edges to mesh points. If
          None this is disabled.
        """
    def __getstate__(self) -> tuple: ...
    @overload
    def __init__(self, filename: str) -> None: ...
    @overload
    def __init__(self) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def _visualizationData(self) -> dict: ...
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x00000166FFF90B28>, '__dict__': <attribute '__dict__' of 'netgen.libngpy._stl.STLGeometry' objects>, '__doc__': None, '__module__': 'netgen.libngpy._stl', '__getstate__': <instancemethod __getstate__ at 0x00000166FFF90B58>, '__setstate__': <instancemethod __setstate__ at 0x00000166FFF90BB8>, '_visualizationData': <instancemethod _visualizationData at 0x00000166FFF90C18>, 'GenerateMesh': <instancemethod GenerateMesh at 0x00000166FFF90C78>})
    pass
def LoadSTLGeometry(arg0: str) -> STLGeometry:
    pass
