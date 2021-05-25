"""pybind geom2d module"""
import netgen.libngpy._geom2d
import typing
import netgen.libngpy._meshing
__all__  = [
"Spline",
"SplineGeometry"
]
class Spline():
    """
    Spline of a SplineGeometry object
    """
    def EndPoint(self) -> netgen.libngpy._meshing.Point2d: ...
    def GetNormal(self, arg0: float) -> netgen.libngpy._meshing.Vec2d: ...
    def StartPoint(self) -> netgen.libngpy._meshing.Point2d: ...
    @property
    def bc(self) -> int:
        """
        :type: int
        """
    @property
    def leftdom(self) -> int:
        """
        :type: int
        """
    @leftdom.setter
    def leftdom(self, arg1: int) -> None:
        pass
    @property
    def rightdom(self) -> int:
        """
        :type: int
        """
    @rightdom.setter
    def rightdom(self, arg1: int) -> None:
        pass
    pass
class SplineGeometry(netgen.libngpy._meshing.NetgenGeometry):
    """
    a 2d boundary representation geometry model by lines and splines
    """
    def AddCurve(self, func: object, leftdomain: int = 1, rightdomain: int = 0, bc: object = <netgen.libngpy._meshing.NGDummyArgument object at 0x104aff330>, maxh: float = 1e+99) -> None: 
        """
        Curve is given as parametrization on the interval [0,1]
        """
    def Append(self, point_indices: list, leftdomain: int = 1, rightdomain: int = 0, bc: typing.Optional[Union[int, str]] = None, copy: typing.Optional[int] = None, maxh: float = 1e+99, hpref: float = 0, hprefleft: float = 0, hprefright: float = 0) -> int: ...
    def AppendPoint(self, x: float, y: float, maxh: float = 1e+99, hpref: float = 0, name: str = '') -> int: ...
    def AppendSegment(self, point_indices: list, leftdomain: int = 1, rightdomain: int = 0) -> None: ...
    def Draw(self) -> None: ...
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
        """
    def GetBCName(self, arg0: int) -> str: ...
    def GetNDomains(self) -> int: ...
    def GetNPoints(self) -> int: ...
    def GetNSplines(self) -> int: ...
    def GetPoint(self, arg0: int) -> netgen.libngpy._meshing.Point2d: ...
    def GetSpline(self, arg0: int) -> Spline: ...
    def Load(self, arg0: str) -> None: ...
    def PlotData(self) -> tuple: ...
    def PointData(self) -> tuple: ...
    def Print(self) -> None: ...
    def SegmentData(self) -> tuple: ...
    def SetDomainMaxH(self, arg0: int, arg1: float) -> None: ...
    def SetMaterial(self, arg0: int, arg1: str) -> None: ...
    def __getstate__(self) -> tuple: ...
    @typing.overload
    def __init__(self) -> None: ...
    @typing.overload
    def __init__(self, arg0: str) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def _visualizationData(self) -> dict: ...
    pass
