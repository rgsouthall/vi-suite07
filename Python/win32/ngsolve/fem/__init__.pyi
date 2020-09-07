"""pybind fem"""
import ngsolve.fem
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import ngsolve.ngstd
import numpy
import ET
import ngsolve.bla
import pyngcore
__all__  = [
"BFI",
"BSpline",
"BaseMappedIntegrationPoint",
"CoefficientFunction",
"DifferentialOperator",
"ET",
"ElementTopology",
"ElementTransformation",
"FiniteElement",
"HCurlFE",
"HDivDivFE",
"HDivFE",
"IntegrationPoint",
"IntegrationRule",
"LFI",
"MeshPoint",
"MixedFE",
"NODE_TYPE",
"Parameter",
"ParameterC",
"ScalarFE",
"SpecialCFCreator",
" ",
"BlockBFI",
"BlockLFI",
"CacheCF",
"Cof",
"CompoundBFI",
"CompoundLFI",
"Conj",
"CoordCF",
"Cross",
"Det",
"GenerateL2ElementCode",
"H1FE",
"Id",
"IfPos",
"Inv",
"L2FE",
"LoggingCF",
"SetPMLParameters",
"Skew",
"Sym",
"Trace",
"VoxelCoefficient",
"acos",
"asin",
"atan",
"atan2",
"ceil",
"cos",
"cosh",
"exp",
"floor",
"log",
"pow",
"sin",
"sinh",
"sqrt",
"tan",
"CELL",
"EDGE",
"ELEMENT",
"FACE",
"FACET",
"HEX",
"POINT",
"PRISM",
"PYRAMID",
"QUAD",
"SEGM",
"TET",
"TRIG",
"VERTEX",
"specialcf"
]
class BFI():
    """
    Bilinear Form Integrator

    Parameters:

    name : string
      Name of the bilinear form integrator.

    py_coef : object
      CoefficientFunction of the bilinear form.

    dim : int
      dimension of the bilinear form integrator

    imag : bool
      Multiplies BFI with 1J

    filename : string
      filename 

    kwargs : kwargs
      For a description of the possible kwargs have a look a bit further down.
    """
    def CalcElementMatrix(self, fel: FiniteElement, trafo: ElementTransformation, heapsize: int = 10000, complex: bool = False) -> object: 
        """
         
        Calculate element matrix of a specific element.

        Parameters:

        fel : ngsolve.fem.FiniteElement
          input finite element

        trafo : ngsolve.fem.ElementTransformation
          input element transformation

        heapsize : int
          input heapsize

        complex : bool
          input complex
        """
    def Evaluator(self, name: str) -> DifferentialOperator: 
        """
        Returns requested evaluator

        Parameters:

        name : string
          input name of requested evaluator
        """
    def GetDefinedOn(self) -> pyngcore.BitArray: 
        """
        Returns a BitArray where the bilinear form is defined on
        """
    def SetDefinedOnElements(self, bitarray: pyngcore.BitArray) -> None: 
        """
         
        Set the elements on which the bilinear form is defined on.

        Parameters:

        bitarray : ngsolve.ngstd.BitArray
          input bitarray
        """
    def SetIntegrationRule(self, et: ET, intrule: IntegrationRule) -> BFI: 
        """
         
        Set integration rule of the bilinear form.

        Parameters:

        et : ngsolve.fem.Element_Type
          input element type

        intrule : ngsolve.fem.Integrationrule
          input integration rule
        """
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __init__(self, name: str = '', coef: object, dim: int = -1, imag: bool = False, filename: str = '', **kwargs) -> None: ...
    def __initialize__(self, **kwargs) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    @property
    def simd_evaluate(self) -> bool:
        """
        SIMD evaluate ?

        :type: bool
        """
    @simd_evaluate.setter
    def simd_evaluate(self, arg1: bool) -> None:
        """
        SIMD evaluate ?
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE184798>, '__dict__': <attribute '__dict__' of 'ngsolve.fem.BFI' objects>, '__doc__': '\nBilinear Form Integrator\n\nParameters:\n\nname : string\n  Name of the bilinear form integrator.\n\npy_coef : object\n  CoefficientFunction of the bilinear form.\n\ndim : int\n  dimension of the bilinear form integrator\n\nimag : bool\n  Multiplies BFI with 1J\n\nfilename : string\n  filename \n\nkwargs : kwargs\n  For a description of the possible kwargs have a look a bit further down.\n\n', '__module__': 'ngsolve.fem', '__flags_doc__': <staticmethod object at 0x0000026DEE186148>, '__special_treated_flags__': <staticmethod object at 0x0000026DEE181888>, '__initialize__': <instancemethod __initialize__ at 0x0000026DEE184858>, 'simd_evaluate': <property object at 0x0000026DEE188E58>, '__str__': <instancemethod __str__ at 0x0000026DEE184918>, 'Evaluator': <instancemethod Evaluator at 0x0000026DEE184978>, 'GetDefinedOn': <instancemethod GetDefinedOn at 0x0000026DEE1849D8>, 'SetDefinedOnElements': <instancemethod SetDefinedOnElements at 0x0000026DEE184A38>, 'SetIntegrationRule': <instancemethod SetIntegrationRule at 0x0000026DEE184A98>, 'CalcElementMatrix': <instancemethod CalcElementMatrix at 0x0000026DEE184AF8>})
    pass
class BSpline():
    """
    BSpline of arbitrary order

    Parameters:

    order : int
      order of the BSpline

    knots : list
      list of float

    vals : list
      list of float
    """
    def Differentiate(self) -> BSpline: 
        """
        Differentiate the BSpline
        """
    def Integrate(self) -> BSpline: 
        """
        Integrate the BSpline
        """
    @overload
    def __call__(self, arg0: float) -> float: ...
    @overload
    def __call__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    def __init__(self, order: int, knots: list, vals: list) -> None: 
        """
        B-Spline of a certain order, provide knot and value vectors
        """
    def __str__(self) -> str: ...
    pass
class BaseMappedIntegrationPoint():
    def __init__(self, arg0: MeshPoint) -> None: ...
    def __str__(self) -> str: ...
    @property
    def elementid(self) -> ngfem::ElementId:
        """
        Element ID of the mapped integration point

        :type: ngfem::ElementId
        """
    @property
    def jacobi(self) -> ngsolve.bla.FlatMatrixD:
        """
        jacobian of the mapped integration point

        :type: ngsolve.bla.FlatMatrixD
        """
    @property
    def measure(self) -> float:
        """
        Measure of the mapped integration point 

        :type: float
        """
    @property
    def point(self) -> ngsolve.bla.FlatVectorD:
        """
        Point of the mapped integration point

        :type: ngsolve.bla.FlatVectorD
        """
    @property
    def trafo(self) -> ngfem::ElementTransformation:
        """
        Transformation of the mapped integration point

        :type: ngfem::ElementTransformation
        """
    pass
class CoefficientFunction():
    """
    A CoefficientFunction (CF) is some function defined on a mesh.
    Examples are coordinates x, y, z, domain-wise constants, solution-fields, ...
    CFs can be combined by mathematical operations (+,-,sin(), ...) to form new CFs
    Parameters:

    val : can be one of the following:

      scalar (float or complex):
        Creates a constant CoefficientFunction with value val

      tuple of scalars or CoefficientFunctions:
        Creates a vector or matrix valued CoefficientFunction, use dims=(h,w)
        for matrix valued CF
      list of scalars or CoefficientFunctions:
        Creates a domain-wise CF, use with generator expressions and mesh.GetMaterials()
        and mesh.GetBoundaries()
    """
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> CoefficientFunction: 
        """
        Compile list of individual steps, experimental improvement for deep trees

        Parameters:

        realcompile : bool
          True -> Compile to C++ code

        maxderiv : int
          input maximal derivative

        wait : bool
          True -> Waits until the previous Compile call is finished before start compiling
        """
    def Derive(self, variable: CoefficientFunction, direction: CoefficientFunction = 1.0) -> CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: CoefficientFunction, direction: CoefficientFunction = None) -> CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: CoefficientFunction = 1.0) -> CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def InnerProduct(self, cf: CoefficientFunction) -> CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Norm(self) -> CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Other(self) -> CoefficientFunction: 
        """
        Evaluate on other element, as needed for DG jumps
        """
    @overload
    def __add__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __add__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __call__(self, mip: BaseMappedIntegrationPoint) -> object: 
        """
        evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)
        """
    @overload
    def __call__(self, arg0: numpy.ndarray[MeshPoint]) -> array: ...
    @overload
    def __getitem__(self, comp: int) -> CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, components: tuple) -> CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __init__(self, arg0: dict) -> None: 
        """
        Construct a CoefficientFunction from either one of
          a scalar (float or complex)
          a tuple of scalars and or CFs to define a vector-valued CF
             use dims=(h,w) to define matrix-valued CF
          a list of scalars and or CFs to define a domain-wise CF
        """
    @overload
    def __init__(self, coef: object, dims: Optional[tuple] = None) -> None: ...
    @overload
    def __mul__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    @overload
    def __mul__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __mul__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    def __neg__(self) -> CoefficientFunction: ...
    @overload
    def __pow__(self, exponent: int) -> CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, arg0: CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> CoefficientFunction: ...
    def __rsub__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    @overload
    def __sub__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def dim(self) -> int:
        """
        number of components of CF

        :type: int
        """
    @property
    def dims(self) -> pyngcore.Array_I_S:
        """
        shape of CF:  (dim) for vector, (h,w) for matrix

        :type: pyngcore.Array_I_S
        """
    @dims.setter
    def dims(self, arg1: tuple) -> None:
        """
        shape of CF:  (dim) for vector, (h,w) for matrix
        """
    @property
    def imag(self) -> CoefficientFunction:
        """
        imaginary part of CF

        :type: CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def real(self) -> CoefficientFunction:
        """
        real part of CF

        :type: CoefficientFunction
        """
    @property
    def trans(self) -> CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: CoefficientFunction
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE18B498>, '__dict__': <attribute '__dict__' of 'ngsolve.fem.CoefficientFunction' objects>, '__doc__': 'A CoefficientFunction (CF) is some function defined on a mesh.\nExamples are coordinates x, y, z, domain-wise constants, solution-fields, ...\nCFs can be combined by mathematical operations (+,-,sin(), ...) to form new CFs\nParameters:\n\nval : can be one of the following:\n\n  scalar (float or complex):\n    Creates a constant CoefficientFunction with value val\n\n  tuple of scalars or CoefficientFunctions:\n    Creates a vector or matrix valued CoefficientFunction, use dims=(h,w)\n    for matrix valued CF\n  list of scalars or CoefficientFunctions:\n    Creates a domain-wise CF, use with generator expressions and mesh.GetMaterials()\n    and mesh.GetBoundaries()\n', '__module__': 'ngsolve.fem', '__str__': <instancemethod __str__ at 0x0000026DEE18B4C8>, '__call__': <instancemethod __call__ at 0x0000026DEE18F108>, 'dim': <property object at 0x0000026DEE18C4F8>, 'dims': <property object at 0x0000026DEE18C5E8>, 'is_complex': <property object at 0x0000026DEE18C688>, '__getitem__': <instancemethod __getitem__ at 0x0000026DEE18B678>, '__add__': <instancemethod __add__ at 0x0000026DEE18B6A8>, '__radd__': <instancemethod __radd__ at 0x0000026DEE18B738>, '__sub__': <instancemethod __sub__ at 0x0000026DEE18B798>, '__rsub__': <instancemethod __rsub__ at 0x0000026DEE18B7C8>, '__mul__': <instancemethod __mul__ at 0x0000026DEE18BB58>, '__pow__': <instancemethod __pow__ at 0x0000026DEE18B888>, 'InnerProduct': <instancemethod InnerProduct at 0x0000026DEE18B8E8>, 'Norm': <instancemethod Norm at 0x0000026DEE18B948>, 'Eig': <instancemethod Eig at 0x0000026DEE18B9A8>, 'Other': <instancemethod Other at 0x0000026DEE18BA08>, 'Derive': <instancemethod Derive at 0x0000026DEE18BA68>, 'Diff': <instancemethod Diff at 0x0000026DEE18BAC8>, 'DiffShape': <instancemethod DiffShape at 0x0000026DEE18BB28>, '__rmul__': <instancemethod __rmul__ at 0x0000026DEE18BBB8>, '__truediv__': <instancemethod __truediv__ at 0x0000026DEE18BBE8>, '__rtruediv__': <instancemethod __rtruediv__ at 0x0000026DEE18BC78>, '__neg__': <instancemethod __neg__ at 0x0000026DEE18BCA8>, 'trans': <property object at 0x0000026DEE18CD18>, 'real': <property object at 0x0000026DEE18CDB8>, 'imag': <property object at 0x0000026DEE18CEA8>, 'Freeze': <instancemethod Freeze at 0x0000026DEE18BD98>, 'Compile': <instancemethod Compile at 0x0000026DEE18BDF8>, 'data': <property object at 0x0000026DEE18D048>, '__getstate__': <instancemethod __getstate__ at 0x0000026DEE18BE88>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE18BEE8>})
    pass
class DifferentialOperator():
    def __timing__(self, arg0: FiniteElement, arg1: ElementTransformation, arg2: IntegrationRule) -> List[Tuple[str, float]]: ...
    pass
class ET():
    """
    Enumeration of all supported element types.

    Members:

      POINT

      SEGM

      TRIG

      QUAD

      TET

      PRISM

      PYRAMID

      HEX
    """
    def __init__(self, arg0: int) -> None: ...
    def __int__(self) -> int: ...
    @property
    def name(self) -> str:
        """
        (self: handle) -> str

        :type: str
        """
    HEX: ngsolve.fem.ET # value = ET.HEX
    POINT: ngsolve.fem.ET # value = ET.POINT
    PRISM: ngsolve.fem.ET # value = ET.PRISM
    PYRAMID: ngsolve.fem.ET # value = ET.PYRAMID
    QUAD: ngsolve.fem.ET # value = ET.QUAD
    SEGM: ngsolve.fem.ET # value = ET.SEGM
    TET: ngsolve.fem.ET # value = ET.TET
    TRIG: ngsolve.fem.ET # value = ET.TRIG
    __entries: dict # value = {'POINT': (ET.POINT, None), 'SEGM': (ET.SEGM, None), 'TRIG': (ET.TRIG, None), 'QUAD': (ET.QUAD, None), 'TET': (ET.TET, None), 'PRISM': (ET.PRISM, None), 'PYRAMID': (ET.PYRAMID, None), 'HEX': (ET.HEX, None)}
    __members__: dict # value = {'POINT': ET.POINT, 'SEGM': ET.SEGM, 'TRIG': ET.TRIG, 'QUAD': ET.QUAD, 'TET': ET.TET, 'PRISM': ET.PRISM, 'PYRAMID': ET.PYRAMID, 'HEX': ET.HEX}
    pass
class ElementTopology():
    """
    Element Topology

    Parameters:

    et : ngsolve.fem.ET
      input element type
    """
    def __init__(self, et: ET) -> None: ...
    @property
    def name(self) -> str:
        """
        Name of the element topology

        :type: str
        """
    @property
    def vertices(self) -> list:
        """
        Vertices of the element topology

        :type: list
        """
    pass
class ElementTransformation():
    @overload
    def __call__(self, ip: IntegrationPoint) -> BaseMappedIntegrationPoint: ...
    @overload
    def __call__(self, arg0: IntegrationRule) -> numpy.ndarray[MeshPoint]: ...
    @overload
    def __call__(self, x: float, y: float = 0, z: float = 0) -> BaseMappedIntegrationPoint: ...
    def __init__(self, et: ET = ET.TRIG, vertices: list) -> None: ...
    @property
    def VB(self) -> ngfem::VorB:
        """
        :type: ngfem::VorB
        """
    @property
    def curved(self) -> bool:
        """
        Is mapping non-affine ?

        :type: bool
        """
    @property
    def elementid(self) -> ngfem::ElementId:
        """
        Element ID of the element transformation

        :type: ngfem::ElementId
        """
    @property
    def spacedim(self) -> int:
        """
        Space dimension of the element transformation

        :type: int
        """
    pass
class FiniteElement():
    """
    any finite element
    """
    def __str__(self) -> str: ...
    def __timing__(self) -> List[Tuple[str, float]]: ...
    @property
    def classname(self) -> str:
        """
        name of element family

        :type: str
        """
    @property
    def dim(self) -> int:
        """
        spatial dimension of element

        :type: int
        """
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom of element

        :type: int
        """
    @property
    def order(self) -> int:
        """
        maximal polynomial order of element

        :type: int
        """
    @property
    def type(self) -> ET:
        """
        geometric type of element

        :type: ET
        """
    pass
class HCurlFE(FiniteElement):
    """
    an H(curl) finite element
    """
    def CalcCurlShape(self, mip: ngfem::BaseMappedIntegrationPoint) -> ngsolve.bla.MatrixD: ...
    @overload
    def CalcShape(self, x: float, y: float = 0.0, z: float = 0.0) -> ngsolve.bla.MatrixD: ...
    @overload
    def CalcShape(self, mip: ngfem::BaseMappedIntegrationPoint) -> ngsolve.bla.MatrixD: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> List[Tuple[str, float]]: ...
    @property
    def classname(self) -> str:
        """
        name of element family

        :type: str
        """
    @property
    def dim(self) -> int:
        """
        spatial dimension of element

        :type: int
        """
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom of element

        :type: int
        """
    @property
    def order(self) -> int:
        """
        maximal polynomial order of element

        :type: int
        """
    @property
    def type(self) -> ET:
        """
        geometric type of element

        :type: ET
        """
    pass
class HDivDivFE(FiniteElement):
    """
    an H(div div) finite element
    """
    def CalcDivShape(self, x: float, y: float = 0.0, z: float = 0.0) -> ngsolve.bla.MatrixD: ...
    def CalcShape(self, x: float, y: float = 0.0, z: float = 0.0) -> ngsolve.bla.MatrixD: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> List[Tuple[str, float]]: ...
    @property
    def classname(self) -> str:
        """
        name of element family

        :type: str
        """
    @property
    def dim(self) -> int:
        """
        spatial dimension of element

        :type: int
        """
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom of element

        :type: int
        """
    @property
    def order(self) -> int:
        """
        maximal polynomial order of element

        :type: int
        """
    @property
    def type(self) -> ET:
        """
        geometric type of element

        :type: ET
        """
    pass
class HDivFE(FiniteElement):
    """
    an H(div) finite element
    """
    def CalcDivShape(self, x: float, y: float = 0.0, z: float = 0.0) -> ngsolve.bla.VectorD: ...
    def CalcShape(self, x: float, y: float = 0.0, z: float = 0.0) -> ngsolve.bla.MatrixD: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> List[Tuple[str, float]]: ...
    @property
    def classname(self) -> str:
        """
        name of element family

        :type: str
        """
    @property
    def dim(self) -> int:
        """
        spatial dimension of element

        :type: int
        """
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom of element

        :type: int
        """
    @property
    def order(self) -> int:
        """
        maximal polynomial order of element

        :type: int
        """
    @property
    def type(self) -> ET:
        """
        geometric type of element

        :type: ET
        """
    pass
class IntegrationPoint():
    @property
    def point(self) -> tuple:
        """
        Integration point coordinates as tuple, has always x,y and z component, which do not have meaning in lesser dimensions

        :type: tuple
        """
    @property
    def weight(self) -> float:
        """
        Weight of the integration point

        :type: float
        """
    pass
class IntegrationRule():
    """
    Integration rule

    2 __init__ overloads


    1)

    Parameters:

    element type : ngsolve.fem.ET
      input element type

    order : int
      input order of integration rule


    2)

    Parameters:

    points : list
      input list of integration points

    weights : list
      input list of integration weights
    """
    def Integrate(self, func: object) -> object: 
        """
        Integrates a given function
        """
    def __getitem__(self, nr: int) -> IntegrationPoint: 
        """
        Return integration point at given position
        """
    @overload
    def __init__(self, element type: ET, order: int) -> None: ...
    @overload
    def __init__(self, points: list, weights: list) -> None: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    @property
    def points(self) -> list:
        """
        Points of IntegrationRule as tuple

        :type: list
        """
    @property
    def weights(self) -> list:
        """
        Weights of IntegrationRule

        :type: list
        """
    pass
class LFI():
    """
    Linear Form Integrator

    Parameters:

    name : string
      Name of the linear form integrator.

    dim : int
      dimension of the linear form integrator

    coef : object
      CoefficientFunction of the bilinear form.

    definedon : object
      input region where the linear form is defined on

    imag : bool
      Multiplies LFI with 1J

    flags : ngsolve.ngstd.Flags
      input flags

    definedonelem : object
      input definedonelem
    """
    @overload
    def CalcElementVector(self, fel: FiniteElement, trafo: ElementTransformation, vec: ngsolve.bla.FlatVectorD, lh: ngsolve.ngstd.LocalHeap) -> None: ...
    @overload
    def CalcElementVector(self, fel: FiniteElement, trafo: ElementTransformation, heapsize: int = 10000, complex: bool = False) -> object: ...
    def GetDefinedOn(self) -> pyngcore.BitArray: 
        """
        Reterns regions where the lienar form integrator is defined on.
        """
    def SetDefinedOnElements(self, ba: pyngcore.BitArray) -> None: 
        """
        Set the elements on which the linear form integrator is defined on

        Parameters:

        ba : ngsolve.ngstd.BitArray
          input bit array ( 1-> defined on, 0 -> not defoned on)
        """
    def SetIntegrationRule(self, et: ET, ir: IntegrationRule) -> LFI: 
        """
        Set a different integration rule for elements of type et

        Parameters:

        et : ngsolve.fem.ET
          input element type

        ir : ngsolve.fem.IntegrationRule
          input integration rule
        """
    def __init__(self, name: str = 0, dim: int = -1, coef: object, definedon: object = <ngsolve.ngstd.DummyArgument>, imag: bool = False, flags: pyngcore.Flags = {}, definedonelements: object = <ngsolve.ngstd.DummyArgument>) -> None: ...
    def __str__(self) -> str: ...
    @property
    def simd_evaluate(self) -> bool:
        """
        SIMD evaluate ?

        :type: bool
        """
    @simd_evaluate.setter
    def simd_evaluate(self, arg1: bool) -> None:
        """
        SIMD evaluate ?
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE184BB8>, '__dict__': <attribute '__dict__' of 'ngsolve.fem.LFI' objects>, '__doc__': '\nLinear Form Integrator\n\nParameters:\n\nname : string\n  Name of the linear form integrator.\n\ndim : int\n  dimension of the linear form integrator\n\ncoef : object\n  CoefficientFunction of the bilinear form.\n\ndefinedon : object\n  input region where the linear form is defined on\n\nimag : bool\n  Multiplies LFI with 1J\n\nflags : ngsolve.ngstd.Flags\n  input flags\n\ndefinedonelem : object\n  input definedonelem\n\n', '__module__': 'ngsolve.fem', '__str__': <instancemethod __str__ at 0x0000026DEE184C18>, 'simd_evaluate': <property object at 0x0000026DEE189548>, 'GetDefinedOn': <instancemethod GetDefinedOn at 0x0000026DEE184CD8>, 'SetDefinedOnElements': <instancemethod SetDefinedOnElements at 0x0000026DEE184D38>, 'SetIntegrationRule': <instancemethod SetIntegrationRule at 0x0000026DEE184D98>, 'CalcElementVector': <instancemethod CalcElementVector at 0x0000026DEE184E28>})
    pass
class MeshPoint():
    @property
    def mesh(self) -> ngcomp::MeshAccess:
        """
        :type: ngcomp::MeshAccess
        """
    @property
    def nr(self) -> int:
        """
        :type: int
        """
    @property
    def pnt(self) -> tuple:
        """
        Gives coordinates of point on reference triangle. One can create a MappedIntegrationPoint using the ngsolve.fem.BaseMappedIntegrationPoint constructor. For physical coordinates the coordinate CoefficientFunctions x,y,z can be evaluated in the MeshPoint

        :type: tuple
        """
    @property
    def vb(self) -> ngfem::VorB:
        """
        :type: ngfem::VorB
        """
    pass
class MixedFE(FiniteElement):
    """
    pair of finite elements for trial and test-functions
    """
    def __init__(self, arg0: FiniteElement, arg1: FiniteElement) -> None: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> List[Tuple[str, float]]: ...
    @property
    def classname(self) -> str:
        """
        name of element family

        :type: str
        """
    @property
    def dim(self) -> int:
        """
        spatial dimension of element

        :type: int
        """
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom of element

        :type: int
        """
    @property
    def order(self) -> int:
        """
        maximal polynomial order of element

        :type: int
        """
    @property
    def type(self) -> ET:
        """
        geometric type of element

        :type: ET
        """
    pass
class NODE_TYPE():
    """
    Enumeration of all supported node types.

    Members:

      VERTEX

      EDGE

      FACE

      CELL

      ELEMENT

      FACET
    """
    def __init__(self, arg0: int) -> None: ...
    def __int__(self) -> int: ...
    @property
    def name(self) -> str:
        """
        (self: handle) -> str

        :type: str
        """
    CELL: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.CELL
    EDGE: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.EDGE
    ELEMENT: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.ELEMENT
    FACE: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.FACE
    FACET: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.FACET
    VERTEX: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.VERTEX
    __entries: dict # value = {'VERTEX': (NODE_TYPE.VERTEX, None), 'EDGE': (NODE_TYPE.EDGE, None), 'FACE': (NODE_TYPE.FACE, None), 'CELL': (NODE_TYPE.CELL, None), 'ELEMENT': (NODE_TYPE.ELEMENT, None), 'FACET': (NODE_TYPE.FACET, None)}
    __members__: dict # value = {'VERTEX': NODE_TYPE.VERTEX, 'EDGE': NODE_TYPE.EDGE, 'FACE': NODE_TYPE.FACE, 'CELL': NODE_TYPE.CELL, 'ELEMENT': NODE_TYPE.ELEMENT, 'FACET': NODE_TYPE.FACET}
    pass
class Parameter(CoefficientFunction):
    """
    CoefficientFunction with a modifiable value

    Parameters:

    value : float
      Parameter value
    """
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> CoefficientFunction: 
        """
        Compile list of individual steps, experimental improvement for deep trees

        Parameters:

        realcompile : bool
          True -> Compile to C++ code

        maxderiv : int
          input maximal derivative

        wait : bool
          True -> Waits until the previous Compile call is finished before start compiling
        """
    def Derive(self, variable: CoefficientFunction, direction: CoefficientFunction = 1.0) -> CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: CoefficientFunction, direction: CoefficientFunction = None) -> CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: CoefficientFunction = 1.0) -> CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def Get(self) -> float: 
        """
        return parameter value
        """
    def InnerProduct(self, cf: CoefficientFunction) -> CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Norm(self) -> CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Other(self) -> CoefficientFunction: 
        """
        Evaluate on other element, as needed for DG jumps
        """
    def Set(self, value: float) -> None: 
        """
        Modify parameter value.

        Parameters:

        value : double
          input scalar  
        """
    @overload
    def __add__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __add__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __call__(self, mip: BaseMappedIntegrationPoint) -> object: 
        """
        evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)
        """
    @overload
    def __call__(self, arg0: numpy.ndarray[MeshPoint]) -> array: ...
    @overload
    def __getitem__(self, comp: int) -> CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, components: tuple) -> CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, value: float) -> None: 
        """
        Construct a ParameterCF from a scalar
        """
    @overload
    def __mul__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    @overload
    def __mul__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __mul__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    def __neg__(self) -> CoefficientFunction: ...
    @overload
    def __pow__(self, exponent: int) -> CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, arg0: CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> CoefficientFunction: ...
    def __rsub__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    @overload
    def __sub__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def dim(self) -> int:
        """
        number of components of CF

        :type: int
        """
    @property
    def dims(self) -> pyngcore.Array_I_S:
        """
        shape of CF:  (dim) for vector, (h,w) for matrix

        :type: pyngcore.Array_I_S
        """
    @dims.setter
    def dims(self, arg1: tuple) -> None:
        """
        shape of CF:  (dim) for vector, (h,w) for matrix
        """
    @property
    def imag(self) -> CoefficientFunction:
        """
        imaginary part of CF

        :type: CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def real(self) -> CoefficientFunction:
        """
        real part of CF

        :type: CoefficientFunction
        """
    @property
    def trans(self) -> CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: CoefficientFunction
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE18F138>, '__dict__': <attribute '__dict__' of 'ngsolve.fem.Parameter' objects>, '__doc__': '\nCoefficientFunction with a modifiable value\n\nParameters:\n\nvalue : float\n  Parameter value\n\n', '__module__': 'ngsolve.fem', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE18F198>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE18F1F8>, 'Set': <instancemethod Set at 0x0000026DEE18F258>, 'Get': <instancemethod Get at 0x0000026DEE18F2B8>})
    pass
class ParameterC(CoefficientFunction):
    """
    CoefficientFunction with a modifiable complex value

    Parameters:

    value : complex
      Parameter value
    """
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> CoefficientFunction: 
        """
        Compile list of individual steps, experimental improvement for deep trees

        Parameters:

        realcompile : bool
          True -> Compile to C++ code

        maxderiv : int
          input maximal derivative

        wait : bool
          True -> Waits until the previous Compile call is finished before start compiling
        """
    def Derive(self, variable: CoefficientFunction, direction: CoefficientFunction = 1.0) -> CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: CoefficientFunction, direction: CoefficientFunction = None) -> CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: CoefficientFunction = 1.0) -> CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def Get(self) -> float: 
        """
        return parameter value
        """
    def InnerProduct(self, cf: CoefficientFunction) -> CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Norm(self) -> CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Other(self) -> CoefficientFunction: 
        """
        Evaluate on other element, as needed for DG jumps
        """
    def Set(self, value: complex) -> None: 
        """
        Modify parameter value.

        Parameters:

        value : complex
          input scalar
        """
    @overload
    def __add__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __add__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __call__(self, mip: BaseMappedIntegrationPoint) -> object: 
        """
        evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)
        """
    @overload
    def __call__(self, arg0: numpy.ndarray[MeshPoint]) -> array: ...
    @overload
    def __getitem__(self, comp: int) -> CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, components: tuple) -> CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, value: complex) -> None: 
        """
        Construct a ParameterCF from a scalar
        """
    @overload
    def __mul__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    @overload
    def __mul__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __mul__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    def __neg__(self) -> CoefficientFunction: ...
    @overload
    def __pow__(self, exponent: int) -> CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, arg0: CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> CoefficientFunction: ...
    def __rsub__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    @overload
    def __sub__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: CoefficientFunction) -> CoefficientFunction: ...
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def dim(self) -> int:
        """
        number of components of CF

        :type: int
        """
    @property
    def dims(self) -> pyngcore.Array_I_S:
        """
        shape of CF:  (dim) for vector, (h,w) for matrix

        :type: pyngcore.Array_I_S
        """
    @dims.setter
    def dims(self, arg1: tuple) -> None:
        """
        shape of CF:  (dim) for vector, (h,w) for matrix
        """
    @property
    def imag(self) -> CoefficientFunction:
        """
        imaginary part of CF

        :type: CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def real(self) -> CoefficientFunction:
        """
        real part of CF

        :type: CoefficientFunction
        """
    @property
    def trans(self) -> CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: CoefficientFunction
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE18F318>, '__dict__': <attribute '__dict__' of 'ngsolve.fem.ParameterC' objects>, '__doc__': '\nCoefficientFunction with a modifiable complex value\n\nParameters:\n\nvalue : complex\n  Parameter value\n\n', '__module__': 'ngsolve.fem', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE18F378>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE18F3D8>, 'Set': <instancemethod Set at 0x0000026DEE18F438>, 'Get': <instancemethod Get at 0x0000026DEE18F498>})
    pass
class ScalarFE(FiniteElement):
    """
    a scalar-valued finite element
    """
    @overload
    def CalcDShape(self, mip: ngfem::BaseMappedIntegrationPoint) -> ngsolve.bla.MatrixD: 
        """
        Computes derivative of the shape in an integration point.

        Parameters:

        mip : ngsolve.BaseMappedIntegrationPoint
          input mapped integration point




        Computes derivative of the shape in an integration point.

        Parameters:

        x : double
          input x value

        y : double
          input y value

        z : double
          input z value
        """
    @overload
    def CalcDShape(self, x: float, y: float = 0.0, z: float = 0.0) -> ngsolve.bla.MatrixD: ...
    @overload
    def CalcShape(self, x: float, y: float = 0.0, z: float = 0.0) -> ngsolve.bla.VectorD: 
        """
        Parameters:

        x : double
          input x value

        y : double
          input y value

        z : double
          input z value




        Parameters:

        mip : ngsolve.BaseMappedIntegrationPoint
          input mapped integration point
        """
    @overload
    def CalcShape(self, mip: ngfem::BaseMappedIntegrationPoint) -> ngsolve.bla.VectorD: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> List[Tuple[str, float]]: ...
    @property
    def classname(self) -> str:
        """
        name of element family

        :type: str
        """
    @property
    def dim(self) -> int:
        """
        spatial dimension of element

        :type: int
        """
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom of element

        :type: int
        """
    @property
    def order(self) -> int:
        """
        maximal polynomial order of element

        :type: int
        """
    @property
    def type(self) -> ET:
        """
        geometric type of element

        :type: ET
        """
    pass
class SpecialCFCreator():
    def JacobianMatrix(self, dim: int) -> ngfem::CoefficientFunction: 
        """
        Jacobian matrix of transformation to physical element
        space-dimension must be provided
        """
    def Weingarten(self, dim: int) -> ngfem::CoefficientFunction: 
        """
        Weingarten tensor 
        space-dimension must be provided
        """
    def normal(self, dim: int) -> ngfem::CoefficientFunction: 
        """
        depending on contents: normal-vector to geometry or element
        space-dimension must be provided
        """
    def tangential(self, dim: int) -> ngfem::CoefficientFunction: 
        """
        depending on contents: tangential-vector to element
        space-dimension must be provided
        """
    def xref(self, dim: int) -> ngfem::CoefficientFunction: 
        """
        element reference-coordinates
        """
    @property
    def mesh_size(self) -> ngfem::CoefficientFunction:
        """
        local mesh-size (approximate element diameter) as CF

        :type: ngfem::CoefficientFunction
        """
    pass
def  (x: object) -> object:
    """
     (x: object) -> object

    Passes value through
    """
def BlockBFI(bfi: BFI = 0, dim: int = 2, comp: int = 0) -> ngfem::BlockBilinearFormIntegrator:
    """
    Block Bilinear Form Integrator

    Parameters:

    bfi : ngsolve.fem.BFI
      input bilinear form integrator

    dim : int
      input dimension of block bilinear form integrator

    comp : int
      input comp
    """
def BlockLFI(lfi: LFI = 0, dim: int = 2, comp: int = 0) -> LFI:
    """
    Block Linear Form Integrator

    Parameters:

    lfi : ngsolve.fem.LFI
      input bilinear form integrator

    dim : int
      input dimension of block linear form integrator

    comp : int
      input comp
    """
def CacheCF(cf: CoefficientFunction) -> CoefficientFunction:
    pass
def Cof(arg0: CoefficientFunction) -> CoefficientFunction:
    pass
def CompoundBFI(bfi: BFI = 0, comp: int = 0) -> ngfem::CompoundBilinearFormIntegrator:
    """
    Compound Bilinear Form Integrator

    Parameters:

    bfi : ngsolve.fem.BFI
      input bilinear form integrator

    comp : int
      input component
    """
def CompoundLFI(lfi: LFI = 0, comp: int = 0) -> LFI:
    """
    Compound Linear Form Integrator

    Parameters:

    lfi : ngsolve.fem.LFI
      input linear form integrator

    comp : int
      input component
    """
def Conj(arg0: CoefficientFunction) -> CoefficientFunction:
    """
    complex-conjugate
    """
def CoordCF(direction: int) -> ngfem::CoefficientFunction:
    """
    CoefficientFunction for x, y, z.

    Parameters:

    direction : int
      input direction
    """
def Cross(arg0: CoefficientFunction, arg1: CoefficientFunction) -> CoefficientFunction:
    pass
def Det(arg0: CoefficientFunction) -> CoefficientFunction:
    pass
def GenerateL2ElementCode(arg0: int) -> str:
    pass
def H1FE(et: ET, order: int) -> None:
    """
    Creates an H1 finite element of given geometric shape and polynomial order.

    Parameters:

    et : ngsolve.fem.ET
      input element type

    order : int
      input polynomial order
    """
def Id(arg0: int) -> CoefficientFunction:
    """
    Identity matrix of given dimension
    """
def IfPos(c1: ngfem::CoefficientFunction, then_obj: object, else_obj: object) -> ngfem::CoefficientFunction:
    """
    Returns new CoefficientFunction with values then_obj if c1 is positive and else_obj else.

    Parameters:

    c1 : ngsolve.CoefficientFunction
      Indicator function

    then_obj : object
      Values of new CF if c1 is positive, object must be implicitly convertible to
      ngsolve.CoefficientFunction. See help(CoefficientFunction ) for information.

    else_obj : object
      Values of new CF if c1 is not positive, object must be implicitly convertible to
      ngsolve.CoefficientFunction. See help(CoefficientFunction ) for information.
    """
def Inv(arg0: CoefficientFunction) -> CoefficientFunction:
    pass
def L2FE(et: ET, order: int) -> None:
    """
    Creates an L2 finite element of given geometric shape and polynomial order.

    Parameters:

    et : ngsolve.fem.ET
      input element type

    order : int
      input polynomial order
    """
def LoggingCF(cf: CoefficientFunction, logfile: str = 'stdout') -> CoefficientFunction:
    pass
def SetPMLParameters(rad: float = 1, alpha: float = 1) -> None:
    """
    Parameters:

    rad : double
      input radius of PML

    alpha : double
      input damping factor of PML
    """
def Skew(arg0: CoefficientFunction) -> CoefficientFunction:
    pass
def Sym(arg0: CoefficientFunction) -> CoefficientFunction:
    pass
def Trace(arg0: CoefficientFunction) -> CoefficientFunction:
    pass
def VoxelCoefficient(start: tuple, end: tuple, values: array, linear: bool = True, trafocf: object = <ngsolve.ngstd.DummyArgument>) -> CoefficientFunction:
    """
    CoefficientFunction defined on a grid.

    Start and end mark the cartesian boundary of domain. The function will be continued by a constant function outside of this box. Inside a cartesian grid will be created by the dimensions of the numpy input array 'values'. This array must have the dimensions of the mesh and the values stored as:
    x1y1z1, x2y1z1, ..., xNy1z1, x1y2z1, ...

    If linear is True the function will be interpolated linearly between the values. Otherwise the nearest voxel value is taken.
    """
def acos(x: object) -> object:
    """
    Inverse cosine in radians
    """
def asin(x: object) -> object:
    """
    Inverse sine in radians
    """
def atan(x: object) -> object:
    """
    Inverse tangent in radians
    """
def atan2(y: object, x: object) -> object:
    """
    Four quadrant inverse tangent in radians
    """
def ceil(x: object) -> object:
    """
    Round to next greater integer
    """
def cos(x: object) -> object:
    """
    Cosine of argument in radians
    """
def cosh(x: object) -> object:
    """
    Hyperbolic cosine of argument in radians
    """
def exp(x: object) -> object:
    """
    Exponential function
    """
def floor(x: object) -> object:
    """
    Round to next lower integer
    """
def log(x: object) -> object:
    """
    Logarithm function
    """
def pow(x: object, y: object) -> object:
    """
    Power function
    """
def sin(x: object) -> object:
    """
    Sine of argument in radians
    """
def sinh(x: object) -> object:
    """
    Hyperbolic sine of argument in radians
    """
def sqrt(x: object) -> object:
    """
    Square root function
    """
def tan(x: object) -> object:
    """
    Tangent of argument in radians
    """
CELL: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.CELL
EDGE: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.EDGE
ELEMENT: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.ELEMENT
FACE: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.FACE
FACET: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.FACET
HEX: ngsolve.fem.ET # value = ET.HEX
POINT: ngsolve.fem.ET # value = ET.POINT
PRISM: ngsolve.fem.ET # value = ET.PRISM
PYRAMID: ngsolve.fem.ET # value = ET.PYRAMID
QUAD: ngsolve.fem.ET # value = ET.QUAD
SEGM: ngsolve.fem.ET # value = ET.SEGM
TET: ngsolve.fem.ET # value = ET.TET
TRIG: ngsolve.fem.ET # value = ET.TRIG
VERTEX: ngsolve.fem.NODE_TYPE # value = NODE_TYPE.VERTEX
specialcf: ngsolve.fem.SpecialCFCreator # value = <ngsolve.fem.SpecialCFCreator object at 0x0000026DEE186470>
