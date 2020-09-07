import ngsolve.krylovspace
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import i
import n
import g
import typing
import pyngcore
import ngsolve.ngstd
import ngsolve.bla
import ngsolve.la
import l
import o
__all__  = [
"BaseMatrix",
"BaseVector",
"BitArray",
"CGSolver",
"Preconditioner",
"Projector",
"CG",
"GMRes",
"InnerProduct",
"Matrix",
"MinRes",
"Norm",
"PreconditionedRichardson",
"QMR",
"TimeFunction",
"Vector",
"_GetStatus",
"_PushStatus",
"_SetThreadPercentage",
"log",
"sqrt",
"Callable",
"Optional"
]
class BaseMatrix():
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def CreateColVector(self) -> BaseVector: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def GetInverseType(self) -> str: ...
    def Inverse(self, freedofs: pyngcore.BitArray = None, inverse: str = '') -> BaseMatrix: 
        """
        Calculate inverse of sparse matrix
        Parameters:

        freedofs : BitArray
          If set, invert only the rows/columns the matrix defined by the bit array, otherwise invert the whole matrix

        inverse : string
          Solver to use, allowed values are:
            sparsecholesky - internal solver of NGSolve for symmetric matrices
            umfpack        - solver by Suitesparse/UMFPACK (if NGSolve was configured with USE_UMFPACK=ON)
            pardiso        - PARDISO, either provided by libpardiso (USE_PARDISO=ON) or Intel MKL (USE_MKL=ON).
                             If neither Pardiso nor Intel MKL was linked at compile-time, NGSolve will look
                             for libmkl_rt in LD_LIBRARY_PATH (Unix) or PATH (Windows) at run-time.
        """
    def Mult(self, x: BaseVector, y: BaseVector) -> None: ...
    @overload
    def MultAdd(self, value: float, x: BaseVector, y: BaseVector) -> None: ...
    @overload
    def MultAdd(self, value: complex, x: BaseVector, y: BaseVector) -> None: ...
    @overload
    def MultScale(self, value: float, x: BaseVector, y: BaseVector) -> None: ...
    @overload
    def MultScale(self, value: complex, x: BaseVector, y: BaseVector) -> None: ...
    @overload
    def MultTrans(self, value: complex, x: BaseVector, y: BaseVector) -> None: ...
    @overload
    def MultTrans(self, value: float, x: BaseVector, y: BaseVector) -> None: ...
    @overload
    def MultTransAdd(self, value: float, x: BaseVector, y: BaseVector) -> None: ...
    @overload
    def MultTransAdd(self, value: complex, x: BaseVector, y: BaseVector) -> None: ...
    def Update(self) -> None: 
        """
        Update matrix
        """
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    def __init__(self) -> None: ...
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __str__(self) -> str: ...
    def __sub__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @property
    def H(self) -> BaseMatrix:
        """
        Return conjugate transpose of matrix (WIP, only partially supported)

        :type: BaseMatrix
        """
    @property
    def T(self) -> BaseMatrix:
        """
        Return transpose of matrix

        :type: BaseMatrix
        """
    @property
    def height(self) -> int:
        """
        Height of the matrix

        :type: int
        """
    @property
    def local_mat(self) -> BaseMatrix:
        """
        :type: BaseMatrix
        """
    @property
    def nze(self) -> int:
        """
        number of non-zero elements

        :type: int
        """
    @property
    def width(self) -> int:
        """
        Width of the matrix

        :type: int
        """
    pass
class BaseVector():
    def Add(self, vec: BaseVector, value: object) -> None: ...
    def Assign(self, vec: BaseVector, value: object) -> None: ...
    def Copy(self) -> BaseVector: 
        """
        creates a new vector of same type, copy contents
        """
    def CreateVector(self) -> BaseVector: 
        """
        creates a new vector of same type, contents is undefined
        """
    def Cumulate(self) -> None: ...
    def Distribute(self) -> None: ...
    def FV(self) -> object: ...
    def GetParallelStatus(self) -> PARALLEL_STATUS: ...
    def InnerProduct(self, other: BaseVector, conjugate: bool = True) -> object: ...
    def Norm(self) -> float: 
        """
        Calculate Norm
        """
    def Range(self, from: int, to: int) -> BaseVector: 
        """
        Return values from given range
        """
    def Reshape(self, width: int) -> ngsolve.bla.FlatMatrixD: ...
    def SetParallelStatus(self, stat: PARALLEL_STATUS) -> None: ...
    def SetRandom(self) -> None: ...
    def __add__(self, arg0: ngla::DynamicVectorExpression) -> ngla::DynamicVectorExpression: ...
    @overload
    def __getitem__(self, arg0: ngsolve.ngstd.IntRange) -> BaseVector: 
        """
        Return value at given position

        Return values at given position
        """
    @overload
    def __getitem__(self, ind: int) -> object: ...
    @overload
    def __getitem__(self, inds: slice) -> BaseVector: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __iadd__(self, vec: BaseVector) -> BaseVector: ...
    @overload
    def __iadd__(self, arg0: ngla::DynamicVectorExpression) -> BaseVector: ...
    @overload
    def __imul__(self, value: complex) -> BaseVector: ...
    @overload
    def __imul__(self, value: float) -> BaseVector: ...
    def __init__(self, size: int, complex: bool = False, entrysize: int = 1) -> None: ...
    @overload
    def __isub__(self, vec: BaseVector) -> BaseVector: ...
    @overload
    def __isub__(self, arg0: ngla::DynamicVectorExpression) -> BaseVector: ...
    @overload
    def __itruediv__(self, value: float) -> BaseVector: ...
    @overload
    def __itruediv__(self, value: complex) -> BaseVector: ...
    def __len__(self) -> int: ...
    def __neg__(self) -> ngla::DynamicVectorExpression: ...
    def __repr__(self) -> str: ...
    @overload
    def __rmul__(self, arg0: complex) -> ngla::DynamicVectorExpression: ...
    @overload
    def __rmul__(self, arg0: float) -> ngla::DynamicVectorExpression: ...
    @overload
    def __setitem__(self, ind: int, value: float) -> None: 
        """
        Set value at given position

        Set value at given position

        Set value at given positions

        Set value at given positions

        Set value for range of indices

        Set value for range of indices
        """
    @overload
    def __setitem__(self, inds: slice, value: complex) -> None: ...
    @overload
    def __setitem__(self, range: ngsolve.ngstd.IntRange, value: float) -> None: ...
    @overload
    def __setitem__(self, ind: int, value: complex) -> None: ...
    @overload
    def __setitem__(self, range: ngsolve.ngstd.IntRange, vec: BaseVector) -> None: ...
    @overload
    def __setitem__(self, inds: slice, value: float) -> None: ...
    @overload
    def __setitem__(self, range: ngsolve.ngstd.IntRange, value: complex) -> None: ...
    @overload
    def __setitem__(self, ind: int, vec: ngsolve.bla.FlatVectorD) -> None: ...
    @overload
    def __setitem__(self, ind: int, vec: ngsolve.bla.FlatVectorC) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    def __sub__(self, arg0: ngla::DynamicVectorExpression) -> ngla::DynamicVectorExpression: ...
    @property
    def data(self) -> BaseVector:
        """
        :type: BaseVector
        """
    @data.setter
    def data(self, arg1: ngla::DynamicVectorExpression) -> None:
        pass
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def local_vec(self) -> BaseVector:
        """
        :type: BaseVector
        """
    @property
    def size(self) -> int:
        """
        :type: int
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE16C9A8>, '__dict__': <attribute '__dict__' of 'ngsolve.la.BaseVector' objects>, '__doc__': None, '__module__': 'ngsolve.la', 'local_vec': <property object at 0x0000026DEE16F908>, '__getstate__': <instancemethod __getstate__ at 0x0000026DEE16CA38>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE16CA98>, '__str__': <instancemethod __str__ at 0x0000026DEE16CAF8>, '__repr__': <instancemethod __repr__ at 0x0000026DEE16CB58>, 'size': <property object at 0x0000026DEE16FC78>, '__len__': <instancemethod __len__ at 0x0000026DEE16CBE8>, 'is_complex': <property object at 0x0000026DEE16FD68>, 'CreateVector': <instancemethod CreateVector at 0x0000026DEE16CC78>, 'Copy': <instancemethod Copy at 0x0000026DEE16CCD8>, 'Assign': <instancemethod Assign at 0x0000026DEE16CD38>, 'Add': <instancemethod Add at 0x0000026DEE16CD98>, '__getitem__': <instancemethod __getitem__ at 0x0000026DEE16CDF8>, '__setitem__': <instancemethod __setitem__ at 0x0000026DEE16CE58>, '__iadd__': <instancemethod __iadd__ at 0x0000026DEE171108>, '__isub__': <instancemethod __isub__ at 0x0000026DEE171168>, '__imul__': <instancemethod __imul__ at 0x0000026DEE16CFA8>, '__itruediv__': <instancemethod __itruediv__ at 0x0000026DEE171048>, 'data': <property object at 0x0000026DEE1701D8>, '__add__': <instancemethod __add__ at 0x0000026DEE1710D8>, '__sub__': <instancemethod __sub__ at 0x0000026DEE171138>, '__neg__': <instancemethod __neg__ at 0x0000026DEE171198>, '__rmul__': <instancemethod __rmul__ at 0x0000026DEE171228>, 'InnerProduct': <instancemethod InnerProduct at 0x0000026DEE171258>, 'Norm': <instancemethod Norm at 0x0000026DEE1712B8>, 'Range': <instancemethod Range at 0x0000026DEE171318>, 'FV': <instancemethod FV at 0x0000026DEE171378>, 'Reshape': <instancemethod Reshape at 0x0000026DEE1713D8>, 'SetRandom': <instancemethod SetRandom at 0x0000026DEE171438>, 'Distribute': <instancemethod Distribute at 0x0000026DEE171498>, 'Cumulate': <instancemethod Cumulate at 0x0000026DEE1714F8>, 'GetParallelStatus': <instancemethod GetParallelStatus at 0x0000026DEE171558>, 'SetParallelStatus': <instancemethod SetParallelStatus at 0x0000026DEE1715B8>})
    pass
class BitArray():
    @overload
    def Clear(self) -> None: 
        """
        Clear all bits

        Clear bit at given position
        """
    @overload
    def Clear(self, i: int) -> None: ...
    def NumSet(self) -> int: ...
    @overload
    def Set(self) -> None: 
        """
        Set all bits

        Set bit at given position
        """
    @overload
    def Set(self, i: int) -> None: ...
    def __and__(self, arg0: BitArray) -> BitArray: ...
    def __getitem__(self, pos: int) -> bool: 
        """
        Returns bit from given position
        """
    def __getstate__(self) -> tuple: ...
    def __iand__(self, arg0: BitArray) -> BitArray: ...
    @overload
    def __init__(self, vec: List[bool]) -> None: ...
    @overload
    def __init__(self, ba: BitArray) -> None: ...
    @overload
    def __init__(self, n: int) -> None: ...
    def __invert__(self) -> BitArray: ...
    def __ior__(self, arg0: BitArray) -> BitArray: ...
    def __len__(self) -> int: ...
    def __or__(self, arg0: BitArray) -> BitArray: ...
    @overload
    def __setitem__(self, range: ngcore::T_Range<unsigned __int64>, value: bool) -> None: 
        """
        Clear/Set bit at given position

        Clear/Set bit at given positions

        copy BitArray

        Set value for range of indices
        """
    @overload
    def __setitem__(self, inds: slice, value: bool) -> None: ...
    @overload
    def __setitem__(self, pos: int, value: bool) -> None: ...
    @overload
    def __setitem__(self, inds: slice, ba: BitArray) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    pass
class CGSolver(ngsolve.la.BaseMatrix):
    def AsVector(self) -> ngsolve.la.BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def CreateColVector(self) -> ngsolve.la.BaseVector: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> ngsolve.la.BaseVector: ...
    def GetInverseType(self) -> str: ...
    def Inverse(self, freedofs: pyngcore.BitArray = None, inverse: str = '') -> BaseMatrix: 
        """
        Calculate inverse of sparse matrix
        Parameters:

        freedofs : BitArray
          If set, invert only the rows/columns the matrix defined by the bit array, otherwise invert the whole matrix

        inverse : string
          Solver to use, allowed values are:
            sparsecholesky - internal solver of NGSolve for symmetric matrices
            umfpack        - solver by Suitesparse/UMFPACK (if NGSolve was configured with USE_UMFPACK=ON)
            pardiso        - PARDISO, either provided by libpardiso (USE_PARDISO=ON) or Intel MKL (USE_MKL=ON).
                             If neither Pardiso nor Intel MKL was linked at compile-time, NGSolve will look
                             for libmkl_rt in LD_LIBRARY_PATH (Unix) or PATH (Windows) at run-time.
        """
    @overload
    def MultAdd(self, value: float, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
    @overload
    def MultAdd(self, value: complex, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
    @overload
    def MultScale(self, value: float, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
    @overload
    def MultScale(self, value: complex, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
    @overload
    def MultTrans(self, value: float, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
    @overload
    def MultTrans(self, value: complex, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
    @overload
    def MultTransAdd(self, value: float, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
    @overload
    def MultTransAdd(self, value: complex, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    def __iadd__(self, mat: ngsolve.la.BaseMatrix) -> None: ...
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __mul__(self, arg0: ngsolve.la.MultiVector) -> ngsolve.la.MultiVectorExpr: ...
    @overload
    def __mul__(self, arg0: ngsolve.la.BaseVector) -> ngla::DynamicVectorExpression: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __str__(self) -> str: ...
    def __sub__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @property
    def H(self) -> BaseMatrix:
        """
        Return conjugate transpose of matrix (WIP, only partially supported)

        :type: BaseMatrix
        """
    @property
    def T(self) -> BaseMatrix:
        """
        Return transpose of matrix

        :type: BaseMatrix
        """
    @property
    def height(self) -> int:
        """
        Height of the matrix

        :type: int
        """
    @property
    def local_mat(self) -> BaseMatrix:
        """
        :type: BaseMatrix
        """
    @property
    def nze(self) -> int:
        """
        number of non-zero elements

        :type: int
        """
    @property
    def width(self) -> int:
        """
        Width of the matrix

        :type: int
        """
    __dict__: mappingproxy # value = mappingproxy({'__module__': 'ngsolve.krylovspace', '__init__': <function CGSolver.__init__ at 0x0000026DF1CF3F78>, 'Height': <function CGSolver.Height at 0x0000026DF1CFE048>, 'Width': <function CGSolver.Width at 0x0000026DF1CFE0D8>, 'IsComplex': <function CGSolver.IsComplex at 0x0000026DF1CFE168>, 'Mult': <function CGSolver.Mult at 0x0000026DF1CFE1F8>, 'Update': <function CGSolver.Update at 0x0000026DF1CFE288>, 'Solve': <function TimeFunction.<locals>.retfunc at 0x0000026DF1CFE3A8>, '__dict__': <attribute '__dict__' of 'CGSolver' objects>, '__doc__': None})
    pass
class Preconditioner():
    pass
class Projector():
    pass
def InnerProduct(x: object, y: object, **kwargs) -> object:
    """
    Compute InnerProduct
    """
@overload
def Matrix(arg0: List[List[complex]]) -> ngsolve.bla.MatrixC:
    """
    Creates a matrix of given height and width.

    Parameters:

    height : int
      input height

    width : int
      input width

    complex : bool
      input complex values
    """
@overload
def Matrix(arg0: List[List[float]]) -> ngsolve.bla.MatrixD:
    pass
@overload
def Matrix(height: int, width: int, complex: bool = False) -> object:
    pass
def Norm(x: object) -> object:
    """
    Compute Norm
    """
@overload
def Vector(arg0: List[float]) -> ngsolve.bla.VectorD:
    """
    Parameters:

    length : int
      input length

    complex : bool
      input complex values
    """
@overload
def Vector(length: int, complex: bool = False) -> object:
    pass
@overload
def Vector(arg0: List[complex]) -> ngsolve.bla.VectorC:
    pass
def _GetStatus() -> tuple:
    pass
def _PushStatus(arg0: str) -> None:
    pass
def _SetThreadPercentage(arg0: float) -> None:
    pass
def sqrt(x: object) -> object:
    """
    Square root function
    """
Callable: typing._VariadicGenericAlias # value = typing.Callable
Optional: typing._SpecialForm # value = typing.Optional
