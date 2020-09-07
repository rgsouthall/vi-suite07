import ngsolve.eigenvalues
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import i
import y
import c
import s
import p
import pyngcore
import ngsolve.bla
import ngsolve.la
__all__  = [
"IdentityMatrix",
"MultiVector",
"Projector",
"InnerProduct",
"Matrix",
"Norm",
"Orthogonalize",
"PINVIT",
"PINVIT1",
"Vector",
"sqrt"
]
class IdentityMatrix(BaseMatrix):
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
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(self, size: int, complex: bool = False) -> None: ...
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
class MultiVector(MultiVectorExpr):
    def Append(self, arg0: BaseVector) -> None: ...
    def AppendOrthogonalize(self, vec: BaseVector, ipmat: ngla::BaseMatrix = None) -> None: 
        """
        assumes that existing vectors are orthogonal, and orthogonalize new vector against existing vectors
        """
    def Expand(self, arg0: int) -> None: 
        """
        deprecated, use Extend instead
        """
    def Extend(self, arg0: int) -> None: ...
    @overload
    def InnerProduct(self, other: MultiVectorExpr, conjugate: bool = True) -> object: ...
    @overload
    def InnerProduct(self, other: MultiVector, conjugate: bool = True) -> object: ...
    def Orthogonalize(self, ipmat: ngla::BaseMatrix = None) -> None: ...
    def Scale(self, arg0: ngsolve.bla.VectorD) -> MultiVectorExpr: ...
    def __add__(self, arg0: MultiVectorExpr) -> MultiVectorExpr: ...
    @overload
    def __getitem__(self, arg0: int) -> BaseVector: ...
    @overload
    def __getitem__(self, arg0: slice) -> MultiVector: ...
    @overload
    def __getitem__(self, arg0: List[int]) -> MultiVector: ...
    @overload
    def __init__(self, arg0: int, arg1: int, arg2: bool) -> None: ...
    @overload
    def __init__(self, arg0: BaseVector, arg1: int) -> None: ...
    def __len__(self) -> int: ...
    @overload
    def __mul__(self, arg0: ngsolve.bla.VectorC) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: ngsolve.bla.MatrixC) -> MultiVectorExpr: ...
    @overload
    def __mul__(self, arg0: ngsolve.bla.VectorD) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: ngsolve.bla.MatrixD) -> MultiVectorExpr: ...
    def __neg__(self) -> MultiVectorExpr: ...
    @overload
    def __setitem__(self, arg0: int, arg1: float) -> None: ...
    @overload
    def __setitem__(self, arg0: List[int], arg1: MultiVector) -> MultiVector: ...
    @overload
    def __setitem__(self, arg0: List[int], arg1: MultiVectorExpr) -> None: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: MultiVectorExpr) -> None: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: MultiVector) -> None: ...
    @overload
    def __setitem__(self, arg0: int, arg1: complex) -> None: ...
    @overload
    def __setitem__(self, arg0: int, arg1: ngla::DynamicVectorExpression) -> None: ...
    def __sub__(self, arg0: MultiVectorExpr) -> MultiVectorExpr: ...
    @property
    def data(self) -> MultiVector:
        """
        :type: MultiVector
        """
    @data.setter
    def data(self, arg1: MultiVectorExpr) -> None:
        pass
    pass
class Projector(BaseMatrix):
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
    def Project(self, arg0: BaseVector) -> None: 
        """
        project vector inline
        """
    def Update(self) -> None: 
        """
        Update matrix
        """
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    def __init__(self, mask: pyngcore.BitArray, range: bool) -> None: 
        """
        Linear operator projecting to true/false bits of BitArray mask, depending on argument range
        """
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
def InnerProduct(x: object, y: object, **kwargs) -> object:
    """
    Computes InnerProduct of given objects
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
