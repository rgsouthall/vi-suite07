"""pybind la"""
import ngsolve.la
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import pyngcore
import ngsolve.bla
import ngsolve.ngstd
__all__  = [
"BaseMatrix",
"BaseSparseMatrix",
"BaseVector",
"BlockMatrix",
"BlockSmoother",
"BlockVector",
"ConstEBEMatrix",
"DynamicVectorExpression",
"Embedding",
"EmbeddingTranspose",
"IdentityMatrix",
"KrylovSpaceSolver",
"MultiVectorExpr",
"MultiVector",
"PARALLEL_STATUS",
"ParallelDofs",
"PermutationMatrix",
"Projector",
"QMRSolverC",
"QMRSolverD",
"Real2ComplexMatrix",
"S_BaseMatrixC",
"S_BaseMatrixD",
"Smoother",
"SparseFactorization",
"SparseCholesky_d",
"SparseCholesky_c",
"SparseMatrixDynamic",
"SparseMatrixclass ngbla::Mat<2,2,class std::complex<double> >",
"SparseMatrixclass ngbla::Mat<2,2,double>",
"SparseMatrixclass ngbla::Mat<3,3,class std::complex<double> >",
"SparseMatrixclass ngbla::Mat<3,3,double>",
"SparseMatrixclass std::complex<double>",
"SparseMatrixdouble",
"SparseMatrixVariableBlocks",
"SparseMatrixSymmetricclass ngbla::Mat<2,2,class std::complex<double> >",
"SparseMatrixSymmetricclass ngbla::Mat<2,2,double>",
"SparseMatrixSymmetricclass ngbla::Mat<3,3,class std::complex<double> >",
"SparseMatrixSymmetricclass ngbla::Mat<3,3,double>",
"SparseMatrixSymmetricclass std::complex<double>",
"SparseMatrixSymmetricdouble",
"ArnoldiSolver",
"CGSolver",
"ChebyshevIteration",
"CreateParallelVector",
"CreateVVector",
"DoArchive",
"EigenValues_Preconditioner",
"GMRESSolver",
"InnerProduct",
"ParallelMatrix",
"QMRSolver",
"Sum"
]
class BaseMatrix():
    pass
class BaseSparseMatrix(BaseMatrix):
    """
    sparse matrix of any type
    """
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
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
    pass
class BlockMatrix(BaseMatrix):
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
    def __getitem__(self, inds: tuple) -> BaseMatrix: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    def __init__(self, mats: List[List[BaseMatrix]]) -> None: 
        """
        Make BlockMatrix with given array of matrices
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
    def col_nblocks(self) -> int:
        """
        :type: int
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
    def row_nblocks(self) -> int:
        """
        :type: int
        """
    @property
    def width(self) -> int:
        """
        Width of the matrix

        :type: int
        """
    pass
class BlockSmoother(BaseMatrix):
    """
    block Jacobi and block Gauss-Seidel smoothing
    """
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
    def Smooth(self, x: BaseVector, b: BaseVector, steps: int = 1) -> None: 
        """
        performs steps block-Gauss-Seidel iterations for the linear system A x = b
        """
    def SmoothBack(self, x: BaseVector, b: BaseVector, steps: int = 1) -> None: 
        """
        performs steps block-Gauss-Seidel iterations for the linear system A x = b in reverse order
        """
    def Update(self) -> None: 
        """
        Update matrix
        """
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    def __iadd__(self, mat: BaseMatrix) -> None: ...
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
class BlockVector(BaseVector):
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
    def __getitem__(self, ind: int) -> BaseVector: 
        """
        Return block at given position
        """
    def __getstate__(self) -> tuple: ...
    @overload
    def __iadd__(self, vec: BaseVector) -> BaseVector: ...
    @overload
    def __iadd__(self, arg0: ngla::DynamicVectorExpression) -> BaseVector: ...
    @overload
    def __imul__(self, value: complex) -> BaseVector: ...
    @overload
    def __imul__(self, value: float) -> BaseVector: ...
    def __init__(self, vecs: List[BaseVector]) -> None: 
        """
        Makes BlockVector by given array of vectors
        """
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
    def nblocks(self) -> int:
        """
        number of blocks in BlockVector

        :type: int
        """
    @property
    def size(self) -> int:
        """
        :type: int
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE171648>, '__dict__': <attribute '__dict__' of 'ngsolve.la.BlockVector' objects>, '__doc__': None, '__module__': 'ngsolve.la', '__getitem__': <instancemethod __getitem__ at 0x0000026DEE1716A8>, 'nblocks': <property object at 0x0000026DEE170A48>})
    pass
class ConstEBEMatrix(BaseMatrix):
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
    def __init__(self, h: int, w: int, matrix: ngsolve.bla.MatrixD, col_ind: list, row_ind: list) -> None: ...
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
class DynamicVectorExpression():
    def __add__(self, arg0: DynamicVectorExpression) -> DynamicVectorExpression: ...
    def __init__(self, arg0: BaseVector) -> None: ...
    def __neg__(self) -> DynamicVectorExpression: ...
    @overload
    def __rmul__(self, arg0: complex) -> DynamicVectorExpression: ...
    @overload
    def __rmul__(self, arg0: float) -> DynamicVectorExpression: ...
    def __sub__(self, arg0: DynamicVectorExpression) -> DynamicVectorExpression: ...
    pass
class Embedding(BaseMatrix):
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
    def __init__(self, height: int, range: ngsolve.ngstd.IntRange) -> None: 
        """
        Linear operator embedding a shorter vector into a longer vector
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
    def T(self) -> ngla::EmbeddingTranspose:
        """
        Return transpose of matrix

        :type: ngla::EmbeddingTranspose
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
class EmbeddingTranspose(BaseMatrix):
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
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    def __rmatmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
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
class IdentityMatrix():
    pass
class KrylovSpaceSolver(BaseMatrix):
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
    def GetSteps(self) -> int: ...
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
class MultiVectorExpr():
    def Scale(self, arg0: ngsolve.bla.VectorD) -> MultiVectorExpr: ...
    def __add__(self, arg0: MultiVectorExpr) -> MultiVectorExpr: ...
    def __neg__(self) -> MultiVectorExpr: ...
    def __sub__(self, arg0: MultiVectorExpr) -> MultiVectorExpr: ...
    pass
class MultiVector():
    pass
class PARALLEL_STATUS():
    """
    enum of possible parallel statuses

    Members:

      DISTRIBUTED

      CUMULATED

      NOT_PARALLEL
    """
    def __init__(self, arg0: int) -> None: ...
    def __int__(self) -> int: ...
    @property
    def name(self) -> str:
        """
        (self: handle) -> str

        :type: str
        """
    CUMULATED: ngsolve.la.PARALLEL_STATUS # value = PARALLEL_STATUS.CUMULATED
    DISTRIBUTED: ngsolve.la.PARALLEL_STATUS # value = PARALLEL_STATUS.DISTRIBUTED
    NOT_PARALLEL: ngsolve.la.PARALLEL_STATUS # value = PARALLEL_STATUS.NOT_PARALLEL
    __entries: dict # value = {'DISTRIBUTED': (PARALLEL_STATUS.DISTRIBUTED, None), 'CUMULATED': (PARALLEL_STATUS.CUMULATED, None), 'NOT_PARALLEL': (PARALLEL_STATUS.NOT_PARALLEL, None)}
    __members__: dict # value = {'DISTRIBUTED': PARALLEL_STATUS.DISTRIBUTED, 'CUMULATED': PARALLEL_STATUS.CUMULATED, 'NOT_PARALLEL': PARALLEL_STATUS.NOT_PARALLEL}
    pass
class ParallelDofs():
    def Dof2Proc(self, dof: int) -> pyngcore.FlatArray_I_S: ...
    def ExchangeProcs(self) -> pyngcore.FlatArray_I_S: ...
    def Proc2Dof(self, proc: int) -> pyngcore.FlatArray_I_S: ...
    @property
    def ndofglobal(self) -> int:
        """
        number of global degrees of freedom

        :type: int
        """
    @property
    def ndoflocal(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    pass
class PermutationMatrix(BaseMatrix):
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
    def __init__(self, w: int, ind: List[int]) -> None: ...
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
class Projector():
    pass
class QMRSolverC(BaseMatrix):
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
class QMRSolverD(BaseMatrix):
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
class Real2ComplexMatrix(BaseMatrix):
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
    def __init__(self, arg0: BaseMatrix) -> None: ...
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
class S_BaseMatrixC(BaseMatrix):
    """
    base sparse matrix
    """
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
class S_BaseMatrixD(BaseMatrix):
    """
    base sparse matrix
    """
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
class Smoother(BaseMatrix):
    """
    Jacobi and Gauss-Seidel smoothing
    """
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
    def Smooth(self, x: BaseVector, b: BaseVector) -> None: 
        """
        performs one step Gauss-Seidel iteration for the linear system A x = b
        """
    def SmoothBack(self, x: BaseVector, b: BaseVector) -> None: 
        """
        performs one step Gauss-Seidel iteration for the linear system A x = b in reverse order
        """
    def Update(self) -> None: 
        """
        Update matrix
        """
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    def __iadd__(self, mat: BaseMatrix) -> None: ...
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
class SparseFactorization(BaseMatrix):
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
    def Smooth(self, arg0: BaseVector, arg1: BaseVector) -> None: 
        """
        perform smoothing step (needs non-symmetric storage so symmetric sparse matrix)
        """
    def Update(self) -> None: 
        """
        Update matrix
        """
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    def __iadd__(self, mat: BaseMatrix) -> None: ...
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
class SparseCholesky_d(SparseFactorization, BaseMatrix):
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
    def Smooth(self, arg0: BaseVector, arg1: BaseVector) -> None: 
        """
        perform smoothing step (needs non-symmetric storage so symmetric sparse matrix)
        """
    def Update(self) -> None: 
        """
        Update matrix
        """
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    def __iadd__(self, mat: BaseMatrix) -> None: ...
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
class SparseCholesky_c(SparseFactorization, BaseMatrix):
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
    def Smooth(self, arg0: BaseVector, arg1: BaseVector) -> None: 
        """
        perform smoothing step (needs non-symmetric storage so symmetric sparse matrix)
        """
    def Update(self) -> None: 
        """
        Update matrix
        """
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    def __iadd__(self, mat: BaseMatrix) -> None: ...
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
class SparseMatrixDynamic(BaseMatrix):
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
    def __init__(self, arg0: BaseMatrix) -> None: ...
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
class SparseMatrixclass ngbla::Mat<2,2,class std::complex<double> >(BaseSparseMatrix, S_BaseMatrixC, BaseMatrix):
    """
    a sparse matrix in CSR storage
    """
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self,2,class std::complex<double> >) -> object: ...
    def CSR(self,2,class std::complex<double> >) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self,2,class std::complex<double> >, pos: tuple) -> ngsolve.bla.Mat2C: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self,2,class std::complex<double> >, pos: tuple, value: ngsolve.bla.Mat2C) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixclass ngbla::Mat<2,2,double>(BaseSparseMatrix, S_BaseMatrixD, BaseMatrix):
    """
    a sparse matrix in CSR storage
    """
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self,2,double>) -> object: ...
    def CSR(self,2,double>) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self,2,double>, pos: tuple) -> ngsolve.bla.Mat2D: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self,2,double>, pos: tuple, value: ngsolve.bla.Mat2D) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixclass ngbla::Mat<3,3,class std::complex<double> >(BaseSparseMatrix, S_BaseMatrixC, BaseMatrix):
    """
    a sparse matrix in CSR storage
    """
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self,3,class std::complex<double> >) -> object: ...
    def CSR(self,3,class std::complex<double> >) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self,3,class std::complex<double> >, pos: tuple) -> ngsolve.bla.Mat3C: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self,3,class std::complex<double> >, pos: tuple, value: ngsolve.bla.Mat3C) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixclass ngbla::Mat<3,3,double>(BaseSparseMatrix, S_BaseMatrixD, BaseMatrix):
    """
    a sparse matrix in CSR storage
    """
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self,3,double>) -> object: ...
    def CSR(self,3,double>) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self,3,double>, pos: tuple) -> ngsolve.bla.Mat3D: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self,3,double>, pos: tuple, value: ngsolve.bla.Mat3D) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixclass std::complex<double>(BaseSparseMatrix, S_BaseMatrixC, BaseMatrix):
    """
    a sparse matrix in CSR storage
    """
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self) -> object: ...
    def CSR(self) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self, pos: tuple) -> complex: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self, pos: tuple, value: complex) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixdouble(BaseSparseMatrix, S_BaseMatrixD, BaseMatrix):
    """
    a sparse matrix in CSR storage
    """
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self) -> object: ...
    def CSR(self) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self, pos: tuple) -> float: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self, pos: tuple, value: float) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixVariableBlocks(BaseMatrix):
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
    def __init__(self, arg0: BaseMatrix) -> None: ...
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
class SparseMatrixSymmetricclass ngbla::Mat<2,2,class std::complex<double> >(SparseMatrixclass ngbla::Mat<2,2,class std::complex<double> >, BaseSparseMatrix, S_BaseMatrixC, BaseMatrix):
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self,2,class std::complex<double> >) -> object: ...
    def CSR(self,2,class std::complex<double> >) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self,2,class std::complex<double> >, pos: tuple) -> ngsolve.bla.Mat2C: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self,2,class std::complex<double> >, pos: tuple, value: ngsolve.bla.Mat2C) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixSymmetricclass ngbla::Mat<2,2,double>(SparseMatrixclass ngbla::Mat<2,2,double>, BaseSparseMatrix, S_BaseMatrixD, BaseMatrix):
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self,2,double>) -> object: ...
    def CSR(self,2,double>) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self,2,double>, pos: tuple) -> ngsolve.bla.Mat2D: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self,2,double>, pos: tuple, value: ngsolve.bla.Mat2D) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixSymmetricclass ngbla::Mat<3,3,class std::complex<double> >(SparseMatrixclass ngbla::Mat<3,3,class std::complex<double> >, BaseSparseMatrix, S_BaseMatrixC, BaseMatrix):
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self,3,class std::complex<double> >) -> object: ...
    def CSR(self,3,class std::complex<double> >) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self,3,class std::complex<double> >, pos: tuple) -> ngsolve.bla.Mat3C: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self,3,class std::complex<double> >, pos: tuple, value: ngsolve.bla.Mat3C) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixSymmetricclass ngbla::Mat<3,3,double>(SparseMatrixclass ngbla::Mat<3,3,double>, BaseSparseMatrix, S_BaseMatrixD, BaseMatrix):
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self,3,double>) -> object: ...
    def CSR(self,3,double>) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self,3,double>, pos: tuple) -> ngsolve.bla.Mat3D: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self,3,double>, pos: tuple, value: ngsolve.bla.Mat3D) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixSymmetricclass std::complex<double>(SparseMatrixclass std::complex<double>, BaseSparseMatrix, S_BaseMatrixC, BaseMatrix):
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self) -> object: ...
    def CSR(self) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self, pos: tuple) -> complex: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self, pos: tuple, value: complex) -> None: 
        """
        Set value at given position
        """
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
class SparseMatrixSymmetricdouble(SparseMatrixdouble, BaseSparseMatrix, S_BaseMatrixD, BaseMatrix):
    def AsVector(self) -> BaseVector: 
        """
        Interprets the matrix values as a vector
        """
    def COO(self) -> object: ...
    def CSR(self) -> object: ...
    def CreateBlockSmoother(self, blocks: object, parallel: bool = False) -> ngla::BaseBlockJacobiPrecond: ...
    def CreateColVector(self) -> BaseVector: ...
    @staticmethod
    def CreateFromCOO(indi: list, indj: list, values: list, h: int, w: int) -> ngla::SparseMatrixTM<double>: ...
    @staticmethod
    def CreateFromElmat(col_ind: list, row_ind: list, matrices: list, h: int, w: int) -> SparseMatrixdouble: ...
    def CreateMatrix(self) -> BaseMatrix: 
        """
        Create matrix of same dimension and same sparsestructure
        """
    def CreateRowVector(self) -> BaseVector: ...
    def CreateSmoother(self, freedofs: pyngcore.BitArray = None) -> ngla::BaseJacobiPrecond: ...
    def CreateTranspose(self) -> ngla::SparseMatrixTM<double>: 
        """
        Return transposed matrix
        """
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
    def __getitem__(self, pos: tuple) -> float: 
        """
        Return value at given position
        """
    def __iadd__(self, mat: BaseMatrix) -> None: ...
    @overload
    def __matmul__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @overload
    def __matmul__(self, mat: SparseMatrixdouble) -> ngla::SparseMatrixTM<double>: ...
    @overload
    def __mul__(self, arg0: BaseVector) -> ngla::DynamicVectorExpression: ...
    @overload
    def __mul__(self, arg0: MultiVector) -> MultiVectorExpr: ...
    def __neg__(self) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: float) -> BaseMatrix: ...
    @overload
    def __rmul__(self, value: complex) -> BaseMatrix: ...
    def __setitem__(self, pos: tuple, value: float) -> None: 
        """
        Set value at given position
        """
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
def ArnoldiSolver(mata: BaseMatrix, matm: BaseMatrix, freedofs: pyngcore.BitArray, vecs: list, shift: complex = <ngsolve.ngstd.DummyArgument>) -> ngsolve.bla.VectorC:
    """
    Shift-and-invert Arnoldi eigenvalue solver

    Solves the generalized linear EVP A*u = M*lam*u using an Arnoldi iteration for the 
    shifted EVP (A-shift*M)^(-1)*M*u = lam*u with a Krylow space of dimension 2*len(vecs)+1.
    len(vecs) eigenpairs with the closest eigenvalues to the shift are returned.

    Parameters:

    mata : ngsolve.la.BaseMatrix
      matrix A

    matm : ngsolve.la.BaseMatrix
      matrix M

    freedofs : nsolve.ngstd.BitArray
      correct degrees of freedom

    vecs : list
      list of BaseVectors for writing eigenvectors

    shift : object
      complex or real shift
    """
def CGSolver(mat: BaseMatrix, pre: BaseMatrix, complex: bool = False, printrates: bool = True, precision: float = 1e-08, maxsteps: int = 200) -> KrylovSpaceSolver:
    """
    A CG Solver.

    Parameters:

    mat : ngsolve.la.BaseMatrix
      input matrix 

    pre : ngsolve.la.BaseMatrix
      input preconditioner matrix

    complex : bool
      input complex, if not set it is deduced from matrix type

    printrates : bool
      input printrates

    precision : float
      input requested precision. CGSolver stops if precision is reached.

    maxsteps : int
      input maximal steps. CGSolver stops after this steps.
    """
def ChebyshevIteration(mat: BaseMatrix = None, pre: BaseMatrix = None, steps: int = 3, lam_min: float = 1, lam_max: float = 1) -> BaseMatrix:
    pass
def CreateParallelVector(pardofs: ParallelDofs) -> ngla::BaseVector:
    pass
def CreateVVector(size: int, complex: bool = False, entrysize: int = 1) -> ngla::BaseVector:
    pass
def DoArchive(arg0: ngsolve.ngstd.Archive, arg1: BaseMatrix) -> ngsolve.ngstd.Archive:
    pass
def EigenValues_Preconditioner(mat: BaseMatrix, pre: BaseMatrix, tol: float = 1e-10) -> ngsolve.bla.VectorD:
    """
    Calculate eigenvalues of pre * mat, where pre and mat are positive definite matrices.
    The typical usecase of this function is to calculate the condition number of a preconditioner.It uses the Lanczos algorithm and bisection for the tridiagonal matrix
    """
def GMRESSolver(mat: BaseMatrix, pre: BaseMatrix, printrates: bool = True, precision: float = 1e-08, maxsteps: int = 200) -> KrylovSpaceSolver:
    """
    A General Minimal Residuum (GMRES) Solver.

    Parameters:

    mat : ngsolve.la.BaseMatrix
      input matrix 

    pre : ngsolve.la.BaseMatrix
      input preconditioner matrix

    printrates : bool
      input printrates

    precision : float
      input requested precision. GMRESSolver stops if precision is reached.

    maxsteps : int
      input maximal steps. GMRESSolver stops after this steps.
    """
def InnerProduct(x: object, y: object, **kwargs) -> object:
    """
    Computes InnerProduct of given objects
    """
@overload
def ParallelMatrix(mat: object = None, row_pardofs: object = None, col_pardofs: object = None, op: object = None) -> None:
    pass
@overload
def ParallelMatrix(mat: object = None, pardofs: object = None, op: object = None) -> None:
    pass
def QMRSolver(mat: BaseMatrix, pre: BaseMatrix, printrates: bool = True, precision: float = 1e-08, maxsteps: int = 200) -> KrylovSpaceSolver:
    """
    A Quasi Minimal Residuum (QMR) Solver.

    Parameters:

    mat : ngsolve.la.BaseMatrix
      input matrix 

    pre : ngsolve.la.BaseMatrix
      input preconditioner matrix

    printrates : bool
      input printrates

    precision : float
      input requested precision. QMRSolver stops if precision is reached.

    maxsteps : int
      input maximal steps. QMRSolver stops after this steps.
    """
def Sum(arg0: DynamicVectorExpression, arg1: DynamicVectorExpression) -> DynamicVectorExpression:
    pass
