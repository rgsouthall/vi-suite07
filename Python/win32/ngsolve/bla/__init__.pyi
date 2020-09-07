"""pybind bla"""
import ngsolve.bla
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
__all__  = [
"FlatMatrixC",
"FlatMatrixD",
"FlatVectorC",
"FlatVectorD",
"Mat2C",
"Mat2D",
"Mat3C",
"Mat3D",
"MatrixC",
"MatrixD",
"SliceVectorC",
"SliceVectorD",
"Vec1D",
"Vec2D",
"Vec3D",
"VectorC",
"VectorD",
"CheckPerformance",
"InnerProduct",
"Matrix",
"Norm",
"Vector",
"__timing__"
]
class FlatMatrixC():
    def Height(self) -> int: 
        """
        Return height of matrix

        Returns height of the matrix
        """
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Width(self) -> int: 
        """
        Return width of matrix

        Returns width of the matrix
        """
    @overload
    def __add__(self, mat: FlatMatrixC) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __add__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __getitem__(self, arg0: tuple) -> object: ...
    @overload
    def __getitem__(self, arg0: int) -> VectorC: ...
    @overload
    def __getitem__(self, arg0: slice) -> ngbla::Matrix<std::complex<double>,1>: ...
    def __iadd__(self, arg0: FlatMatrixC) -> FlatMatrixC: ...
    def __imul__(self, arg0: complex) -> FlatMatrixC: ...
    def __isub__(self, arg0: FlatMatrixC) -> FlatMatrixC: ...
    def __len__(self) -> int: 
        """
        Return height of matrix
        """
    @overload
    def __mul__(self, vec: FlatVectorC) -> VectorC: ...
    @overload
    def __mul__(self, values: complex) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __mul__(self, mat: FlatMatrixC) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __mul__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __mul__(self, value: float) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __mul__(self, vec: FlatVectorD) -> VectorC: ...
    def __neg__(self) -> ngbla::Matrix<std::complex<double>,1>: ...
    def __radd__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    def __repr__(self) -> str: ...
    @overload
    def __rmul__(self, value: float) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __rmul__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __rmul__(self, value: complex) -> ngbla::Matrix<std::complex<double>,1>: ...
    def __rsub__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: complex) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: FlatVectorC) -> None: ...
    @overload
    def __setitem__(self, arg0: int, arg1: VectorC) -> None: ...
    @overload
    def __setitem__(self, arg0: int, arg1: complex) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: complex) -> None: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: FlatMatrixC) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: FlatMatrixC) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, mat: FlatMatrixC) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __sub__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @property
    def A(self) -> VectorC:
        """
        Returns matrix as vector

        :type: VectorC
        """
    @property
    def C(self) -> ngbla::Matrix<std::complex<double>,1>:
        """
        Return conjugate matrix

        :type: ngbla::Matrix<std::complex<double>,1>
        """
    @property
    def H(self) -> ngbla::Matrix<std::complex<double>,1>:
        """
        Return conjugate and transposed matrix

        :type: ngbla::Matrix<std::complex<double>,1>
        """
    @property
    def I(self) -> ngbla::Matrix<std::complex<double>,1>:
        """
        :type: ngbla::Matrix<std::complex<double>,1>
        """
    @property
    def T(self) -> ngbla::Matrix<std::complex<double>,1>:
        """
        Return transpose of matrix

        :type: ngbla::Matrix<std::complex<double>,1>
        """
    @property
    def diag(self) -> VectorC:
        """
        :type: VectorC
        """
    @diag.setter
    def diag(self, arg1: FlatVectorC) -> None:
        pass
    @property
    def h(self) -> int:
        """
        Height of the matrix

        :type: int
        """
    @property
    def w(self) -> int:
        """
        Width of the matrix

        :type: int
        """
    pass
class FlatMatrixD():
    def Height(self) -> int: 
        """
        Return height of matrix
        """
    def Inverse(self, arg0: FlatMatrixD) -> None: ...
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Width(self) -> int: 
        """
        Return width of matrix
        """
    def __add__(self, mat: FlatMatrixD) -> ngbla::Matrix<double,1>: ...
    @overload
    def __getitem__(self, arg0: tuple) -> object: ...
    @overload
    def __getitem__(self, arg0: int) -> VectorD: ...
    @overload
    def __getitem__(self, arg0: slice) -> ngbla::Matrix<double,1>: ...
    def __iadd__(self, arg0: FlatMatrixD) -> FlatMatrixD: ...
    def __imul__(self, arg0: float) -> FlatMatrixD: ...
    def __isub__(self, arg0: FlatMatrixD) -> FlatMatrixD: ...
    def __len__(self) -> int: 
        """
        Return height of matrix
        """
    @overload
    def __mul__(self, mat: FlatMatrixD) -> ngbla::Matrix<double,1>: ...
    @overload
    def __mul__(self, values: float) -> ngbla::Matrix<double,1>: ...
    @overload
    def __mul__(self, vec: FlatVectorD) -> VectorD: ...
    def __neg__(self) -> ngbla::Matrix<double,1>: ...
    def __repr__(self) -> str: ...
    def __rmul__(self, value: float) -> ngbla::Matrix<double,1>: ...
    @overload
    def __setitem__(self, arg0: int, arg1: VectorD) -> None: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: FlatMatrixD) -> None: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: float) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: float) -> None: ...
    @overload
    def __setitem__(self, arg0: int, arg1: float) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: FlatVectorD) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: FlatMatrixD) -> None: ...
    def __str__(self) -> str: ...
    def __sub__(self, mat: FlatMatrixD) -> ngbla::Matrix<double,1>: ...
    @property
    def A(self) -> VectorD:
        """
        Returns matrix as vector

        :type: VectorD
        """
    @property
    def I(self) -> ngbla::Matrix<double,1>:
        """
        :type: ngbla::Matrix<double,1>
        """
    @property
    def T(self) -> ngbla::Matrix<double,1>:
        """
        return transpose of matrix

        :type: ngbla::Matrix<double,1>
        """
    @property
    def diag(self) -> VectorD:
        """
        :type: VectorD
        """
    @diag.setter
    def diag(self, arg1: FlatVectorD) -> None:
        pass
    @property
    def h(self) -> int:
        """
        Height of the matrix

        :type: int
        """
    @property
    def w(self) -> int:
        """
        Width of the matrix

        :type: int
        """
    pass
class FlatVectorC():
    def Get(self, pos: int) -> complex: 
        """
        Return value at given position
        """
    def InnerProduct(self, y: FlatVectorC) -> complex: 
        """
        Returns InnerProduct with other object
        """
    def Norm(self) -> float: 
        """
        Returns L2-norm
        """
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Range(self, arg0: int, arg1: int) -> FlatVectorC: ...
    def Set(self, pos: int, value: complex) -> None: 
        """
        Set value at given position
        """
    def __add__(self, vec: FlatVectorC) -> ngbla::Vector<std::complex<double> >: ...
    @overload
    def __getitem__(self, ind: list) -> ngbla::Vector<std::complex<double> >: 
        """
        Return value at given position

        Return values at given positions

        Return values at given positions
        """
    @overload
    def __getitem__(self, pos: int) -> complex: ...
    @overload
    def __getitem__(self, inds: slice) -> ngbla::Vector<std::complex<double> >: ...
    def __iadd__(self, arg0: FlatVectorC) -> FlatVectorC: ...
    @overload
    def __imul__(self, arg0: complex) -> FlatVectorC: ...
    @overload
    def __imul__(self, arg0: float) -> FlatVectorC: ...
    def __init__(self, arg0: int, arg1: complex) -> None: ...
    def __isub__(self, arg0: FlatVectorC) -> FlatVectorC: ...
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: 
        """
        Return length of the array
        """
    def __mul__(self, value: complex) -> ngbla::Vector<std::complex<double> >: ...
    def __neg__(self) -> ngbla::Vector<std::complex<double> >: ...
    def __repr__(self) -> str: ...
    def __rmul__(self, value: complex) -> ngbla::Vector<std::complex<double> >: ...
    @overload
    def __setitem__(self, pos: int, value: complex) -> None: 
        """
        Set value at given position

        Set values at given positions

        Set value at given positions
        """
    @overload
    def __setitem__(self, inds: slice, rv: FlatVectorC) -> None: ...
    @overload
    def __setitem__(self, inds: slice, value: complex) -> None: ...
    def __str__(self) -> str: ...
    def __sub__(self, vec: FlatVectorC) -> ngbla::Vector<std::complex<double> >: ...
    pass
class FlatVectorD():
    def Get(self, pos: int) -> float: 
        """
        Return value at given position
        """
    def InnerProduct(self, y: FlatVectorD) -> float: 
        """
        Returns InnerProduct with other object
        """
    def Norm(self) -> float: 
        """
        Returns L2-norm
        """
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Range(self, arg0: int, arg1: int) -> FlatVectorD: ...
    def Set(self, pos: int, value: float) -> None: 
        """
        Set value at given position
        """
    def __add__(self, vec: FlatVectorD) -> ngbla::Vector<double>: ...
    @overload
    def __getitem__(self, inds: slice) -> ngbla::Vector<double>: 
        """
        Return value at given position

        Return values at given positions

        Return values at given positions
        """
    @overload
    def __getitem__(self, ind: list) -> ngbla::Vector<double>: ...
    @overload
    def __getitem__(self, pos: int) -> float: ...
    def __iadd__(self, arg0: FlatVectorD) -> FlatVectorD: ...
    def __imul__(self, arg0: float) -> FlatVectorD: ...
    def __init__(self, arg0: int, arg1: float) -> None: ...
    def __isub__(self, arg0: FlatVectorD) -> FlatVectorD: ...
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: 
        """
        Return length of the array
        """
    def __mul__(self, value: float) -> ngbla::Vector<double>: ...
    def __neg__(self) -> ngbla::Vector<double>: ...
    def __repr__(self) -> str: ...
    def __rmul__(self, value: float) -> ngbla::Vector<double>: ...
    @overload
    def __setitem__(self, inds: slice, rv: FlatVectorD) -> None: 
        """
        Set value at given position

        Set values at given positions

        Set value at given positions
        """
    @overload
    def __setitem__(self, pos: int, value: float) -> None: ...
    @overload
    def __setitem__(self, inds: slice, value: float) -> None: ...
    def __str__(self) -> str: ...
    def __sub__(self, vec: FlatVectorD) -> ngbla::Vector<double>: ...
    pass
class Mat2C():
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def __getitem__(self, arg0: tuple) -> complex: ...
    pass
class Mat2D():
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def __getitem__(self, arg0: tuple) -> float: ...
    pass
class Mat3C():
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def __getitem__(self, arg0: tuple) -> complex: ...
    pass
class Mat3D():
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def __getitem__(self, arg0: tuple) -> float: ...
    pass
class MatrixC(FlatMatrixC):
    def Height(self) -> int: 
        """
        Return height of matrix

        Returns height of the matrix
        """
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Width(self) -> int: 
        """
        Return width of matrix

        Returns width of the matrix
        """
    @overload
    def __add__(self, mat: FlatMatrixC) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __add__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __getitem__(self, arg0: tuple) -> object: ...
    @overload
    def __getitem__(self, arg0: int) -> VectorC: ...
    @overload
    def __getitem__(self, arg0: slice) -> ngbla::Matrix<std::complex<double>,1>: ...
    def __iadd__(self, arg0: MatrixC) -> MatrixC: ...
    def __imul__(self, arg0: complex) -> MatrixC: ...
    def __init__(self, n: int, m: int) -> None: 
        """
        Makes matrix of dimension n x m
        """
    def __isub__(self, arg0: MatrixC) -> MatrixC: ...
    def __len__(self) -> int: 
        """
        Return height of matrix
        """
    @overload
    def __mul__(self, vec: FlatVectorC) -> VectorC: ...
    @overload
    def __mul__(self, values: complex) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __mul__(self, mat: FlatMatrixC) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __mul__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __mul__(self, value: float) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __mul__(self, vec: FlatVectorD) -> VectorC: ...
    def __neg__(self) -> ngbla::Matrix<std::complex<double>,1>: ...
    def __radd__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    def __repr__(self) -> str: ...
    @overload
    def __rmul__(self, value: float) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __rmul__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __rmul__(self, value: complex) -> ngbla::Matrix<std::complex<double>,1>: ...
    def __rsub__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: complex) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: FlatVectorC) -> None: ...
    @overload
    def __setitem__(self, arg0: int, arg1: VectorC) -> None: ...
    @overload
    def __setitem__(self, arg0: int, arg1: complex) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: complex) -> None: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: FlatMatrixC) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: FlatMatrixC) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, mat: FlatMatrixC) -> ngbla::Matrix<std::complex<double>,1>: ...
    @overload
    def __sub__(self, mat: FlatMatrixD) -> ngbla::Matrix<std::complex<double>,1>: ...
    @property
    def A(self) -> VectorC:
        """
        Returns matrix as vector

        :type: VectorC
        """
    @property
    def C(self) -> ngbla::Matrix<std::complex<double>,1>:
        """
        Return conjugate matrix

        :type: ngbla::Matrix<std::complex<double>,1>
        """
    @property
    def H(self) -> ngbla::Matrix<std::complex<double>,1>:
        """
        Return conjugate and transposed matrix

        :type: ngbla::Matrix<std::complex<double>,1>
        """
    @property
    def I(self) -> ngbla::Matrix<std::complex<double>,1>:
        """
        :type: ngbla::Matrix<std::complex<double>,1>
        """
    @property
    def T(self) -> ngbla::Matrix<std::complex<double>,1>:
        """
        Return transpose of matrix

        :type: ngbla::Matrix<std::complex<double>,1>
        """
    @property
    def diag(self) -> VectorC:
        """
        :type: VectorC
        """
    @diag.setter
    def diag(self, arg1: FlatVectorC) -> None:
        pass
    @property
    def h(self) -> int:
        """
        Height of the matrix

        :type: int
        """
    @property
    def w(self) -> int:
        """
        Width of the matrix

        :type: int
        """
    pass
class MatrixD(FlatMatrixD):
    def Height(self) -> int: 
        """
        Return height of matrix
        """
    def Inverse(self, arg0: FlatMatrixD) -> None: ...
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Width(self) -> int: 
        """
        Return width of matrix
        """
    def __add__(self, mat: FlatMatrixD) -> ngbla::Matrix<double,1>: ...
    @overload
    def __getitem__(self, arg0: tuple) -> object: ...
    @overload
    def __getitem__(self, arg0: int) -> VectorD: ...
    @overload
    def __getitem__(self, arg0: slice) -> ngbla::Matrix<double,1>: ...
    def __iadd__(self, arg0: MatrixD) -> MatrixD: ...
    def __imul__(self, arg0: float) -> MatrixD: ...
    def __init__(self, n: int, m: int) -> None: 
        """
        Makes matrix of dimension n x m
        """
    def __isub__(self, arg0: MatrixD) -> MatrixD: ...
    def __len__(self) -> int: 
        """
        Return height of matrix
        """
    @overload
    def __mul__(self, mat: FlatMatrixD) -> ngbla::Matrix<double,1>: ...
    @overload
    def __mul__(self, values: float) -> ngbla::Matrix<double,1>: ...
    @overload
    def __mul__(self, vec: FlatVectorD) -> VectorD: ...
    def __neg__(self) -> ngbla::Matrix<double,1>: ...
    def __repr__(self) -> str: ...
    def __rmul__(self, value: float) -> ngbla::Matrix<double,1>: ...
    @overload
    def __setitem__(self, arg0: int, arg1: VectorD) -> None: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: FlatMatrixD) -> None: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: float) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: float) -> None: ...
    @overload
    def __setitem__(self, arg0: int, arg1: float) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: FlatVectorD) -> None: ...
    @overload
    def __setitem__(self, arg0: tuple, arg1: FlatMatrixD) -> None: ...
    def __str__(self) -> str: ...
    def __sub__(self, mat: FlatMatrixD) -> ngbla::Matrix<double,1>: ...
    @property
    def A(self) -> VectorD:
        """
        Returns matrix as vector

        :type: VectorD
        """
    @property
    def I(self) -> ngbla::Matrix<double,1>:
        """
        :type: ngbla::Matrix<double,1>
        """
    @property
    def T(self) -> ngbla::Matrix<double,1>:
        """
        return transpose of matrix

        :type: ngbla::Matrix<double,1>
        """
    @property
    def diag(self) -> VectorD:
        """
        :type: VectorD
        """
    @diag.setter
    def diag(self, arg1: FlatVectorD) -> None:
        pass
    @property
    def h(self) -> int:
        """
        Height of the matrix

        :type: int
        """
    @property
    def w(self) -> int:
        """
        Width of the matrix

        :type: int
        """
    pass
class SliceVectorC():
    def Get(self, pos: int) -> complex: 
        """
        Return value at given position
        """
    def InnerProduct(self, y: SliceVectorC) -> complex: 
        """
        Returns InnerProduct with other object
        """
    def Norm(self) -> float: 
        """
        Returns L2-norm
        """
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Range(self, arg0: int, arg1: int) -> SliceVectorC: ...
    def Set(self, pos: int, value: complex) -> None: 
        """
        Set value at given position
        """
    def __add__(self, vec: SliceVectorC) -> ngbla::Vector<std::complex<double> >: ...
    @overload
    def __getitem__(self, ind: list) -> ngbla::Vector<std::complex<double> >: 
        """
        Return value at given position

        Return values at given positions

        Return values at given positions
        """
    @overload
    def __getitem__(self, inds: slice) -> ngbla::Vector<std::complex<double> >: ...
    @overload
    def __getitem__(self, pos: int) -> complex: ...
    def __iadd__(self, arg0: SliceVectorC) -> SliceVectorC: ...
    @overload
    def __imul__(self, arg0: complex) -> SliceVectorC: ...
    @overload
    def __imul__(self, arg0: float) -> SliceVectorC: ...
    def __isub__(self, arg0: SliceVectorC) -> SliceVectorC: ...
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: 
        """
        Return length of the array
        """
    def __mul__(self, value: complex) -> ngbla::Vector<std::complex<double> >: ...
    def __neg__(self) -> ngbla::Vector<std::complex<double> >: ...
    def __repr__(self) -> str: ...
    def __rmul__(self, value: complex) -> ngbla::Vector<std::complex<double> >: ...
    @overload
    def __setitem__(self, inds: slice, rv: SliceVectorC) -> None: 
        """
        Set value at given position

        Set values at given positions

        Set value at given positions
        """
    @overload
    def __setitem__(self, pos: int, value: complex) -> None: ...
    @overload
    def __setitem__(self, inds: slice, value: complex) -> None: ...
    def __str__(self) -> str: ...
    def __sub__(self, vec: SliceVectorC) -> ngbla::Vector<std::complex<double> >: ...
    pass
class SliceVectorD():
    def Get(self, pos: int) -> float: 
        """
        Return value at given position
        """
    def InnerProduct(self, y: SliceVectorD) -> float: 
        """
        Returns InnerProduct with other object
        """
    def Norm(self) -> float: 
        """
        Returns L2-norm
        """
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Range(self, arg0: int, arg1: int) -> SliceVectorD: ...
    def Set(self, pos: int, value: float) -> None: 
        """
        Set value at given position
        """
    def __add__(self, vec: SliceVectorD) -> ngbla::Vector<double>: ...
    @overload
    def __getitem__(self, ind: list) -> ngbla::Vector<double>: 
        """
        Return value at given position

        Return values at given positions

        Return values at given positions
        """
    @overload
    def __getitem__(self, pos: int) -> float: ...
    @overload
    def __getitem__(self, inds: slice) -> ngbla::Vector<double>: ...
    def __iadd__(self, arg0: SliceVectorD) -> SliceVectorD: ...
    def __imul__(self, arg0: float) -> SliceVectorD: ...
    def __isub__(self, arg0: SliceVectorD) -> SliceVectorD: ...
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: 
        """
        Return length of the array
        """
    def __mul__(self, value: float) -> ngbla::Vector<double>: ...
    def __neg__(self) -> ngbla::Vector<double>: ...
    def __repr__(self) -> str: ...
    def __rmul__(self, value: float) -> ngbla::Vector<double>: ...
    @overload
    def __setitem__(self, inds: slice, rv: SliceVectorD) -> None: 
        """
        Set value at given position

        Set values at given positions

        Set value at given positions
        """
    @overload
    def __setitem__(self, pos: int, value: float) -> None: ...
    @overload
    def __setitem__(self, inds: slice, value: float) -> None: ...
    def __str__(self) -> str: ...
    def __sub__(self, vec: SliceVectorD) -> ngbla::Vector<double>: ...
    pass
class Vec1D():
    def Get(self, pos: int) -> float: 
        """
        Return value at given position
        """
    def InnerProduct(self, y: Vec1D) -> float: 
        """
        Returns InnerProduct with other object
        """
    def Norm(self) -> float: 
        """
        Returns L2-norm
        """
    def __add__(self, vec: Vec1D) -> Vec1D: ...
    @overload
    def __getitem__(self, pos: int) -> float: 
        """
        Return values at given positions

        Return values at given positions

        Return value at given position
        """
    @overload
    def __getitem__(self, inds: slice) -> Vec1D: ...
    @overload
    def __getitem__(self, ind: list) -> Vec1D: ...
    def __mul__(self, value: float) -> Vec1D: ...
    def __neg__(self) -> Vec1D: ...
    def __rmul__(self, value: float) -> Vec1D: ...
    @overload
    def __setitem__(self, inds: slice, value: float) -> None: 
        """
        Set values at given positions

        Set value at given positions
        """
    @overload
    def __setitem__(self, inds: slice, rv: Vec1D) -> None: ...
    def __sub__(self, vec: Vec1D) -> Vec1D: ...
    pass
class Vec2D():
    def Get(self, pos: int) -> float: 
        """
        Return value at given position
        """
    def InnerProduct(self, y: Vec2D) -> float: 
        """
        Returns InnerProduct with other object
        """
    def Norm(self) -> float: 
        """
        Returns L2-norm
        """
    def __add__(self, vec: Vec2D) -> Vec2D: ...
    @overload
    def __getitem__(self, ind: list) -> Vec2D: 
        """
        Return values at given positions

        Return values at given positions

        Return value at given position
        """
    @overload
    def __getitem__(self, inds: slice) -> Vec2D: ...
    @overload
    def __getitem__(self, pos: int) -> float: ...
    def __mul__(self, value: float) -> Vec2D: ...
    def __neg__(self) -> Vec2D: ...
    def __rmul__(self, value: float) -> Vec2D: ...
    @overload
    def __setitem__(self, inds: slice, value: float) -> None: 
        """
        Set values at given positions

        Set value at given positions
        """
    @overload
    def __setitem__(self, inds: slice, rv: Vec2D) -> None: ...
    def __sub__(self, vec: Vec2D) -> Vec2D: ...
    pass
class Vec3D():
    def Get(self, pos: int) -> float: 
        """
        Return value at given position
        """
    def InnerProduct(self, y: Vec3D) -> float: 
        """
        Returns InnerProduct with other object
        """
    def Norm(self) -> float: 
        """
        Returns L2-norm
        """
    def __add__(self, vec: Vec3D) -> Vec3D: ...
    @overload
    def __getitem__(self, inds: slice) -> Vec3D: 
        """
        Return values at given positions

        Return values at given positions

        Return value at given position
        """
    @overload
    def __getitem__(self, pos: int) -> float: ...
    @overload
    def __getitem__(self, ind: list) -> Vec3D: ...
    def __mul__(self, value: float) -> Vec3D: ...
    def __neg__(self) -> Vec3D: ...
    def __rmul__(self, value: float) -> Vec3D: ...
    @overload
    def __setitem__(self, inds: slice, value: float) -> None: 
        """
        Set values at given positions

        Set value at given positions
        """
    @overload
    def __setitem__(self, inds: slice, rv: Vec3D) -> None: ...
    def __sub__(self, vec: Vec3D) -> Vec3D: ...
    pass
class VectorC(FlatVectorC):
    def Get(self, pos: int) -> complex: 
        """
        Return value at given position
        """
    def InnerProduct(self, y: FlatVectorC) -> complex: 
        """
        Returns InnerProduct with other object
        """
    def Norm(self) -> float: 
        """
        Returns L2-norm
        """
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Range(self, arg0: int, arg1: int) -> FlatVectorC: ...
    def Set(self, pos: int, value: complex) -> None: 
        """
        Set value at given position
        """
    def __add__(self, vec: FlatVectorC) -> ngbla::Vector<std::complex<double> >: ...
    @overload
    def __getitem__(self, ind: list) -> ngbla::Vector<std::complex<double> >: 
        """
        Return value at given position

        Return values at given positions

        Return values at given positions
        """
    @overload
    def __getitem__(self, pos: int) -> complex: ...
    @overload
    def __getitem__(self, inds: slice) -> ngbla::Vector<std::complex<double> >: ...
    def __iadd__(self, arg0: VectorC) -> VectorC: ...
    def __imul__(self, arg0: complex) -> VectorC: ...
    def __init__(self, arg0: int) -> None: ...
    def __isub__(self, arg0: VectorC) -> VectorC: ...
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: 
        """
        Return length of the array
        """
    def __mul__(self, value: complex) -> ngbla::Vector<std::complex<double> >: ...
    def __neg__(self) -> ngbla::Vector<std::complex<double> >: ...
    def __repr__(self) -> str: ...
    def __rmul__(self, value: complex) -> ngbla::Vector<std::complex<double> >: ...
    @overload
    def __setitem__(self, pos: int, value: complex) -> None: 
        """
        Set value at given position

        Set values at given positions

        Set value at given positions
        """
    @overload
    def __setitem__(self, inds: slice, rv: FlatVectorC) -> None: ...
    @overload
    def __setitem__(self, inds: slice, value: complex) -> None: ...
    def __str__(self) -> str: ...
    def __sub__(self, vec: FlatVectorC) -> ngbla::Vector<std::complex<double> >: ...
    pass
class VectorD(FlatVectorD):
    def Get(self, pos: int) -> float: 
        """
        Return value at given position
        """
    def InnerProduct(self, y: FlatVectorD) -> float: 
        """
        Returns InnerProduct with other object
        """
    def Norm(self) -> float: 
        """
        Returns L2-norm
        """
    def NumPy(self) -> object: 
        """
        Return NumPy object
        """
    def Range(self, arg0: int, arg1: int) -> FlatVectorD: ...
    def Set(self, pos: int, value: float) -> None: 
        """
        Set value at given position
        """
    def __add__(self, vec: FlatVectorD) -> ngbla::Vector<double>: ...
    @overload
    def __getitem__(self, inds: slice) -> ngbla::Vector<double>: 
        """
        Return value at given position

        Return values at given positions

        Return values at given positions
        """
    @overload
    def __getitem__(self, ind: list) -> ngbla::Vector<double>: ...
    @overload
    def __getitem__(self, pos: int) -> float: ...
    def __iadd__(self, arg0: VectorD) -> VectorD: ...
    def __imul__(self, arg0: float) -> VectorD: ...
    def __init__(self, arg0: int) -> None: ...
    def __isub__(self, arg0: VectorD) -> VectorD: ...
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: 
        """
        Return length of the array
        """
    def __mul__(self, value: float) -> ngbla::Vector<double>: ...
    def __neg__(self) -> ngbla::Vector<double>: ...
    def __repr__(self) -> str: ...
    def __rmul__(self, value: float) -> ngbla::Vector<double>: ...
    @overload
    def __setitem__(self, inds: slice, rv: FlatVectorD) -> None: 
        """
        Set value at given position

        Set values at given positions

        Set value at given positions
        """
    @overload
    def __setitem__(self, pos: int, value: float) -> None: ...
    @overload
    def __setitem__(self, inds: slice, value: float) -> None: ...
    def __str__(self) -> str: ...
    def __sub__(self, vec: FlatVectorD) -> ngbla::Vector<double>: ...
    pass
def CheckPerformance(n: int, m: int, k: int) -> None:
    pass
def InnerProduct(x: object, y: object, **kwargs) -> object:
    """
    Compute InnerProduct
    """
@overload
def Matrix(arg0: List[List[complex]]) -> MatrixC:
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
def Matrix(height: int, width: int, complex: bool = False) -> object:
    pass
@overload
def Matrix(arg0: List[List[float]]) -> MatrixD:
    pass
def Norm(x: object) -> object:
    """
    Compute Norm
    """
@overload
def Vector(arg0: List[float]) -> VectorD:
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
def Vector(arg0: List[complex]) -> VectorC:
    pass
def __timing__(what: int, n: int, m: int, k: int, lapack: bool = False) -> List[Tuple[str, float]]:
    pass
