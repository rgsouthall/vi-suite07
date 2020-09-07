"""module for perfectly matched layers"""
import ngsolve.comp.pml
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import ngsolve.fem
import ngsolve.bla
import ngsolve.ngstd
__all__  = [
"PML",
"BrickRadial",
"Cartesian",
"Compound",
"Custom",
"HalfSpace",
"Radial"
]
class PML():
    """
    Base PML object

    can only be created by generator functions. Use PML(x, [y, z]) to evaluate the scaling.
    """
    def __add__(self, pml: PML) -> PML: ...
    @staticmethod
    def __call__(*args) -> ngsolve.bla.VectorC: 
        """
        map a point
        """
    def __str__(self) -> str: ...
    @staticmethod
    def call_jacobian(*args) -> ngsolve.bla.MatrixC: 
        """
        evaluate PML jacobian at point x, [y, z]
        """
    @property
    def Det_CF(self) -> ngsolve.fem.CoefficientFunction:
        """
        the determinant of the jacobian as coefficient function

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def JacInv_CF(self) -> ngsolve.fem.CoefficientFunction:
        """
        the inverse of the jacobian as coefficient function

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def Jac_CF(self) -> ngsolve.fem.CoefficientFunction:
        """
        the jacobian of the PML as coefficient function

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def PML_CF(self) -> ngsolve.fem.CoefficientFunction:
        """
        the scaling as coefficient function

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def dim(self) -> int:
        """
        dimension

        :type: int
        """
    pass
def BrickRadial(mins: object, maxs: object, origin: object = (0.0, 0.0, 0.0), alpha: complex = 1j) -> PML:
    """
    BrickRadial(mins: object, maxs: object, origin: object = (0.0, 0.0, 0.0), alpha: complex = 1j) -> ngsolve.comp.pml.PML

    radial pml on a brick

    mins, maxs and origin are given as tuples/lists
    """
def Cartesian(mins: object, maxs: object, alpha: complex = 1j) -> PML:
    """
    cartesian pml transformation

    mins and maxs are tuples/lists determining the dimension
    """
def Compound(pml1: PML, pml2: PML, dims1: object = <ngsolve.ngstd.DummyArgument>, dims2: object = <ngsolve.ngstd.DummyArgument>) -> PML:
    """
    tensor product of two pml transformations

            dimensions are optional, given as tuples/lists and start with 1
    """
def Custom(trafo: ngsolve.fem.CoefficientFunction, jac: ngsolve.fem.CoefficientFunction) -> PML:
    """
    custom pml transformation

    trafo and jac are coefficient functions of the scaling and the jacobian
    """
def HalfSpace(point: object, normal: object, alpha: complex = 1j) -> PML:
    """
    half space pml

    scales orthogonal to specified plane in direction of normal point and normal are given as tuples/lists determining the dimension
    """
def Radial(origin: object, rad: float = 1, alpha: complex = 1j) -> PML:
    """
    radial pml transformation

    origin is a list/tuple with as many entries as dimenson
    """
