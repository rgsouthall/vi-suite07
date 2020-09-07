"""
NGSolve
=======

A high order finite element library

Modules:
ngsolve.bla .... simple vectors and matrices
ngsolve.fem .... finite elements and integrators
ngsolve.comp ... function spaces, forms
"""
import ngsolve
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import ngsolve.comp
import netgen.libngpy._meshing
import n
import ngstd
import ngsolve.fem
import g
import comp
import VorB
import la
import bla
import e
import t
import pyngcore
import fem
__all__  = [
"BFI",
"BSpline",
"BaseMatrix",
"BaseVector",
"BilinearForm",
"BitArray",
"BlockMatrix",
"BlockVector",
"COUPLING_TYPE",
"CoefficientFunction",
"FESpace",
"ConstEBEMatrix",
"ContactBoundary",
"Discontinuous",
"ET",
"ElementId",
"Embedding",
"Compress",
"FacetFESpace",
"FacetSurface",
"GridFunction",
"H1",
"HCurl",
"HCurlCurl",
"HCurlDiv",
"HDiv",
"HDivDiv",
"HDivDivSurface",
"HDivSurface",
"IdentityMatrix",
"IntRange",
"IntegrationRule",
"L2",
"LFI",
"LinearForm",
"Mesh",
"Preconditioner",
"MultiVector",
"NodalFESpace",
"NodeId",
"NormalFacetFESpace",
"NumProc",
"NumberSpace",
"ORDER_POLICY",
"PARALLEL_STATUS",
"PDE",
"Parameter",
"ParameterC",
"Periodic",
"PermutationMatrix",
"MultiGridPreconditioner",
"Projector",
"Region",
"SurfaceL2",
"TangentialFacetFESpace",
"TaskManager",
"Timer",
"Timing",
"VTKOutput",
"Variation",
"VectorFacetFESpace",
"VectorFacetSurface",
"VectorH1",
"VectorL2",
"VectorNodalFESpace",
"VectorSurfaceL2",
"ArnoldiSolver",
"BVP",
"BlockBFI",
"BlockLFI",
"BoundaryFromVolumeCF",
"CGSolver",
"CacheCF",
"CalcFlux",
"Cof",
"CompoundBFI",
"CompoundLFI",
"CompressCompound",
"Conj",
"CreateVVector",
"Cross",
"Det",
"Draw",
"DrawFlux",
"GMRESSolver",
"Grad",
"Id",
"IfPos",
"InnerProduct",
"Integrate",
"Interpolate",
"Inv",
"MPI_Init",
"Matrix",
"Norm",
"Normalize",
"OuterProduct",
"ParallelMatrix",
"PatchwiseSolve",
"PyCof",
"PyCross",
"PyDet",
"PyId",
"PyInv",
"PySkew",
"PySym",
"PyTrace",
"QMRSolver",
"Redraw",
"SetHeapSize",
"SetNumThreads",
"SetTestoutFile",
"SetVisualization",
"Skew",
"Sym",
"SymbolicBFI",
"SymbolicEnergy",
"SymbolicLFI",
"TimeFunction",
"Timers",
"Trace",
"Vector",
"VoxelCoefficient",
"_add_flags_doc",
"_jupyter_nbextension_paths",
"acos",
"asin",
"atan",
"atan2",
"ceil",
"cos",
"cosh",
"curl",
"div",
"exp",
"floor",
"grad",
"log",
"pow",
"sin",
"sinh",
"sqrt",
"tan",
"bla",
"bvp",
"comp",
"eigenvalues",
"fem",
"krylovspace",
"la",
"ngsolve",
"ngstd",
"nonlinearsolvers",
"pml",
"solve",
"solvers",
"timing",
"utils",
"BBBND",
"BBND",
"BND",
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
"VOL",
"__version__",
"ds",
"dx",
"mpi_world",
"ngsglobals",
"specialcf",
"x",
"y",
"z"
]
class BFI():
    pass
class BSpline():
    pass
class BaseMatrix():
    pass
class BaseVector():
    pass
class BilinearForm():
    pass
class BitArray():
    pass
class BlockMatrix():
    pass
class BlockVector():
    pass
class COUPLING_TYPE():
    pass
class CoefficientFunction():
    pass
class FESpace():
    pass
class ConstEBEMatrix():
    pass
class ContactBoundary():
    pass
class Discontinuous():
    pass
class ET():
    pass
class ElementId():
    pass
class Embedding():
    pass
class Compress():
    pass
class FacetFESpace():
    pass
class FacetSurface():
    pass
class GridFunction():
    pass
class H1():
    pass
class HCurl():
    pass
class HCurlCurl():
    pass
class HCurlDiv():
    pass
class HDiv():
    pass
class HDivDiv():
    pass
class HDivDivSurface():
    pass
class HDivSurface():
    pass
class IdentityMatrix():
    pass
class IntRange():
    pass
class IntegrationRule():
    pass
class L2():
    pass
class LFI():
    pass
class LinearForm():
    pass
class Mesh():
    pass
class Preconditioner():
    pass
class MultiVector():
    pass
class NodalFESpace():
    pass
class NodeId():
    pass
class NormalFacetFESpace():
    pass
class NumProc():
    pass
class NumberSpace():
    pass
class ORDER_POLICY():
    pass
class PARALLEL_STATUS():
    pass
class PDE():
    pass
class Parameter():
    pass
class ParameterC():
    pass
class Periodic():
    pass
class PermutationMatrix():
    pass
class MultiGridPreconditioner():
    pass
class Projector():
    pass
class Region():
    pass
class SurfaceL2():
    pass
class TangentialFacetFESpace():
    pass
class TaskManager():
    pass
class Timer():
    pass
class Timing():
    pass
class VTKOutput():
    pass
class Variation():
    pass
class VectorFacetFESpace():
    pass
class VectorFacetSurface():
    pass
class VectorH1():
    pass
class VectorL2():
    pass
class VectorNodalFESpace():
    pass
class VectorSurfaceL2():
    pass
def ArnoldiSolver(mata: BaseMatrix, matm: BaseMatrix, freedofs: pyngcore.BitArray, vecs: list, shift: complex = <ngstd.DummyArgument>) -> bla.VectorC:
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
def BVP(bf: comp.BilinearForm, lf: comp.LinearForm, gf: comp.GridFunction, pre: comp.Preconditioner, maxsteps: int = 100, prec: float = 1e-08) -> comp.NumProc:
    """
    Solves the given boundary value problem: bf * gf = lf, non homogeneous boundary conditions
    on gf are respected (they must be set in advance). If eliminate_internal is set for the
    bf, then static condensation of inner bubbles is used.

    Parameters:

    bf : ngsolve.comp.BilinearForm
      input bilinear form as the right hand side of the equation

    lf : ngsolve.comp.LinearForm
      input linear form as the left hand side of the equation

    gf : ngsolve.comp.GridFunction
      input GridFunction where the solution is saved

    pre : ngsolve.comp.Preconditioner
      input Preconditioner for the problem

    maxsteps : int
      input maximal steps. After the maximal step is reached, the computations stop.

    prec : float
      input precision of the residuum. if it is reached the computations stop.
    """
def BlockBFI(bfi: fem.BFI = 0, dim: int = 2, comp: int = 0) -> ngfem::BlockBilinearFormIntegrator:
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
def BlockLFI(lfi: fem.LFI = 0, dim: int = 2, comp: int = 0) -> fem.LFI:
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
def BoundaryFromVolumeCF(vol_cf: fem.CoefficientFunction) -> fem.CoefficientFunction:
    """
    Allows the evaluation of volumetric functions on the boundary.

    When evaluated on a boundary element, this function searches for the associated
    volume element, transforms the local coordinates, and evaluates the function in the
    volume. A typical use case is to visualize L2-functions, or mechanical stresses at
    the boundary.

    It is different from the boundary Trace()-operator. The trace provides a function
    which is defined by boundary degrees of freedom only. E.g. the trace of an H(div)
    function is only the normal component, while the BoundaryFromVolumeCF gives the
    whole function. Obviously, the Trace() function is cheaper to evaluate.

    If called on an interface, it evaluates from one side (which one is not specified).
    If the function is only defined on one side, this side will be taken. One can use
    a domain-wise CF to define a function only locally:
    uloc = CoefficientFunction( [None, None, u, None] )
    """
def CGSolver(mat: BaseMatrix, pre: BaseMatrix, complex: bool = False, printrates: bool = True, precision: float = 1e-08, maxsteps: int = 200) -> la.KrylovSpaceSolver:
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
def CacheCF(cf: fem.CoefficientFunction) -> fem.CoefficientFunction:
    pass
def CalcFlux(pde: comp.PDE, bf: comp.BilinearForm, gf: comp.GridFunction, flux: comp.GridFunction, applyd: bool = False) -> comp.NumProc:
    """
    Calculate Flux

    Parameters:

    pde : ngsolve.comp.PDE
      input pde

    bf : ngsolve.comp.BilinearForm
      input bilinear form

    gf : ngsolve.comp.GridFunction
      input GridFunction where the solution is saved

    flux : ngsolve.comp.GridFunction
      input GridFunction where the flux is saved

    applyd : bool
      input applyd
    """
def Cof(arg0: fem.CoefficientFunction) -> fem.CoefficientFunction:
    pass
def CompoundBFI(bfi: fem.BFI = 0, comp: int = 0) -> ngfem::CompoundBilinearFormIntegrator:
    """
    Compound Bilinear Form Integrator

    Parameters:

    bfi : ngsolve.fem.BFI
      input bilinear form integrator

    comp : int
      input component
    """
def CompoundLFI(lfi: fem.LFI = 0, comp: int = 0) -> fem.LFI:
    """
    Compound Linear Form Integrator

    Parameters:

    lfi : ngsolve.fem.LFI
      input linear form integrator

    comp : int
      input component
    """
def CompressCompound(fespace: comp.FESpace, active_dofs: object = <ngstd.DummyArgument>) -> comp.FESpace:
    pass
def Conj(arg0: fem.CoefficientFunction) -> fem.CoefficientFunction:
    """
    complex-conjugate
    """
def CreateVVector(size: int, complex: bool = False, entrysize: int = 1) -> ngla::BaseVector:
    pass
def Cross(arg0: fem.CoefficientFunction, arg1: fem.CoefficientFunction) -> fem.CoefficientFunction:
    pass
def Det(arg0: fem.CoefficientFunction) -> fem.CoefficientFunction:
    pass
@overload
def Draw(gf: comp.GridFunction, sd: int = 2, autoscale: bool = True, min: float = 0.0, max: float = 1.0, **kwargs) -> None:
    """
    Parameters:

    cf : ngsolve.comp.CoefficientFunction
      input CoefficientFunction to draw

    mesh : ngsolve.comp.Mesh
      input mesh

    name : string
      input name

    sd : int
      input subdivisions

    autoscale : bool
      input autscale

    min : float
      input minimum value. Need autoscale = false

    max : float
      input maximum value. Need autoscale = false

    draw_vol : bool
      input draw volume

    draw_surf : bool
      input draw surface




    Parameters:

    gf : ngsolve.comp.GridFunction
      input GridFunction to draw

    sd : int
      input subdivisions

    autoscale : bool
      input autscale

    min : float
      input minimum value. Need autoscale = false

    max : float
      input maximum value. Need autoscale = false
    """
@overload
def Draw(cf: fem.CoefficientFunction, mesh: comp.Mesh, name: str, sd: int = 2, autoscale: bool = True, min: float = 0.0, max: float = 1.0, draw_vol: bool = True, draw_surf: bool = True, reset: bool = False, **kwargs) -> None:
    pass
@overload
def Draw(mesh: comp.Mesh, **kwargs) -> None:
    pass
@overload
def Draw(arg0: object) -> None:
    pass
def DrawFlux(bf: comp.BilinearForm, gf: comp.GridFunction, label: str = 'flux', applyd: bool = False, useall: bool = False) -> comp.NumProc:
    """
    draw Flux

    Parameters:


    bf : ngsolve.comp.BilinearForm
      input bilinear form

    gf : ngsolve.comp.GridFunction
      input GridFunction where the flux is saved

    label : string
      input name of the flux

    applyd : bool
      input applyd

    useall : bool
      input useall
    """
def GMRESSolver(mat: BaseMatrix, pre: BaseMatrix, printrates: bool = True, precision: float = 1e-08, maxsteps: int = 200) -> la.KrylovSpaceSolver:
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
def Id(arg0: int) -> fem.CoefficientFunction:
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
def InnerProduct(x: object, y: object, **kwargs) -> object:
    """
    Compute InnerProduct
    """
@overload
def Integrate(cf: fem.CoefficientFunction, mesh: Union[comp.Mesh, comp.Region], VOL_or_BND: comp.VorB = VorB.VOL, order: int = 5, definedon: comp.Region = None, region_wise: bool = False, element_wise: bool = False) -> object:
    """
    Parameters
    ----------

    cf: ngsolve.CoefficientFunction
      Function to be integrated. Can be vector valued, then the result is an array. If you want to integrate
      a lot of functions on the same domain, it will be faster to put them into a vector valued function,
      NGSolve will then be able to use parallelization and SIMD vectorization more efficiently.

    mesh: ngsolve.Mesh
      The mesh to be integrated on.

    VOL_or_BND: ngsolve.VorB = VOL
      Co-dimension to be integrated on. Historically this could be volume (VOL) or boundary (BND). If your mesh
      contains co-dim 2 elements this can now be BBND (edges in 3d) as well.

    order: int = 5
      Integration order, polynomials up to this order will be integrated exactly.

    definedon: ngsolve.Region
      Region to be integrated on. Such region can be created with mesh.Boundaries('bcname') or mesh.Materials('matname')
      it will overwrite the VOL_or_BND argument if given.

    region_wise: bool = False
      Integrates region wise on the co-dimension given by VOL_or_BND. Returns results as an array, matching the array
      returned by mesh.GetMaterials() or mesh.GetBoundaries(). Does not support vector valued CoefficientFunctions.

    element_wise: bool = False
      Integrates element wise and returns result in a list. This is typically used for local error estimators.
      Does not support vector valued CoefficientFunctions
    """
@overload
def Integrate(igls: comp.SumOfIntegrals, mesh: comp.Mesh, element_wise: bool = False) -> object:
    pass
def Interpolate(cf: fem.CoefficientFunction, space: comp.FESpace, bonus_intorder: int = 0) -> fem.CoefficientFunction:
    """
    Interpolate a CoefficientFunction into the finite element space.
    The interpolation is canonical interpolation using dual shapes.
    The result is a CoefficientFunction.
    Interpolation is done on the fly for each element, no global GridFunction is allocated.
    """
def Inv(arg0: fem.CoefficientFunction) -> fem.CoefficientFunction:
    pass
def MPI_Init() -> netgen.libngpy._meshing.MPI_Comm:
    pass
@overload
def Matrix(arg0: List[List[float]]) -> bla.MatrixD:
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
def Matrix(arg0: List[List[complex]]) -> bla.MatrixC:
    pass
@overload
def Matrix(height: int, width: int, complex: bool = False) -> object:
    pass
def Norm(x: object) -> object:
    """
    Compute Norm
    """
@overload
def ParallelMatrix(mat: object = None, row_pardofs: object = None, col_pardofs: object = None, op: object = None) -> None:
    pass
@overload
def ParallelMatrix(mat: object = None, pardofs: object = None, op: object = None) -> None:
    pass
def PatchwiseSolve(bf: comp.SumOfIntegrals, lf: comp.SumOfIntegrals, gf: comp.GridFunction) -> None:
    pass
def QMRSolver(mat: BaseMatrix, pre: BaseMatrix, printrates: bool = True, precision: float = 1e-08, maxsteps: int = 200) -> la.KrylovSpaceSolver:
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
def SetHeapSize(size: int) -> None:
    """
    Set a new heapsize.

    Parameters:

    size : int
      input heap size
    """
def SetNumThreads(threads: int) -> None:
    """
    Set number of threads

    Parameters:

    threads : int
      input number of threads
    """
def SetTestoutFile(file: str) -> None:
    """
    Enable some logging into file with given filename

    Parameters:

    file : string
      input file name
    """
def SetVisualization(deformation: object = <ngstd.DummyArgument>, min: object = <ngstd.DummyArgument>, max: object = <ngstd.DummyArgument>, clipnormal: object = <ngstd.DummyArgument>, clipping: object = <ngstd.DummyArgument>) -> None:
    """
    Set visualization options

    Parameters:

    deformation : object
      input deformation

    min : object
      input min

    max : object
      input max

    clipnormal : object
      input clipnormal

    clipping : object
      input clipping
    """
def Skew(arg0: fem.CoefficientFunction) -> fem.CoefficientFunction:
    pass
def Sym(arg0: fem.CoefficientFunction) -> fem.CoefficientFunction:
    pass
def SymbolicBFI(form: fem.CoefficientFunction, VOL_or_BND: comp.VorB = VorB.VOL, element_boundary: bool = False, skeleton: bool = False, definedon: Optional[Union[comp.Region, list]] = None, intrule: fem.IntegrationRule = <fem.IntegrationRule object at 0x0000026DEE4255F0>, bonus_intorder: int = 0, definedonelements: pyngcore.BitArray = None, simd_evaluate: bool = True, element_vb: comp.VorB = VorB.VOL, geom_free: bool = False, deformation: comp.GridFunction = None) -> fem.BFI:
    """
    A symbolic bilinear form integrator, where test and trial functions, CoefficientFunctions, etc. can be used to formulate PDEs in a symbolic way.

    Parameters:

    form : ngsolve.fem.CoefficientFunction
      input the symbolic right hand side form

    VOL_or_BND : ngsolve.comp.VorB
      input VOL, BND, BBND, ...

    element_boundary : bool
      input element_boundary. True -> iterates over all element boundaries, but uses volume transformations

    skeleton : bool
      input skeleton. True -> iterates over all faces, but uses volume transformations

    definedon : object
      input definedon region

    intrule : ngsolve.fem.IntegrationRule
      input integration rule

    bonus_intorder : int
      input additional integration order

    definedonelements : object
      input definedonelements

    simd_evaluate : bool
      input simd_evaluate. True -> tries to use SIMD for faster evaluation

    element_vb : ngsolve.comp.VorB
      input element_vb. Used for skeleton formulation. VOL -> interior faces, BND -> boundary faces

    deformation : ngsolve.comp.GridFunction
      input GridFunction to transform/deform the bilinear form with
    """
def SymbolicEnergy(form: fem.CoefficientFunction, VOL_or_BND: comp.VorB = VorB.VOL, definedon: object = <ngstd.DummyArgument>, element_boundary: bool = False, bonus_intorder: int = 0, definedonelements: object = <ngstd.DummyArgument>, simd_evaluate: bool = True, element_vb: comp.VorB = VorB.VOL, deformation: comp.GridFunction = None) -> fem.BFI:
    """
    A symbolic energy form integrator, where test and trial functions, CoefficientFunctions, etc. can be used to formulate PDEs in a symbolic way.

    Parameters:

    form : ngsolve.fem.CoefficientFunction
      input the symbolic right hand side form

    VOL_or_BND : ngsolve.comp.VorB
      input VOL, BND, BBND, ...

    definedon : object
      input definedon region

    element_boundary : bool
      input element_boundary. True -> iterates over all element boundaries, but uses volume transformations

    bonus_intorder : int
      input additional integration order

    definedonelements : object
      input definedonelements

    simd_evaluate : bool
      input simd_evaluate. True -> tries to use SIMD for faster evaluation

    element_vb : ngsolve.fem.VorB
      input eleemnt VorB

    deformation : ngsolve.comp.GridFunction
      input GridFunction to transform/deform the bilinear form with
    """
def SymbolicLFI(form: fem.CoefficientFunction, VOL_or_BND: comp.VorB = VorB.VOL, element_boundary: bool = False, skeleton: bool = False, definedon: Optional[Union[comp.Region, list]] = None, intrule: fem.IntegrationRule = <fem.IntegrationRule object at 0x0000026DEE4256F0>, bonus_intorder: int = 0, definedonelements: pyngcore.BitArray = None, simd_evaluate: bool = True, element_vb: comp.VorB = VorB.VOL, deformation: comp.GridFunction = None) -> fem.LFI:
    """
    A symbolic linear form integrator, where test and trial functions, CoefficientFunctions, etc. can be used to formulate right hand sides in a symbolic way.

    Parameters:

    form : ngsolve.fem.CoefficientFunction
      input the symbolic right hand side form

    VOL_or_BND : ngsolve.comp.VorB
      input VOL, BND, BBND, ...

    element_boundary : bool
      input element_boundary. True -> iterates over all element boundaries, but uses volume transformations

    skeleton : bool
      input skeleton. True -> iterates over all faces, but uses volume transformations

    definedon : object
      input definedon region

    intrule : ngsolve.fem.IntegrationRule
      input integration rule

    bonus_intorder : int
      input additional integration order

    definedonelements : object
      input BitArray that marks all elements or facets (for skeleton-integrators) that the integrator is applied on

    simd_evaluate : bool
      input simd_evaluate. True -> tries to use SIMD for faster evaluation

    element_vb : ngsolve.fem.VorB
      input element VorB

    deformation : ngsolve.comp.GridFunction
      input GridFunction to transform/deform the linear form with
    """
def Timers() -> list:
    """
    Returns list of timers
    """
def Trace(arg0: fem.CoefficientFunction) -> fem.CoefficientFunction:
    pass
@overload
def Vector(length: int, complex: bool = False) -> object:
    """
    Parameters:

    length : int
      input length

    complex : bool
      input complex values
    """
@overload
def Vector(arg0: List[float]) -> bla.VectorD:
    pass
@overload
def Vector(arg0: List[complex]) -> bla.VectorC:
    pass
def VoxelCoefficient(start: tuple, end: tuple, values: array, linear: bool = True, trafocf: object = <ngstd.DummyArgument>) -> fem.CoefficientFunction:
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
BBBND: ngsolve.comp.VorB # value = VorB.BBBND
BBND: ngsolve.comp.VorB # value = VorB.BBND
BND: ngsolve.comp.VorB # value = VorB.BND
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
VOL: ngsolve.comp.VorB # value = VorB.VOL
__version__ = '6.2.2007'
ds: ngsolve.comp.DifferentialSymbol # value = <ngsolve.comp.DifferentialSymbol object at 0x0000026DEE18E230>
dx: ngsolve.comp.DifferentialSymbol # value = <ngsolve.comp.DifferentialSymbol object at 0x0000026DEE18EB30>
mpi_world: netgen.libngpy._meshing.MPI_Comm # value = <netgen.libngpy._meshing.MPI_Comm object at 0x0000026DEE40EA30>
ngsglobals: ngsolve.comp.GlobalVariables # value = <ngsolve.comp.GlobalVariables object at 0x0000026DEE409B70>
specialcf: ngsolve.fem.SpecialCFCreator # value = <ngsolve.fem.SpecialCFCreator object at 0x0000026DEE186470>
x: ngsolve.fem.CoefficientFunction # value = <ngsolve.fem.CoefficientFunction object at 0x0000026DEE419F48>
y: ngsolve.fem.CoefficientFunction # value = <ngsolve.fem.CoefficientFunction object at 0x0000026DEE419FA8>
z: ngsolve.fem.CoefficientFunction # value = <ngsolve.fem.CoefficientFunction object at 0x0000026DEE430048>
