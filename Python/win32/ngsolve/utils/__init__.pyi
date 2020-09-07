import ngsolve.utils
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import ngsolve.comp
import netgen.libngpy._meshing
import m
import n
import s
import pyngcore
import l
import ngsolve.fem
import v
import c
import g
import p
import VorB
import e
import ngsolve.ngstd
import .
import ngsolve.la
import o
__all__  = [
"FlatArray_enum ngcomp::COUPLING_TYPE_S",
"BFI",
"BSpline",
"BaseMappedIntegrationPoint",
"NGS_Object",
"COUPLING_TYPE",
"CoefficientFunction",
"GridFunction",
"FESpace",
"Compress",
"ContactBoundary",
"DifferentialOperator",
"DifferentialSymbol",
"Discontinuous",
"ProxyFunction",
"ET",
"ElementId",
"ElementRange",
"ElementTopology",
"ElementTransformation",
"CompoundFESpace",
"Ngs_Element",
"FESpaceElementRange",
"FacetFESpace",
"FacetSurface",
"FiniteElement",
"Array_enum ngcomp::COUPLING_TYPE_S",
"GlobalVariables",
"ComponentGridFunction",
"GridFunctionC",
"GridFunctionCoefficientFunction",
"GridFunctionD",
"H1",
"HCurl",
"HCurlCurl",
"HCurlDiv",
"HCurlFE",
"HDiv",
"HDivDiv",
"HDivDivFE",
"HDivDivSurface",
"HDivFE",
"HDivSurface",
"Integral",
"IntegrationPoint",
"IntegrationRule",
"InterpolateProxy",
"L2",
"LFI",
"LinearForm",
"Mesh",
"NodeId",
"MeshNodeRange",
"MeshPoint",
"MixedFE",
"Preconditioner",
"BilinearForm",
"NODE_TYPE",
"FESpaceElement",
"NodalFESpace",
"MeshNode",
"NodeRange",
"NormalFacetFESpace",
"NumProc",
"NumberSpace",
"ORDER_POLICY",
"PDE",
"Parameter",
"ParameterC",
"Periodic",
"MultiGridPreconditioner",
"Prolongation",
"DualProxyFunction",
"PyNumProc",
"Region",
"Reorder",
"ScalarFE",
"SpecialCFCreator",
"SumOfIntegrals",
"SurfaceL2",
"SymbolTable_D",
"SymbolTable_sp_D",
"SymbolTable_sp_class ngcomp::BilinearForm",
"SymbolTable_sp_class ngcomp::FESpace",
"SymbolTable_sp_class ngcomp::GridFunction",
"SymbolTable_sp_class ngcomp::LinearForm",
"SymbolTable_sp_class ngcomp::NumProc",
"SymbolTable_sp_class ngcomp::Preconditioner",
"SymbolTable_sp_class ngfem::CoefficientFunction",
"TangentialFacetFESpace",
"Timer",
"VTKOutput",
"Variation",
"VectorFacetFESpace",
"VectorFacetSurface",
"VectorH1",
"VectorL2",
"VectorNodalFESpace",
"VectorSurfaceL2",
"VorB",
" ",
"BlockBFI",
"BlockLFI",
"BndElementId",
"BoundaryFromVolumeCF",
"CacheCF",
"Cof",
"CompoundBFI",
"CompoundLFI",
"CompressCompound",
"Conj",
"ConstantCF",
"ConvertOperator",
"CoordCF",
"Cross",
"Det",
"DomainConstantCF",
"GenerateL2ElementCode",
"Grad",
"H1FE",
"Id",
"IfPos",
"Integrate",
"Interpolate",
"Inv",
"KSpaceCoeffs",
"L2FE",
"Laplace",
"LoggingCF",
"MPI_Init",
"Mass",
"Neumann",
"Norm",
"Normalize",
"OuterProduct",
"PatchwiseSolve",
"Prolongate",
"ProlongateCoefficientFunction",
"PyCof",
"PyCross",
"PyDet",
"PyId",
"PyInv",
"PySkew",
"PySym",
"PyTrace",
"SetHeapSize",
"SetPMLParameters",
"SetTestoutFile",
"Skew",
"Source",
"Sym",
"SymbolicBFI",
"SymbolicEnergy",
"SymbolicLFI",
"SymbolicTPBFI",
"TensorProductFESpace",
"TensorProductIntegrate",
"TimeFunction",
"Trace",
"Transfer2StdMesh",
"VectorFacet",
"VoxelCoefficient",
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
"ds",
"dx",
"ngsglobals",
"specialcf",
"x",
"y",
"z"
]
class FlatArray_enum ngcomp::COUPLING_TYPE_S():
    pass
class BFI():
    pass
class BSpline():
    pass
class BaseMappedIntegrationPoint():
    pass
class NGS_Object():
    pass
class COUPLING_TYPE():
    pass
class CoefficientFunction():
    pass
class GridFunction():
    pass
class FESpace():
    pass
class Compress():
    pass
class ContactBoundary():
    pass
class DifferentialOperator():
    pass
class DifferentialSymbol():
    pass
class Discontinuous():
    pass
class ProxyFunction():
    pass
class ET():
    pass
class ElementId():
    pass
class ElementRange():
    pass
class ElementTopology():
    pass
class ElementTransformation():
    pass
class CompoundFESpace():
    pass
class Ngs_Element():
    pass
class FESpaceElementRange():
    pass
class FacetFESpace():
    pass
class FacetSurface():
    pass
class FiniteElement():
    pass
class Array_enum ngcomp::COUPLING_TYPE_S():
    pass
class GlobalVariables():
    pass
class ComponentGridFunction():
    pass
class GridFunctionC():
    pass
class GridFunctionCoefficientFunction():
    pass
class GridFunctionD():
    pass
class H1():
    pass
class HCurl():
    pass
class HCurlCurl():
    pass
class HCurlDiv():
    pass
class HCurlFE():
    pass
class HDiv():
    pass
class HDivDiv():
    pass
class HDivDivFE():
    pass
class HDivDivSurface():
    pass
class HDivFE():
    pass
class HDivSurface():
    pass
class Integral():
    pass
class IntegrationPoint():
    pass
class IntegrationRule():
    pass
class InterpolateProxy():
    pass
class L2():
    pass
class LFI():
    pass
class LinearForm():
    pass
class Mesh():
    pass
class NodeId():
    pass
class MeshNodeRange():
    pass
class MeshPoint():
    pass
class MixedFE():
    pass
class Preconditioner():
    pass
class BilinearForm():
    pass
class NODE_TYPE():
    pass
class FESpaceElement():
    pass
class NodalFESpace():
    pass
class MeshNode():
    pass
class NodeRange():
    pass
class NormalFacetFESpace():
    pass
class NumProc():
    pass
class NumberSpace():
    pass
class ORDER_POLICY():
    pass
class PDE():
    pass
class Parameter():
    pass
class ParameterC():
    pass
class Periodic():
    pass
class MultiGridPreconditioner():
    pass
class Prolongation():
    pass
class DualProxyFunction():
    pass
class PyNumProc():
    pass
class Region():
    pass
class Reorder():
    pass
class ScalarFE():
    pass
class SpecialCFCreator():
    pass
class SumOfIntegrals():
    pass
class SurfaceL2():
    pass
class SymbolTable_D():
    pass
class SymbolTable_sp_D():
    pass
class SymbolTable_sp_class ngcomp::BilinearForm():
    pass
class SymbolTable_sp_class ngcomp::FESpace():
    pass
class SymbolTable_sp_class ngcomp::GridFunction():
    pass
class SymbolTable_sp_class ngcomp::LinearForm():
    pass
class SymbolTable_sp_class ngcomp::NumProc():
    pass
class SymbolTable_sp_class ngcomp::Preconditioner():
    pass
class SymbolTable_sp_class ngfem::CoefficientFunction():
    pass
class TangentialFacetFESpace():
    pass
class Timer():
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
class VorB():
    pass
def  (x: object) -> object:
    """
     (x: object) -> object

    Passes value through
    """
def BlockBFI(bfi: ngsolve.fem.BFI = 0, dim: int = 2, comp: int = 0) -> ngfem::BlockBilinearFormIntegrator:
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
def BlockLFI(lfi: ngsolve.fem.LFI = 0, dim: int = 2, comp: int = 0) -> ngsolve.fem.LFI:
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
def BndElementId(nr: int) -> ngsolve.comp.ElementId:
    """
    Creates an element-id for a boundary element

    Parameters:

    nr : int
      input Bnd element number
    """
def BoundaryFromVolumeCF(vol_cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
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
def CacheCF(cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
    pass
def Cof(arg0: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
    pass
def CompoundBFI(bfi: ngsolve.fem.BFI = 0, comp: int = 0) -> ngfem::CompoundBilinearFormIntegrator:
    """
    Compound Bilinear Form Integrator

    Parameters:

    bfi : ngsolve.fem.BFI
      input bilinear form integrator

    comp : int
      input component
    """
def CompoundLFI(lfi: ngsolve.fem.LFI = 0, comp: int = 0) -> ngsolve.fem.LFI:
    """
    Compound Linear Form Integrator

    Parameters:

    lfi : ngsolve.fem.LFI
      input linear form integrator

    comp : int
      input component
    """
def CompressCompound(fespace: ngsolve.comp.FESpace, active_dofs: object = <ngsolve.ngstd.DummyArgument>) -> ngsolve.comp.FESpace:
    pass
def Conj(arg0: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
    """
    complex-conjugate
    """
def ConvertOperator(spacea: ngsolve.comp.FESpace, spaceb: ngsolve.comp.FESpace, trial_proxy: ngsolve.comp.ProxyFunction = None, trial_cf: ngsolve.fem.CoefficientFunction = None, definedon: Optional[ngsolve.comp.Region] = None, vb: ngsolve.comp.VorB = VorB.VOL, range_dofs: pyngcore.BitArray = None, localop: bool = False, parmat: bool = True, use_simd: bool = True, bonus_intorder_ab: int = 0, bonus_intorder_bb: int = 0) -> ngsolve.la.BaseMatrix:
    """
    A conversion operator between FESpaces. Embedding if spacea is a subspace of spaceb, otherwise an interpolation operator defined by element-wise application of dual shapes (and averaging between elements).

    Parameters:

    spacea: ngsolve.comp.FESpace
      the origin space

    spaceb: ngsolve.comp.FESpace
      the goal space

    trial_proxy: ngsolve.comp.ProxyFunction
      (optional) Must be a trial-proxy on spacea. If given, instead of a FE-function funca from spacea, the operator converts trial_proxy(funca) to spaceb.

    trial_proxy: ngsolve.comp.CoefficientFunction
      (optional) Same as trial_proxy, but takes any CoefficientFunction. Use at your own peril.

    definedon: object
      what part of the domain to restrict the operator to

    vb: ngsolve.comp.VorB
      what kind of co-dimension elements to convert on VOL, BND, BBND, ...

    range_dofs: ngsolve.ngstd.BitArray
      Projects out DOFs in the range where range_dofs are not set

    localop: bool
      True -> do not average across MPI boundaries. No effect for non MPI-paralell space. Use carefully!!

    parmat: bool
      If True, returns a ParallelMatrix for MPI-parallel spaces. If False, or for non MPI-parallel spaces, returns a local BaseMatrix.

    use_simd:
      False -> Do not use SIMD for setting up the Matrix. (for debugging purposes).

    bonus_intorder_ab/bb: int
      Bonus integration order for spacea/spaceb and spaceb/spaceb integrals. Can be useful for curved elements. Should only be necessary for
    spacea/spaceb integrals.
    """
def CoordCF(direction: int) -> ngfem::CoefficientFunction:
    """
    CoefficientFunction for x, y, z.

    Parameters:

    direction : int
      input direction
    """
def Cross(arg0: ngsolve.fem.CoefficientFunction, arg1: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
    pass
def Det(arg0: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
    pass
def GenerateL2ElementCode(arg0: int) -> str:
    pass
def H1FE(et: ngsolve.fem.ET, order: int) -> None:
    """
    Creates an H1 finite element of given geometric shape and polynomial order.

    Parameters:

    et : ngsolve.fem.ET
      input element type

    order : int
      input polynomial order
    """
def Id(arg0: int) -> ngsolve.fem.CoefficientFunction:
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
@overload
def Integrate(cf: ngsolve.fem.CoefficientFunction, mesh: Union[ngsolve.comp.Mesh, ngsolve.comp.Region], VOL_or_BND: ngsolve.comp.VorB = VorB.VOL, order: int = 5, definedon: ngsolve.comp.Region = None, region_wise: bool = False, element_wise: bool = False) -> object:
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
def Integrate(igls: ngsolve.comp.SumOfIntegrals, mesh: ngsolve.comp.Mesh, element_wise: bool = False) -> object:
    pass
def Interpolate(cf: ngsolve.fem.CoefficientFunction, space: ngsolve.comp.FESpace, bonus_intorder: int = 0) -> ngsolve.fem.CoefficientFunction:
    """
    Interpolate a CoefficientFunction into the finite element space.
    The interpolation is canonical interpolation using dual shapes.
    The result is a CoefficientFunction.
    Interpolation is done on the fly for each element, no global GridFunction is allocated.
    """
def Inv(arg0: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
    pass
def KSpaceCoeffs(arg0: ngsolve.comp.GridFunction, arg1: ngsolve.comp.GridFunction, arg2: float, arg3: float) -> None:
    pass
def L2FE(et: ngsolve.fem.ET, order: int) -> None:
    """
    Creates an L2 finite element of given geometric shape and polynomial order.

    Parameters:

    et : ngsolve.fem.ET
      input element type

    order : int
      input polynomial order
    """
def LoggingCF(cf: ngsolve.fem.CoefficientFunction, logfile: str = 'stdout') -> ngsolve.fem.CoefficientFunction:
    pass
def MPI_Init() -> netgen.libngpy._meshing.MPI_Comm:
    pass
def Norm(x: object) -> object:
    """
    Compute Norm
    """
def PatchwiseSolve(bf: ngsolve.comp.SumOfIntegrals, lf: ngsolve.comp.SumOfIntegrals, gf: ngsolve.comp.GridFunction) -> None:
    pass
def Prolongate(arg0: ngsolve.comp.GridFunction, arg1: ngsolve.comp.GridFunction) -> None:
    pass
def ProlongateCoefficientFunction(arg0: ngsolve.fem.CoefficientFunction, arg1: int, arg2: ngsolve.comp.FESpace) -> ngsolve.fem.CoefficientFunction:
    pass
def SetHeapSize(size: int) -> None:
    """
    Set a new heapsize.

    Parameters:

    size : int
      input heap size
    """
def SetPMLParameters(rad: float = 1, alpha: float = 1) -> None:
    """
    Parameters:

    rad : double
      input radius of PML

    alpha : double
      input damping factor of PML
    """
def SetTestoutFile(file: str) -> None:
    """
    Enable some logging into file with given filename

    Parameters:

    file : string
      input file name
    """
def Skew(arg0: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
    pass
def Sym(arg0: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
    pass
def SymbolicBFI(form: ngsolve.fem.CoefficientFunction, VOL_or_BND: ngsolve.comp.VorB = VorB.VOL, element_boundary: bool = False, skeleton: bool = False, definedon: Optional[Union[ngsolve.comp.Region, list]] = None, intrule: ngsolve.fem.IntegrationRule = <ngsolve.fem.IntegrationRule object at 0x0000026DEE4255F0>, bonus_intorder: int = 0, definedonelements: pyngcore.BitArray = None, simd_evaluate: bool = True, element_vb: ngsolve.comp.VorB = VorB.VOL, geom_free: bool = False, deformation: ngsolve.comp.GridFunction = None) -> ngsolve.fem.BFI:
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
def SymbolicEnergy(form: ngsolve.fem.CoefficientFunction, VOL_or_BND: ngsolve.comp.VorB = VorB.VOL, definedon: object = <ngsolve.ngstd.DummyArgument>, element_boundary: bool = False, bonus_intorder: int = 0, definedonelements: object = <ngsolve.ngstd.DummyArgument>, simd_evaluate: bool = True, element_vb: ngsolve.comp.VorB = VorB.VOL, deformation: ngsolve.comp.GridFunction = None) -> ngsolve.fem.BFI:
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
def SymbolicLFI(form: ngsolve.fem.CoefficientFunction, VOL_or_BND: ngsolve.comp.VorB = VorB.VOL, element_boundary: bool = False, skeleton: bool = False, definedon: Optional[Union[ngsolve.comp.Region, list]] = None, intrule: ngsolve.fem.IntegrationRule = <ngsolve.fem.IntegrationRule object at 0x0000026DEE4256F0>, bonus_intorder: int = 0, definedonelements: pyngcore.BitArray = None, simd_evaluate: bool = True, element_vb: ngsolve.comp.VorB = VorB.VOL, deformation: ngsolve.comp.GridFunction = None) -> ngsolve.fem.LFI:
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
def SymbolicTPBFI(form: ngsolve.fem.CoefficientFunction, VOL_or_BND: ngsolve.comp.VorB = VorB.VOL, element_boundary: bool = False, skeleton: bool = False, definedon: object = <ngsolve.ngstd.DummyArgument>) -> ngsolve.fem.BFI:
    pass
def TensorProductFESpace(spaces: list, flags: pyngcore.Flags = <pyngcore.Flags object at 0x0000026DEE425BB0>) -> ngsolve.comp.FESpace:
    pass
@overload
def TensorProductIntegrate(gftp: ngsolve.comp.GridFunction, gfx: ngsolve.comp.GridFunction, weight: ngsolve.fem.CoefficientFunction = None) -> None:
    pass
@overload
def TensorProductIntegrate(arg0: ngsolve.comp.GridFunction, arg1: list, arg2: ngsolve.fem.CoefficientFunction) -> float:
    pass
def Trace(arg0: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction:
    pass
@overload
def Transfer2StdMesh(arg0: ngsolve.fem.CoefficientFunction, arg1: ngsolve.comp.GridFunction) -> None:
    pass
@overload
def Transfer2StdMesh(gftp: ngsolve.comp.GridFunction, gfstd: ngsolve.comp.GridFunction) -> None:
    pass
def VoxelCoefficient(start: tuple, end: tuple, values: array, linear: bool = True, trafocf: object = <ngsolve.ngstd.DummyArgument>) -> ngsolve.fem.CoefficientFunction:
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
ds: ngsolve.comp.DifferentialSymbol # value = <ngsolve.comp.DifferentialSymbol object at 0x0000026DEE18E230>
dx: ngsolve.comp.DifferentialSymbol # value = <ngsolve.comp.DifferentialSymbol object at 0x0000026DEE18EB30>
ngsglobals: ngsolve.comp.GlobalVariables # value = <ngsolve.comp.GlobalVariables object at 0x0000026DEE409B70>
specialcf: ngsolve.fem.SpecialCFCreator # value = <ngsolve.fem.SpecialCFCreator object at 0x0000026DEE186470>
x: ngsolve.fem.CoefficientFunction # value = <ngsolve.fem.CoefficientFunction object at 0x0000026DEE419F48>
y: ngsolve.fem.CoefficientFunction # value = <ngsolve.fem.CoefficientFunction object at 0x0000026DEE419FA8>
z: ngsolve.fem.CoefficientFunction # value = <ngsolve.fem.CoefficientFunction object at 0x0000026DEE430048>
