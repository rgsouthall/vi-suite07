"""pybind comp"""
import ngsolve.comp
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import netgen.libngpy._meshing
import pml
import numpy
import VorB
import pyngcore
import ngsolve.ngstd
import ngsolve.la
import ngsolve.fem
__all__  = [
"FlatArray_enum ngcomp::COUPLING_TYPE_S",
"NGS_Object",
"COUPLING_TYPE",
"GridFunction",
"FESpace",
"Compress",
"ContactBoundary",
"DifferentialSymbol",
"Discontinuous",
"ProxyFunction",
"ElementId",
"ElementRange",
"CompoundFESpace",
"Ngs_Element",
"FESpaceElementRange",
"FacetFESpace",
"FacetSurface",
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
"HDiv",
"HDivDiv",
"HDivDivSurface",
"HDivSurface",
"Integral",
"InterpolateProxy",
"L2",
"LinearForm",
"Mesh",
"NodeId",
"MeshNodeRange",
"Preconditioner",
"BilinearForm",
"FESpaceElement",
"NodalFESpace",
"MeshNode",
"NodeRange",
"NormalFacetFESpace",
"NumProc",
"NumberSpace",
"ORDER_POLICY",
"PDE",
"Periodic",
"MultiGridPreconditioner",
"Prolongation",
"DualProxyFunction",
"PyNumProc",
"Region",
"Reorder",
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
"VTKOutput",
"Variation",
"VectorFacetFESpace",
"VectorFacetSurface",
"VectorH1",
"VectorL2",
"VectorNodalFESpace",
"VectorSurfaceL2",
"VorB",
"BndElementId",
"BoundaryFromVolumeCF",
"CompressCompound",
"ConvertOperator",
"Integrate",
"Interpolate",
"KSpaceCoeffs",
"MPI_Init",
"PatchwiseSolve",
"Prolongate",
"ProlongateCoefficientFunction",
"SetHeapSize",
"SetTestoutFile",
"SymbolicBFI",
"SymbolicEnergy",
"SymbolicLFI",
"SymbolicTPBFI",
"TensorProductFESpace",
"TensorProductIntegrate",
"Transfer2StdMesh",
"pml",
"BBBND",
"BBND",
"BND",
"VOL",
"ngsglobals"
]
class FlatArray_enum ngcomp::COUPLING_TYPE_S():
    def NumPy(self) -> object: ...
    def __getitem__(self, arg0: int) -> COUPLING_TYPE: ...
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: ...
    @overload
    def __setitem__(self, arg0: int, arg1: COUPLING_TYPE) -> COUPLING_TYPE: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: COUPLING_TYPE) -> None: ...
    def __str__(self) -> str: ...
    pass
class NGS_Object():
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    pass
class COUPLING_TYPE():
    """
    Enum specifying the coupling type of a degree of freedom, each dof is
    either UNUSED_DOF, LOCAL_DOF, INTERFACE_DOF or WIREBASKET_DOF, other values
    are provided as combinations of these:

    UNUSED_DOF: Dof is not used, i.e the minion dofs in a Periodic finite
        element space.

    LOCAL_DOF: Inner degree of freedom, will be eliminated by static
        condensation and reconstructed afterwards.

    HIDDEN_DOF: Inner degree of freedom, that will be eliminated by static
        condensation and *not* reconstruced afterwards(spares some entries).
        Note: 
         * without static condensation a HIDDEN_DOF is treated as any other
           DOF, e.g. as a LOCAL_DOF
         * To a HIDDEN_DOF the r.h.s. vector must have zero entries.
         * When static condensation is applied (eliminate_hidden/
           eliminate_internal) the block corresponding to HIDDEN_DOFs
           has to be invertible.

    CONDENSABLE_DOF: Inner degree of freedom, that will be eliminated by static
        condensation (LOCAL_DOF or HIDDEN_DOF)

    INTERFACE_DOF: Degree of freedom between two elements, these will not be
        eliminated by static condensation, but not be put into the wirebasket
        system for i.e. a bddc Preconditioner.

    NONWIREBASKET_DOF: Either a LOCAL_DOF or an INTERFACE_DOF

    WIREBASKET_DOF: Degree of freedom coupling with many elements (more than
        one). These will be put into the system for a bddc preconditioner.
        The HCurl space also treats degrees of freedom of badly shaped
        elements as WIREBASKET_DOFs.

    EXTERNAL_DOF: Either INTERFACE_DOF or WIREBASKET_DOF

    VISIBLE_DOF: not UNUSED_DOF or HIDDEN_DOF

    ANY_DOF: Any used dof (LOCAL_DOF or INTERFACE_DOF or WIREBASKET_DOF)



    Members:

      UNUSED_DOF

      HIDDEN_DOF

      LOCAL_DOF

      CONDENSABLE_DOF

      INTERFACE_DOF

      NONWIREBASKET_DOF

      WIREBASKET_DOF

      EXTERNAL_DOF

      VISIBLE_DOF

      ANY_DOF
    """
    def __init__(self, arg0: int) -> None: ...
    def __int__(self) -> int: ...
    @property
    def name(self) -> str:
        """
        (self: handle) -> str

        :type: str
        """
    ANY_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.ANY_DOF
    CONDENSABLE_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.CONDENSABLE_DOF
    EXTERNAL_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.EXTERNAL_DOF
    HIDDEN_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.HIDDEN_DOF
    INTERFACE_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.INTERFACE_DOF
    LOCAL_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.LOCAL_DOF
    NONWIREBASKET_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.NONWIREBASKET_DOF
    UNUSED_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.UNUSED_DOF
    VISIBLE_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.VISIBLE_DOF
    WIREBASKET_DOF: ngsolve.comp.COUPLING_TYPE # value = COUPLING_TYPE.WIREBASKET_DOF
    __entries: dict # value = {'UNUSED_DOF': (COUPLING_TYPE.UNUSED_DOF, None), 'HIDDEN_DOF': (COUPLING_TYPE.HIDDEN_DOF, None), 'LOCAL_DOF': (COUPLING_TYPE.LOCAL_DOF, None), 'CONDENSABLE_DOF': (COUPLING_TYPE.CONDENSABLE_DOF, None), 'INTERFACE_DOF': (COUPLING_TYPE.INTERFACE_DOF, None), 'NONWIREBASKET_DOF': (COUPLING_TYPE.NONWIREBASKET_DOF, None), 'WIREBASKET_DOF': (COUPLING_TYPE.WIREBASKET_DOF, None), 'EXTERNAL_DOF': (COUPLING_TYPE.EXTERNAL_DOF, None), 'VISIBLE_DOF': (COUPLING_TYPE.VISIBLE_DOF, None), 'ANY_DOF': (COUPLING_TYPE.ANY_DOF, None)}
    __members__: dict # value = {'UNUSED_DOF': COUPLING_TYPE.UNUSED_DOF, 'HIDDEN_DOF': COUPLING_TYPE.HIDDEN_DOF, 'LOCAL_DOF': COUPLING_TYPE.LOCAL_DOF, 'CONDENSABLE_DOF': COUPLING_TYPE.CONDENSABLE_DOF, 'INTERFACE_DOF': COUPLING_TYPE.INTERFACE_DOF, 'NONWIREBASKET_DOF': COUPLING_TYPE.NONWIREBASKET_DOF, 'WIREBASKET_DOF': COUPLING_TYPE.WIREBASKET_DOF, 'EXTERNAL_DOF': COUPLING_TYPE.EXTERNAL_DOF, 'VISIBLE_DOF': COUPLING_TYPE.VISIBLE_DOF, 'ANY_DOF': COUPLING_TYPE.ANY_DOF}
    pass
class GridFunction(ngsolve.fem.CoefficientFunction, NGS_Object):
    """
    a field approximated in some finite element space
     Keyword arguments can be:
    multidim: 
     Multidimensional GridFunction
    nested: bool = False
     Generates prolongation matrices for each mesh level and prolongates
     the solution onto the finer grid after a refinement.
    """
    def AddMultiDimComponent(self, arg0: ngsolve.la.BaseVector) -> None: ...
    def CF(self, diffop: ngsolve.fem.DifferentialOperator) -> ngsolve.fem.CoefficientFunction: 
        """
        Parameters:

        diffop : ngsolve.fem.DifferentialOperator
          input differential operator
        """
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> ngsolve.fem.CoefficientFunction: 
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
    def Deriv(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns the canonical derivative of the space behind the GridFunction if possible.
        """
    def Derive(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = None) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> ngsolve.fem.CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def InnerProduct(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Load(self, filename: str, parallel: bool = False) -> None: 
        """
               
        Loads a gridfunction from a file.

        Parameters:

        filename : string
          input file name

        parallel : bool
          input parallel
        """
    def Norm(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Operator(self, name: str, VOL_or_BND: VorB = VorB.VOL) -> GridFunctionCoefficientFunction: 
        """
        Get access to an operator depending on the FESpace.

        Parameters:

        name : string
          input name of the requested operator

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND, ...
        """
    def Operators(self) -> list: 
        """
        returns list of available differential operators
        """
    def Other(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Evaluate on other element, as needed for DG jumps
        """
    def Save(self, filename: str, parallel: bool = False) -> None: 
        """
        Saves the gridfunction into a file.

        Parameters:

        filename : string
          input file name

        parallel : bool
          input parallel
        """
    def Set(self, coefficient: ngsolve.fem.CoefficientFunction, VOL_or_BND: VorB = VorB.VOL, definedon: object = <ngsolve.ngstd.DummyArgument>, dual: bool = False, use_simd: bool = True) -> None: 
        """
        Set values

        Parameters:

        coefficient : ngsolve.fem.CoefficientFunction
          input CF to set

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND, ...

        definedon : object
          input definedon region

        dual : bool
          If set to true dual shapes are used, otherwise local L2-projection is used.
          Default is False.

        use_simd : bool
          If set to false does not use SIMD (for debugging).
        """
    def Trace(self) -> GridFunctionCoefficientFunction: 
        """
        take canonical boundary trace. This function is optional, added for consistency with proxies
        """
    def Update(self) -> None: 
        """
        update vector size to finite element space dimension after mesh refinement
        """
    @overload
    def __add__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __call__(self, *args, **kwargs) -> object: ...
    @overload
    def __call__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0, VOL_or_BND: VorB = VorB.VOL) -> object: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    @overload
    def __getitem__(self, components: tuple) -> ngsolve.fem.CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, comp: int) -> ngsolve.fem.CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, space: FESpace, name: str = 'gfu', autoupdate: bool = False, **kwargs) -> None: 
        """
        creates a gridfunction in finite element space
        """
    @overload
    def __mul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    def __neg__(self) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, exponent: int) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: ngsolve.fem.CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    def __rsub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __sub__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        list of gridfunctions for compound gridfunction

        :type: tuple
        """
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def derivname(self) -> str:
        """
        Name of canonical derivative of the space behind the GridFunction.

        :type: str
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
    def imag(self) -> ngsolve.fem.CoefficientFunction:
        """
        imaginary part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def name(self) -> str:
        """
        Name of the Gridfunction

        :type: str
        """
    @property
    def real(self) -> ngsolve.fem.CoefficientFunction:
        """
        real part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def space(self) -> FESpace:
        """
        the finite element space

        :type: FESpace
        """
    @property
    def trans(self) -> ngsolve.fem.CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def vec(self) -> ngsolve.la.BaseVector:
        """
        coefficient vector

        :type: ngsolve.la.BaseVector
        """
    @property
    def vecs(self) -> list:
        """
        list of coefficient vectors for multi-dim gridfunction

        :type: list
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE41B288>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.GridFunction' objects>, '__doc__': 'a field approximated in some finite element space\n Keyword arguments can be:\nmultidim: \n Multidimensional GridFunction\nnested: bool = False\n Generates prolongation matrices for each mesh level and prolongates\n the solution onto the finer grid after a refinement.\n', '__module__': 'ngsolve.comp', '__flags_doc__': <staticmethod object at 0x0000026DEE41C0C8>, '__getstate__': <instancemethod __getstate__ at 0x0000026DEE41B318>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE41B378>, '__str__': <instancemethod __str__ at 0x0000026DEE41B3D8>, 'space': <property object at 0x0000026DEE41AC28>, 'Update': <instancemethod Update at 0x0000026DEE41B468>, 'Save': <instancemethod Save at 0x0000026DEE41B4C8>, 'Load': <instancemethod Load at 0x0000026DEE41B528>, 'Set': <instancemethod Set at 0x0000026DEE41B588>, 'name': <property object at 0x0000026DEE41AE58>, 'components': <property object at 0x0000026DEE41AEF8>, 'vec': <property object at 0x0000026DEE41D048>, 'vecs': <property object at 0x0000026DEE41D0E8>, 'Deriv': <instancemethod Deriv at 0x0000026DEE41B6A8>, 'Trace': <instancemethod Trace at 0x0000026DEE41B708>, 'Operators': <instancemethod Operators at 0x0000026DEE41B768>, 'Operator': <instancemethod Operator at 0x0000026DEE41B7C8>, 'derivname': <property object at 0x0000026DEE41D2C8>, '__call__': <instancemethod __call__ at 0x0000026DEE41B888>, 'CF': <instancemethod CF at 0x0000026DEE41B8B8>, 'AddMultiDimComponent': <instancemethod AddMultiDimComponent at 0x0000026DEE41B918>})
    pass
class FESpace(NGS_Object):
    """
    Finite Element Space

    Provides the functionality for finite element calculations.

    Some available FESpaces are:

    H1, HCurl, HDiv, L2, FacetFESpace, HDivDiv

    2 __init__ overloads:
      1) To create a registered FESpace
      2) To create a compound FESpace from multiple created FESpaces

    1)

    Parameters:

    type : string
      Type of the finite element space. This parameter is automatically
      set if the space is constructed with a generator function.

    mesh : ngsolve.Mesh
      Mesh on which the finite element space is defined on.

    kwargs : kwargs
      For a description of the possible kwargs have a look a bit further down.

    2)

    Parameters:

    spaces : list of ngsolve.FESpace
      List of the spaces for the compound finite element space

    kwargs : kwargs
      For a description of the possible kwargs have a look a bit further down.


     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __init__(self, type: str, mesh: Mesh, **kwargs) -> None: 
        """
        construct product space (compound-space) from list of component spaces

        allowed types are: 'h1ho', 'l2ho', 'hcurlho', 'hdivho' etc.
        """
    @overload
    def __init__(self, spaces: list, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE40AE88>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.FESpace' objects>, '__doc__': "Finite Element Space\n\nProvides the functionality for finite element calculations.\n\nSome available FESpaces are:\n\nH1, HCurl, HDiv, L2, FacetFESpace, HDivDiv\n\n2 __init__ overloads:\n  1) To create a registered FESpace\n  2) To create a compound FESpace from multiple created FESpaces\n\n1)\n\nParameters:\n\ntype : string\n  Type of the finite element space. This parameter is automatically\n  set if the space is constructed with a generator function.\n\nmesh : ngsolve.Mesh\n  Mesh on which the finite element space is defined on.\n\nkwargs : kwargs\n  For a description of the possible kwargs have a look a bit further down.\n\n2)\n\nParameters:\n\nspaces : list of ngsolve.FESpace\n  List of the spaces for the compound finite element space\n\nkwargs : kwargs\n  For a description of the possible kwargs have a look a bit further down.\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__flags_doc__': <staticmethod object at 0x0000026DEE40E0C8>, '__special_treated_flags__': <staticmethod object at 0x0000026DEE409F08>, '__getstate__': <instancemethod __getstate__ at 0x0000026DEE40AF18>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE40AF78>, 'Update': <instancemethod Update at 0x0000026DEE40AFD8>, 'UpdateDofTables': <instancemethod UpdateDofTables at 0x0000026DEE410078>, 'FinalizeUpdate': <instancemethod FinalizeUpdate at 0x0000026DEE4100D8>, 'HideAllDofs': <instancemethod HideAllDofs at 0x0000026DEE410138>, 'ndof': <property object at 0x0000026DEE40D728>, 'ndofglobal': <property object at 0x0000026DEE40D7C8>, 'dim': <property object at 0x0000026DEE40D8B8>, '__str__': <instancemethod __str__ at 0x0000026DEE410228>, '__timing__': <instancemethod __timing__ at 0x0000026DEE410288>, 'lospace': <property object at 0x0000026DEE40D9F8>, 'loembedding': <property object at 0x0000026DEE40DA98>, 'mesh': <property object at 0x0000026DEE40DB38>, 'globalorder': <property object at 0x0000026DEE40DC28>, 'type': <property object at 0x0000026DEE40DD18>, 'is_complex': <property object at 0x0000026DEE40DDB8>, 'SetDefinedOn': <instancemethod SetDefinedOn at 0x0000026DEE410408>, 'SetOrder': <instancemethod SetOrder at 0x0000026DEE410498>, 'GetOrder': <instancemethod GetOrder at 0x0000026DEE4104C8>, 'Elements': <instancemethod Elements at 0x0000026DEE410528>, 'GetDofNrs': <instancemethod GetDofNrs at 0x0000026DEE4105B8>, 'GetDofs': <instancemethod GetDofs at 0x0000026DEE4105E8>, 'CouplingType': <instancemethod CouplingType at 0x0000026DEE410648>, 'SetCouplingType': <instancemethod SetCouplingType at 0x0000026DEE4106D8>, 'couplingtype': <property object at 0x0000026DEE40F138>, 'GetFE': <instancemethod GetFE at 0x0000026DEE410738>, 'FreeDofs': <instancemethod FreeDofs at 0x0000026DEE410798>, 'ParallelDofs': <instancemethod ParallelDofs at 0x0000026DEE4107F8>, 'Prolongation': <instancemethod Prolongation at 0x0000026DEE410858>, 'Range': <instancemethod Range at 0x0000026DEE4108B8>, 'components': <property object at 0x0000026DEE40F368>, 'TrialFunction': <instancemethod TrialFunction at 0x0000026DEE410948>, 'TestFunction': <instancemethod TestFunction at 0x0000026DEE4109A8>, 'TnT': <instancemethod TnT at 0x0000026DEE410A08>, 'InvM': <instancemethod InvM at 0x0000026DEE410A68>, 'Mass': <instancemethod Mass at 0x0000026DEE410AC8>, 'SolveM': <instancemethod SolveM at 0x0000026DEE410B28>, 'ApplyM': <instancemethod ApplyM at 0x0000026DEE410B88>, 'TraceOperator': <instancemethod TraceOperator at 0x0000026DEE410BE8>, 'ConvertL2Operator': <instancemethod ConvertL2Operator at 0x0000026DEE410C48>, 'GetTrace': <instancemethod GetTrace at 0x0000026DEE410CA8>, 'GetTraceTrans': <instancemethod GetTraceTrans at 0x0000026DEE410D08>, '__eq__': <instancemethod __eq__ at 0x0000026DEE410D68>})
    pass
class Compress(FESpace, NGS_Object):
    """
    Wrapper Finite Element Spaces.
    The compressed fespace is a wrapper around a standard fespace which removes
    certain dofs (e.g. UNUSED_DOFs).

    Parameters:

    fespace : ngsolve.comp.FESpace
        finite element space

    active_dofs : BitArray or None
        don't use the COUPLING_TYPEs of dofs to compress the FESpace, 
        but use a BitArray directly to compress the FESpace

     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    def GetBaseSpace(self) -> FESpace: ...
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    def SetActiveDofs(self, dofs: pyngcore.BitArray) -> None: ...
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, fespace: FESpace, active_dofs: object = <ngsolve.ngstd.DummyArgument>) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416F18>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.Compress' objects>, '__doc__': "Wrapper Finite Element Spaces.\nThe compressed fespace is a wrapper around a standard fespace which removes\ncertain dofs (e.g. UNUSED_DOFs).\n\nParameters:\n\nfespace : ngsolve.comp.FESpace\n    finite element space\n\nactive_dofs : BitArray or None\n    don't use the COUPLING_TYPEs of dofs to compress the FESpace, \n    but use a BitArray directly to compress the FESpace\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', 'SetActiveDofs': <instancemethod SetActiveDofs at 0x0000026DEE416F78>, 'GetBaseSpace': <instancemethod GetBaseSpace at 0x0000026DEE416FD8>, '__getstate__': <instancemethod __getstate__ at 0x0000026DEE41B078>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE41B0D8>})
    pass
class ContactBoundary():
    def AddEnergy(self, arg0: ngsolve.fem.CoefficientFunction) -> None: ...
    def AddIntegrator(self, arg0: ngsolve.fem.CoefficientFunction) -> None: ...
    def Update(self, gf: GridFunction, bf: BilinearForm = None, intorder: int = 4, maxdist: float = 0.0) -> None: 
        """
        Update searchtree for gap function.
        If bf is given add specialelements corresponding to
        integrationrules of order 'intorder' on each master
        element to BilinearForm bf.
        `maxdist` is the maximum distance where this function is accurate.
        If `maxdist` == 0. then 2*meshsize is used.
        """
    def __init__(self, fes: FESpace, master: Region, minion: Region, draw_pairs: bool = False) -> None: 
        """
        Class for managing contact interfaces.
        The created object must be kept alive in python as long as
        operations of it are used!
        """
    @property
    def gap(self) -> ngsolve.fem.CoefficientFunction:
        """
        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def normal(self) -> ngsolve.fem.CoefficientFunction:
        """
        :type: ngsolve.fem.CoefficientFunction
        """
    pass
class DifferentialSymbol():
    def __call__(self, definedon: Optional[Union[Region, str]] = None, element_boundary: bool = False, element_vb: VorB = VorB.VOL, skeleton: bool = False, bonus_intorder: int = 0, intrules: Dict[ngsolve.fem.ET, ngsolve.fem.IntegrationRule] = {}, deformation: GridFunction = None, definedonelements: pyngcore.BitArray = None) -> DifferentialSymbol: ...
    def __init__(self, arg0: VorB) -> None: ...
    pass
class Discontinuous(FESpace, NGS_Object):
    """
    Discontinuous Finite Element Spaces.
    FESpace that splits up all dofs that are shared by several (volume or surface) elements. Every element gets a single copy of that dof. Basis functions become element-local.

    Parameters:

    fespace : ngsolve.comp.FESpace
        finite element space

    BND : boolean or None
        separate across surface elements instead of volume elements (for surface FESpaces)

     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, fespace: FESpace, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416E28>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.Discontinuous' objects>, '__doc__': "Discontinuous Finite Element Spaces.\nFESpace that splits up all dofs that are shared by several (volume or surface) elements. Every element gets a single copy of that dof. Basis functions become element-local.\n\nParameters:\n\nfespace : ngsolve.comp.FESpace\n    finite element space\n\nBND : boolean or None\n    separate across surface elements instead of volume elements (for surface FESpaces)\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp'})
    pass
class ProxyFunction(ngsolve.fem.CoefficientFunction):
    """
    Either FESpace.TrialFunction or FESpace.TestFunction. Is a
    placeholder coefficient function for Symbolic Integrators. The
    integrators will replace it with the basis functions of the finite element space
    when building the system matrices.
    """
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> ngsolve.fem.CoefficientFunction: 
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
    def Deriv(self) -> ProxyFunction: 
        """
        take canonical derivative (grad, curl, div)
        """
    def Derive(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = None) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> ngsolve.fem.CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def InnerProduct(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Norm(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Operator(self, name: str) -> ProxyFunction: 
        """
        Use an additional operator of the finite element space
        """
    def Operators(self) -> list: 
        """
        returns list of available differential operators
        """
    def Other(self, bnd: object = <ngsolve.ngstd.DummyArgument>) -> ProxyFunction: 
        """
        take value from neighbour element (DG)
        """
    def Trace(self) -> ProxyFunction: 
        """
        take canonical boundary trace
        """
    @overload
    def __add__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __call__(self, mip: ngsolve.fem.BaseMappedIntegrationPoint) -> object: 
        """
        evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)
        """
    @overload
    def __call__(self, arg0: numpy.ndarray[ngsolve.fem.MeshPoint]) -> array: ...
    def __diffop__(self) -> ngsolve.fem.DifferentialOperator: ...
    @overload
    def __getitem__(self, components: tuple) -> ngsolve.fem.CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, comp: int) -> ngsolve.fem.CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, arg0: ProxyFunction) -> None: ...
    @overload
    def __mul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    def __neg__(self) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, exponent: int) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: ngsolve.fem.CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    def __rsub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __sub__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def derivname(self) -> str:
        """
        name of the canonical derivative

        :type: str
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
    def imag(self) -> ngsolve.fem.CoefficientFunction:
        """
        imaginary part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def real(self) -> ngsolve.fem.CoefficientFunction:
        """
        real part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def trans(self) -> ngsolve.fem.CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: ngsolve.fem.CoefficientFunction
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE40A828>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.ProxyFunction' objects>, '__doc__': '\nEither FESpace.TrialFunction or FESpace.TestFunction. Is a\nplaceholder coefficient function for Symbolic Integrators. The\nintegrators will replace it with the basis functions of the finite element space\nwhen building the system matrices.\n\n', '__module__': 'ngsolve.comp', 'Deriv': <instancemethod Deriv at 0x0000026DEE40A888>, 'Trace': <instancemethod Trace at 0x0000026DEE40A8E8>, 'Other': <instancemethod Other at 0x0000026DEE40A948>, 'derivname': <property object at 0x0000026DEE40B8B8>, 'Operator': <instancemethod Operator at 0x0000026DEE40A9D8>, 'Operators': <instancemethod Operators at 0x0000026DEE40AA38>, '__diffop__': <instancemethod __diffop__ at 0x0000026DEE40AA98>})
    pass
class ElementId():
    """
    An element identifier containing element number and Volume/Boundary flag

    3 __init__ overloads:

    1)

    Parameters:

    vb : ngsolve.comp.VorB
      input Volume or Boundary (VOL, BND, BBND, BBBND)

    nr : int
      input element number


    2)

    Parameters:

    nr : int
      input element number


    3)

    Parameters:

    el : ngcomp::Ngs_Element
      input Ngs element
    """
    def VB(self) -> VorB: 
        """
        VorB of element
        """
    def __eq__(self, arg0: ElementId) -> bool: ...
    def __hash__(self) -> int: ...
    @overload
    def __init__(self, vb: VorB, nr: int) -> None: ...
    @overload
    def __init__(self, el: ngcomp::Ngs_Element) -> None: ...
    @overload
    def __init__(self, nr: int) -> None: ...
    def __ne__(self, arg0: ElementId) -> bool: ...
    def __str__(self) -> str: ...
    @property
    def nr(self) -> int:
        """
        the element number

        :type: int
        """
    @property
    def valid(self) -> bool:
        """
        is element valid

        :type: bool
        """
    pass
class ElementRange(ngsolve.ngstd.IntRange):
    def __init__(self, mesh: Mesh, vb: VorB, range: ngsolve.ngstd.IntRange) -> None: ...
    def __iter__(self) -> iterator: ...
    def __str__(self) -> str: ...
    @property
    def start(self) -> int:
        """
        :type: int
        """
    @property
    def step(self) -> int:
        """
        :type: int
        """
    @property
    def stop(self) -> int:
        """
        :type: int
        """
    pass
class CompoundFESpace(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <slot wrapper '__init__' of 'ngsolve.comp.CompoundFESpace' objects>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.CompoundFESpace' objects>, '__doc__': "\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE410DC8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE410E28>})
    pass
class Ngs_Element():
    def VB(self) -> VorB: 
        """
        VorB of element
        """
    @property
    def edges(self) -> tuple:
        """
        tuple of global edge numbers

        :type: tuple
        """
    @property
    def elementnode(self) -> NodeId:
        """
        inner node, i.e. cell, face or edge node for 3D/2D/1D

        :type: NodeId
        """
    @property
    def faces(self) -> tuple:
        """
        tuple of global face numbers

        :type: tuple
        """
    @property
    def facets(self) -> tuple:
        """
        tuple of global face, edge or vertex numbers

        :type: tuple
        """
    @property
    def index(self) -> int:
        """
        material or boundary condition index

        :type: int
        """
    @property
    def mat(self) -> str:
        """
        material or boundary condition label

        :type: str
        """
    @property
    def nr(self) -> int:
        """
        the element number

        :type: int
        """
    @property
    def type(self) -> ngsolve.fem.ET:
        """
        geometric shape of element

        :type: ngsolve.fem.ET
        """
    @property
    def valid(self) -> bool:
        """
        is element valid

        :type: bool
        """
    @property
    def vertices(self) -> tuple:
        """
        tuple of global vertex numbers

        :type: tuple
        """
    pass
class FESpaceElementRange(ngsolve.ngstd.IntRange):
    def __iter__(self) -> iterator: ...
    def __str__(self) -> str: ...
    @property
    def start(self) -> int:
        """
        :type: int
        """
    @property
    def step(self) -> int:
        """
        :type: int
        """
    @property
    def stop(self) -> int:
        """
        :type: int
        """
    pass
class FacetFESpace(FESpace, NGS_Object):
    """
    A finite element space living on facets.

    The FacetFESpace provides polynomials on facets, i.e. faces in 3D,
    edges in 2D, and vertices in 1D. The functions are discontinuous from facet to facet.

    Typecal usecases for the FacetFESpace are hybrid mixed and hybrid DG methods.

    The function is only defined on the mesh skeleton. Evaluation inside the element throws
    an exception. Thus, functions from the FacetFESpace can be used only within element_boundary 
    or skeleton expressions. 

    Functions have meaningful boundary-values, which are obtained using the Trace-operator.
    (the trace operator might become redundant in future).

    (coming soon) The FacetFESpace provides variable order, which can be set for FACET-nodes. Alternatively,
    one can use FACE, EDGE, or VERTEX nodes for 3D, 2D, or 1D meshes, respectively.

    The basis is L2-orthogonal on the facets. The highest order basis functions can be duplicated
    for the two neighbouring elements. This allows a simple implementation of the Lehrenfeld-Schoeberl
    'projected jumps' HDG method.

     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    highest_order_dc: bool = False
      Splits highest order facet functions into two which are associated with
      the corresponding neighbors and are local dofs on the corresponding element
      (used to realize projected jumps)
    hide_highest_order_dc: bool = False
      if highest_order_dc is used this flag marks the corresponding local dofs
      as hidden dofs (reduces number of non-zero entries in a matrix). These dofs
      can also be compressed.
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416228>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.FacetFESpace' objects>, '__doc__': "A finite element space living on facets.\n\nThe FacetFESpace provides polynomials on facets, i.e. faces in 3D,\nedges in 2D, and vertices in 1D. The functions are discontinuous from facet to facet.\n\nTypecal usecases for the FacetFESpace are hybrid mixed and hybrid DG methods.\n\nThe function is only defined on the mesh skeleton. Evaluation inside the element throws\nan exception. Thus, functions from the FacetFESpace can be used only within element_boundary \nor skeleton expressions. \n\nFunctions have meaningful boundary-values, which are obtained using the Trace-operator.\n(the trace operator might become redundant in future).\n\n(coming soon) The FacetFESpace provides variable order, which can be set for FACET-nodes. Alternatively,\none can use FACE, EDGE, or VERTEX nodes for 3D, 2D, or 1D meshes, respectively.\n\nThe basis is L2-orthogonal on the facets. The highest order basis functions can be duplicated\nfor the two neighbouring elements. This allows a simple implementation of the Lehrenfeld-Schoeberl\n'projected jumps' HDG method.\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\nhighest_order_dc: bool = False\n  Splits highest order facet functions into two which are associated with\n  the corresponding neighbors and are local dofs on the corresponding element\n  (used to realize projected jumps)\nhide_highest_order_dc: bool = False\n  if highest_order_dc is used this flag marks the corresponding local dofs\n  as hidden dofs (reduces number of non-zero entries in a matrix). These dofs\n  can also be compressed.\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE416288>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE4162E8>, '__flags_doc__': <staticmethod object at 0x0000026DEE414708>})
    pass
class FacetSurface(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416378>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.FacetSurface' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE4163D8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE416438>, '__flags_doc__': <staticmethod object at 0x0000026DEE414808>})
    pass
class Array_enum ngcomp::COUPLING_TYPE_S(FlatArray_enum ngcomp::COUPLING_TYPE_S):
    def NumPy(self) -> object: ...
    def __getitem__(self, arg0: int) -> COUPLING_TYPE: ...
    @overload
    def __init__(self, vec: List[COUPLING_TYPE]) -> None: 
        """
        Makes array of given length

        Makes array with given list of elements
        """
    @overload
    def __init__(self, n: int) -> None: ...
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: ...
    @overload
    def __setitem__(self, arg0: int, arg1: COUPLING_TYPE) -> COUPLING_TYPE: ...
    @overload
    def __setitem__(self, arg0: slice, arg1: COUPLING_TYPE) -> None: ...
    def __str__(self) -> str: ...
    pass
class GlobalVariables():
    @property
    def msg_level(self) -> int:
        """
        message level

        :type: int
        """
    @msg_level.setter
    def msg_level(self, arg1: int) -> None:
        """
        message level
        """
    @property
    def numthreads(self) -> int:
        """
        :type: int
        """
    @numthreads.setter
    def numthreads(self, arg1: int) -> None:
        pass
    @property
    def pajetrace(self) -> str:
        """
        :type: str
        """
    @pajetrace.setter
    def pajetrace(self, arg1: int) -> None:
        pass
    @property
    def testout(self) -> str:
        """
        testout file

        :type: str
        """
    @testout.setter
    def testout(self, arg1: str) -> None:
        """
        testout file
        """
    pass
class ComponentGridFunction(GridFunction, ngsolve.fem.CoefficientFunction, NGS_Object):
    """
     Keyword arguments can be:
    multidim: 
     Multidimensional GridFunction
    nested: bool = False
     Generates prolongation matrices for each mesh level and prolongates
     the solution onto the finer grid after a refinement.
    """
    def AddMultiDimComponent(self, arg0: ngsolve.la.BaseVector) -> None: ...
    def CF(self, diffop: ngsolve.fem.DifferentialOperator) -> ngsolve.fem.CoefficientFunction: 
        """
        Parameters:

        diffop : ngsolve.fem.DifferentialOperator
          input differential operator
        """
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> ngsolve.fem.CoefficientFunction: 
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
    def Deriv(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns the canonical derivative of the space behind the GridFunction if possible.
        """
    def Derive(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = None) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> ngsolve.fem.CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def InnerProduct(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Load(self, filename: str, parallel: bool = False) -> None: 
        """
               
        Loads a gridfunction from a file.

        Parameters:

        filename : string
          input file name

        parallel : bool
          input parallel
        """
    def Norm(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Operator(self, name: str, VOL_or_BND: VorB = VorB.VOL) -> GridFunctionCoefficientFunction: 
        """
        Get access to an operator depending on the FESpace.

        Parameters:

        name : string
          input name of the requested operator

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND, ...
        """
    def Operators(self) -> list: 
        """
        returns list of available differential operators
        """
    def Other(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Evaluate on other element, as needed for DG jumps
        """
    def Save(self, filename: str, parallel: bool = False) -> None: 
        """
        Saves the gridfunction into a file.

        Parameters:

        filename : string
          input file name

        parallel : bool
          input parallel
        """
    def Set(self, coefficient: ngsolve.fem.CoefficientFunction, VOL_or_BND: VorB = VorB.VOL, definedon: object = <ngsolve.ngstd.DummyArgument>, dual: bool = False, use_simd: bool = True) -> None: 
        """
        Set values

        Parameters:

        coefficient : ngsolve.fem.CoefficientFunction
          input CF to set

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND, ...

        definedon : object
          input definedon region

        dual : bool
          If set to true dual shapes are used, otherwise local L2-projection is used.
          Default is False.

        use_simd : bool
          If set to false does not use SIMD (for debugging).
        """
    def Trace(self) -> GridFunctionCoefficientFunction: 
        """
        take canonical boundary trace. This function is optional, added for consistency with proxies
        """
    def Update(self) -> None: 
        """
        update vector size to finite element space dimension after mesh refinement
        """
    @overload
    def __add__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __call__(self, *args, **kwargs) -> object: ...
    @overload
    def __call__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0, VOL_or_BND: VorB = VorB.VOL) -> object: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    @overload
    def __getitem__(self, components: tuple) -> ngsolve.fem.CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, comp: int) -> ngsolve.fem.CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __mul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    def __neg__(self) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, exponent: int) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: ngsolve.fem.CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    def __rsub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __sub__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        list of gridfunctions for compound gridfunction

        :type: tuple
        """
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def derivname(self) -> str:
        """
        Name of canonical derivative of the space behind the GridFunction.

        :type: str
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
    def imag(self) -> ngsolve.fem.CoefficientFunction:
        """
        imaginary part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def name(self) -> str:
        """
        Name of the Gridfunction

        :type: str
        """
    @property
    def real(self) -> ngsolve.fem.CoefficientFunction:
        """
        real part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def space(self) -> FESpace:
        """
        the finite element space

        :type: FESpace
        """
    @property
    def trans(self) -> ngsolve.fem.CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def vec(self) -> ngsolve.la.BaseVector:
        """
        coefficient vector

        :type: ngsolve.la.BaseVector
        """
    @property
    def vecs(self) -> list:
        """
        list of coefficient vectors for multi-dim gridfunction

        :type: list
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <slot wrapper '__init__' of 'ngsolve.comp.ComponentGridFunction' objects>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.ComponentGridFunction' objects>, '__doc__': '\n Keyword arguments can be:\nmultidim: \n Multidimensional GridFunction\nnested: bool = False\n Generates prolongation matrices for each mesh level and prolongates\n the solution onto the finer grid after a refinement.\n', '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE41BAF8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE41BB58>})
    pass
class GridFunctionC(GridFunction, ngsolve.fem.CoefficientFunction, NGS_Object):
    """
     Keyword arguments can be:
    multidim: 
     Multidimensional GridFunction
    nested: bool = False
     Generates prolongation matrices for each mesh level and prolongates
     the solution onto the finer grid after a refinement.
    """
    def AddMultiDimComponent(self, arg0: ngsolve.la.BaseVector) -> None: ...
    def CF(self, diffop: ngsolve.fem.DifferentialOperator) -> ngsolve.fem.CoefficientFunction: 
        """
        Parameters:

        diffop : ngsolve.fem.DifferentialOperator
          input differential operator
        """
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> ngsolve.fem.CoefficientFunction: 
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
    def Deriv(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns the canonical derivative of the space behind the GridFunction if possible.
        """
    def Derive(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = None) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> ngsolve.fem.CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def InnerProduct(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Load(self, filename: str, parallel: bool = False) -> None: 
        """
               
        Loads a gridfunction from a file.

        Parameters:

        filename : string
          input file name

        parallel : bool
          input parallel
        """
    def Norm(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Operator(self, name: str, VOL_or_BND: VorB = VorB.VOL) -> GridFunctionCoefficientFunction: 
        """
        Get access to an operator depending on the FESpace.

        Parameters:

        name : string
          input name of the requested operator

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND, ...
        """
    def Operators(self) -> list: 
        """
        returns list of available differential operators
        """
    def Other(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Evaluate on other element, as needed for DG jumps
        """
    def Save(self, filename: str, parallel: bool = False) -> None: 
        """
        Saves the gridfunction into a file.

        Parameters:

        filename : string
          input file name

        parallel : bool
          input parallel
        """
    def Set(self, coefficient: ngsolve.fem.CoefficientFunction, VOL_or_BND: VorB = VorB.VOL, definedon: object = <ngsolve.ngstd.DummyArgument>, dual: bool = False, use_simd: bool = True) -> None: 
        """
        Set values

        Parameters:

        coefficient : ngsolve.fem.CoefficientFunction
          input CF to set

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND, ...

        definedon : object
          input definedon region

        dual : bool
          If set to true dual shapes are used, otherwise local L2-projection is used.
          Default is False.

        use_simd : bool
          If set to false does not use SIMD (for debugging).
        """
    def Trace(self) -> GridFunctionCoefficientFunction: 
        """
        take canonical boundary trace. This function is optional, added for consistency with proxies
        """
    def Update(self) -> None: 
        """
        update vector size to finite element space dimension after mesh refinement
        """
    @overload
    def __add__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __call__(self, *args, **kwargs) -> object: ...
    @overload
    def __call__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0, VOL_or_BND: VorB = VorB.VOL) -> object: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    @overload
    def __getitem__(self, components: tuple) -> ngsolve.fem.CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, comp: int) -> ngsolve.fem.CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __mul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    def __neg__(self) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, exponent: int) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: ngsolve.fem.CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    def __rsub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __sub__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        list of gridfunctions for compound gridfunction

        :type: tuple
        """
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def derivname(self) -> str:
        """
        Name of canonical derivative of the space behind the GridFunction.

        :type: str
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
    def imag(self) -> ngsolve.fem.CoefficientFunction:
        """
        imaginary part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def name(self) -> str:
        """
        Name of the Gridfunction

        :type: str
        """
    @property
    def real(self) -> ngsolve.fem.CoefficientFunction:
        """
        real part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def space(self) -> FESpace:
        """
        the finite element space

        :type: FESpace
        """
    @property
    def trans(self) -> ngsolve.fem.CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def vec(self) -> ngsolve.la.BaseVector:
        """
        coefficient vector

        :type: ngsolve.la.BaseVector
        """
    @property
    def vecs(self) -> list:
        """
        list of coefficient vectors for multi-dim gridfunction

        :type: list
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <slot wrapper '__init__' of 'ngsolve.comp.GridFunctionC' objects>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.GridFunctionC' objects>, '__doc__': '\n Keyword arguments can be:\nmultidim: \n Multidimensional GridFunction\nnested: bool = False\n Generates prolongation matrices for each mesh level and prolongates\n the solution onto the finer grid after a refinement.\n', '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE41BA38>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE41BA98>})
    pass
class GridFunctionCoefficientFunction(ngsolve.fem.CoefficientFunction):
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> ngsolve.fem.CoefficientFunction: 
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
    def Derive(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = None) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> ngsolve.fem.CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def InnerProduct(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Norm(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Other(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Evaluate on other element, as needed for DG jumps
        """
    def Trace(self) -> GridFunctionCoefficientFunction: 
        """
        take canonical boundary trace.
        """
    @overload
    def __add__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __call__(self, mip: ngsolve.fem.BaseMappedIntegrationPoint) -> object: 
        """
        evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)
        """
    @overload
    def __call__(self, arg0: numpy.ndarray[ngsolve.fem.MeshPoint]) -> array: ...
    @overload
    def __getitem__(self, components: tuple) -> ngsolve.fem.CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, comp: int) -> ngsolve.fem.CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __mul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    def __neg__(self) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, exponent: int) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: ngsolve.fem.CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    def __rsub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __sub__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
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
    def imag(self) -> ngsolve.fem.CoefficientFunction:
        """
        imaginary part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def real(self) -> ngsolve.fem.CoefficientFunction:
        """
        real part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def trans(self) -> ngsolve.fem.CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: ngsolve.fem.CoefficientFunction
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <slot wrapper '__init__' of 'ngsolve.comp.GridFunctionCoefficientFunction' objects>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.GridFunctionCoefficientFunction' objects>, '__doc__': None, '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE41B168>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE41B1C8>, 'Trace': <instancemethod Trace at 0x0000026DEE41B228>})
    pass
class GridFunctionD(GridFunction, ngsolve.fem.CoefficientFunction, NGS_Object):
    """
     Keyword arguments can be:
    multidim: 
     Multidimensional GridFunction
    nested: bool = False
     Generates prolongation matrices for each mesh level and prolongates
     the solution onto the finer grid after a refinement.
    """
    def AddMultiDimComponent(self, arg0: ngsolve.la.BaseVector) -> None: ...
    def CF(self, diffop: ngsolve.fem.DifferentialOperator) -> ngsolve.fem.CoefficientFunction: 
        """
        Parameters:

        diffop : ngsolve.fem.DifferentialOperator
          input differential operator
        """
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> ngsolve.fem.CoefficientFunction: 
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
    def Deriv(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns the canonical derivative of the space behind the GridFunction if possible.
        """
    def Derive(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = None) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> ngsolve.fem.CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def InnerProduct(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Load(self, filename: str, parallel: bool = False) -> None: 
        """
               
        Loads a gridfunction from a file.

        Parameters:

        filename : string
          input file name

        parallel : bool
          input parallel
        """
    def Norm(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Operator(self, name: str, VOL_or_BND: VorB = VorB.VOL) -> GridFunctionCoefficientFunction: 
        """
        Get access to an operator depending on the FESpace.

        Parameters:

        name : string
          input name of the requested operator

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND, ...
        """
    def Operators(self) -> list: 
        """
        returns list of available differential operators
        """
    def Other(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Evaluate on other element, as needed for DG jumps
        """
    def Save(self, filename: str, parallel: bool = False) -> None: 
        """
        Saves the gridfunction into a file.

        Parameters:

        filename : string
          input file name

        parallel : bool
          input parallel
        """
    def Set(self, coefficient: ngsolve.fem.CoefficientFunction, VOL_or_BND: VorB = VorB.VOL, definedon: object = <ngsolve.ngstd.DummyArgument>, dual: bool = False, use_simd: bool = True) -> None: 
        """
        Set values

        Parameters:

        coefficient : ngsolve.fem.CoefficientFunction
          input CF to set

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND, ...

        definedon : object
          input definedon region

        dual : bool
          If set to true dual shapes are used, otherwise local L2-projection is used.
          Default is False.

        use_simd : bool
          If set to false does not use SIMD (for debugging).
        """
    def Trace(self) -> GridFunctionCoefficientFunction: 
        """
        take canonical boundary trace. This function is optional, added for consistency with proxies
        """
    def Update(self) -> None: 
        """
        update vector size to finite element space dimension after mesh refinement
        """
    @overload
    def __add__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __call__(self, *args, **kwargs) -> object: ...
    @overload
    def __call__(self, x: float = 0.0, y: float = 0.0, z: float = 0.0, VOL_or_BND: VorB = VorB.VOL) -> object: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    @overload
    def __getitem__(self, components: tuple) -> ngsolve.fem.CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, comp: int) -> ngsolve.fem.CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __mul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    def __neg__(self) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, exponent: int) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: ngsolve.fem.CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    def __rsub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __sub__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        list of gridfunctions for compound gridfunction

        :type: tuple
        """
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def derivname(self) -> str:
        """
        Name of canonical derivative of the space behind the GridFunction.

        :type: str
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
    def imag(self) -> ngsolve.fem.CoefficientFunction:
        """
        imaginary part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def name(self) -> str:
        """
        Name of the Gridfunction

        :type: str
        """
    @property
    def real(self) -> ngsolve.fem.CoefficientFunction:
        """
        real part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def space(self) -> FESpace:
        """
        the finite element space

        :type: FESpace
        """
    @property
    def trans(self) -> ngsolve.fem.CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def vec(self) -> ngsolve.la.BaseVector:
        """
        coefficient vector

        :type: ngsolve.la.BaseVector
        """
    @property
    def vecs(self) -> list:
        """
        list of coefficient vectors for multi-dim gridfunction

        :type: list
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <slot wrapper '__init__' of 'ngsolve.comp.GridFunctionD' objects>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.GridFunctionD' objects>, '__doc__': '\n Keyword arguments can be:\nmultidim: \n Multidimensional GridFunction\nnested: bool = False\n Generates prolongation matrices for each mesh level and prolongates\n the solution onto the finer grid after a refinement.\n', '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE41B978>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE41B9D8>})
    pass
class H1(FESpace, NGS_Object):
    """
    An H1-conforming finite element space.

    The H1 finite element space consists of continuous and
    elemenet-wise polynomial functions. It uses a hierarchical (=modal)
    basis built from integrated Legendre polynomials on tensor-product elements,
    and Jaboci polynomials on simplicial elements. 

    Boundary values are well defined. The function can be used directly on the
    boundary, using the trace operator is optional.

    The H1 space supports variable order, which can be set individually for edges, 
    faces and cells. 

    Internal degrees of freedom are declared as local dofs and are eliminated 
    if static condensation is on.

    The wirebasket consists of all vertex dofs. Optionally, one can include the 
    first (the quadratic bubble) edge basis function, or all edge basis functions
    into the wirebasket.

     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    wb_withedges: bool = true(3D) / false(2D)
      use lowest-order edge dofs for BDDC wirebasket
    wb_fulledges: bool = false
      use all edge dofs for BDDC wirebasket
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412228>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.H1' objects>, '__doc__': "An H1-conforming finite element space.\n\nThe H1 finite element space consists of continuous and\nelemenet-wise polynomial functions. It uses a hierarchical (=modal)\nbasis built from integrated Legendre polynomials on tensor-product elements,\nand Jaboci polynomials on simplicial elements. \n\nBoundary values are well defined. The function can be used directly on the\nboundary, using the trace operator is optional.\n\nThe H1 space supports variable order, which can be set individually for edges, \nfaces and cells. \n\nInternal degrees of freedom are declared as local dofs and are eliminated \nif static condensation is on.\n\nThe wirebasket consists of all vertex dofs. Optionally, one can include the \nfirst (the quadratic bubble) edge basis function, or all edge basis functions\ninto the wirebasket.\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\nwb_withedges: bool = true(3D) / false(2D)\n  use lowest-order edge dofs for BDDC wirebasket\nwb_fulledges: bool = false\n  use all edge dofs for BDDC wirebasket\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE412288>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE4122E8>, '__flags_doc__': <staticmethod object at 0x0000026DEE40E988>})
    pass
class HCurl(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    nograds: bool = False
      Remove higher order gradients of H1 basis functions from HCurl FESpace
    type1: bool = False
      Use type 1 Nedelec elements
    discontinuous: bool = False
      Create discontinuous HCurl space
    gradientdomains: List[int] = None
      Remove high order gradients from domains where the value is 0.
      This list can be generated for example like this:
      graddoms = [1 if mat == 'iron' else 0 for mat in mesh.GetMaterials()]
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def CreateGradient(self) -> tuple: ...
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE410E88>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.HCurl' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\nnograds: bool = False\n  Remove higher order gradients of H1 basis functions from HCurl FESpace\ntype1: bool = False\n  Use type 1 Nedelec elements\ndiscontinuous: bool = False\n  Create discontinuous HCurl space\ngradientdomains: List[int] = None\n  Remove high order gradients from domains where the value is 0.\n  This list can be generated for example like this:\n  graddoms = [1 if mat == 'iron' else 0 for mat in mesh.GetMaterials()]\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE410EE8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE410F48>, '__flags_doc__': <staticmethod object at 0x0000026DEE40E908>, 'CreateGradient': <instancemethod CreateGradient at 0x0000026DEE410FD8>})
    pass
class HCurlCurl(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    discontinuous: bool = False
      Create discontinuous HCurlCurl space
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412CA8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.HCurlCurl' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\ndiscontinuous: bool = False\n  Create discontinuous HCurlCurl space\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE412D08>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE412D68>, '__flags_doc__': <staticmethod object at 0x0000026DEE414408>})
    pass
class HCurlDiv(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    discontinuous: bool = false
      Create discontinuous HCurlDiv space
    ordertrace: int = -1
      Set order of trace bubbles
    orderinner: int = -1
      Set order of inner nt-bubbles
    GGbubbles: bool = false
      Add GG-bubbles for weak-symmetric formulation
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412B58>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.HCurlDiv' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\ndiscontinuous: bool = false\n  Create discontinuous HCurlDiv space\nordertrace: int = -1\n  Set order of trace bubbles\norderinner: int = -1\n  Set order of inner nt-bubbles\nGGbubbles: bool = false\n  Add GG-bubbles for weak-symmetric formulation\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE412BB8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE412C18>, '__flags_doc__': <staticmethod object at 0x0000026DEE414308>})
    pass
class HDiv(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    RT: bool = False
      RT elements for simplicial elements: P^k subset RT_k subset P^{k+1}
    discontinuous: bool = False
      Create discontinuous HDiv space
    hodivfree: bool = False
      Remove high order element bubbles with non zero divergence
    highest_order_dc: bool = False
      Activates relaxed H(div)-conformity. Allows normal discontinuity of highest order facet basis functions
    hide_all_dofs: bool = False
      Set all used dofs to HIDDEN_DOFs
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def Average(self, vector: ngsolve.la.BaseVector) -> None: ...
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412078>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.HDiv' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\nRT: bool = False\n  RT elements for simplicial elements: P^k subset RT_k subset P^{k+1}\ndiscontinuous: bool = False\n  Create discontinuous HDiv space\nhodivfree: bool = False\n  Remove high order element bubbles with non zero divergence\nhighest_order_dc: bool = False\n  Activates relaxed H(div)-conformity. Allows normal discontinuity of highest order facet basis functions\nhide_all_dofs: bool = False\n  Set all used dofs to HIDDEN_DOFs\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE4120D8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE412138>, '__flags_doc__': <staticmethod object at 0x0000026DEE40E948>, 'Average': <instancemethod Average at 0x0000026DEE4121C8>})
    pass
class HDivDiv(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    discontinuous: bool = False
      Create discontinuous HDivDiv space
    plus: bool = False
      Add additional internal element bubble
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412A08>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.HDivDiv' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\ndiscontinuous: bool = False\n  Create discontinuous HDivDiv space\nplus: bool = False\n  Add additional internal element bubble\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE412A68>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE412AC8>, '__flags_doc__': <staticmethod object at 0x0000026DEE414208>})
    pass
class HDivDivSurface(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    discontinuous: bool = False
      Create discontinuous HDivDiv space
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412DF8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.HDivDivSurface' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\ndiscontinuous: bool = False\n  Create discontinuous HDivDiv space\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE412E58>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE412EB8>, '__flags_doc__': <staticmethod object at 0x0000026DEE4144C8>})
    pass
class HDivSurface(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    discontinuous: bool = False
      Create discontinuous HDivSurface space
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def Average(self, vector: ngsolve.la.BaseVector) -> None: ...
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE4164C8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.HDivSurface' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\ndiscontinuous: bool = False\n  Create discontinuous HDivSurface space\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE416528>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE416588>, '__flags_doc__': <staticmethod object at 0x0000026DEE414908>, 'Average': <instancemethod Average at 0x0000026DEE416618>})
    pass
class Integral():
    @property
    def coef(self) -> ngsolve.fem.CoefficientFunction:
        """
        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def symbol(self) -> DifferentialSymbol:
        """
        :type: DifferentialSymbol
        """
    pass
class InterpolateProxy(ProxyFunction, ngsolve.fem.CoefficientFunction):
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> ngsolve.fem.CoefficientFunction: 
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
    def Deriv(self) -> ProxyFunction: 
        """
        take canonical derivative (grad, curl, div)
        """
    def Derive(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = None) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> ngsolve.fem.CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def InnerProduct(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Norm(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Operator(self, name: str) -> ProxyFunction: 
        """
        Use an additional operator of the finite element space
        """
    def Operators(self) -> list: 
        """
        returns list of available differential operators
        """
    def Other(self, bnd: object = <ngsolve.ngstd.DummyArgument>) -> ProxyFunction: 
        """
        take value from neighbour element (DG)
        """
    def Trace(self) -> ProxyFunction: 
        """
        take canonical boundary trace
        """
    @overload
    def __add__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __call__(self, mip: ngsolve.fem.BaseMappedIntegrationPoint) -> object: 
        """
        evaluate CF at a mapped integrationpoint mip. mip can be generated by calling mesh(x,y,z)
        """
    @overload
    def __call__(self, arg0: numpy.ndarray[ngsolve.fem.MeshPoint]) -> array: ...
    def __diffop__(self) -> ngsolve.fem.DifferentialOperator: ...
    @overload
    def __getitem__(self, components: tuple) -> ngsolve.fem.CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, comp: int) -> ngsolve.fem.CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __mul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    def __neg__(self) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, exponent: int) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: ngsolve.fem.CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    def __rsub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __sub__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def derivname(self) -> str:
        """
        name of the canonical derivative

        :type: str
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
    def imag(self) -> ngsolve.fem.CoefficientFunction:
        """
        imaginary part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def real(self) -> ngsolve.fem.CoefficientFunction:
        """
        real part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def trans(self) -> ngsolve.fem.CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: ngsolve.fem.CoefficientFunction
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <slot wrapper '__init__' of 'ngsolve.comp.InterpolateProxy' objects>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.InterpolateProxy' objects>, '__doc__': None, '__module__': 'ngsolve.comp'})
    pass
class L2(FESpace, NGS_Object):
    """
    An L2-conforming finite element space.

    The L2 finite element space consists of element-wise polynomials,
    which are discontinuous from element to element. It uses an
    L2-orthogonal hierarchical basis which leads to orthogonal
    mass-matrices on non-curved elements.

    Boundary values are not meaningful for an L2 function space.

    The L2 space supports element-wise variable order, which can be set
    for ELEMENT-nodes.

    Per default, all dofs are local dofs and are condensed if static
    condensation is performed. The lowest order can be kept in the
    WIRE_BASKET via the flag 'lowest_order_wb=True'.

    All dofs can be hidden. Then the basis functions don't show up in the
    global system.

     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    all_dofs_together: bool = False
      Change ordering of dofs. If this flag ist set,
      all dofs of an element are ordered successively.
      Otherwise, the lowest order dofs (the constants)
      of all elements are ordered first.
    lowest_order_wb: bool = False
      Keep lowest order dof in WIRE_BASKET
    hide_all_dofs: bool = False
      Set all used dofs to HIDDEN_DOFs
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE4128B8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.L2' objects>, '__doc__': "An L2-conforming finite element space.\n\nThe L2 finite element space consists of element-wise polynomials,\nwhich are discontinuous from element to element. It uses an\nL2-orthogonal hierarchical basis which leads to orthogonal\nmass-matrices on non-curved elements.\n\nBoundary values are not meaningful for an L2 function space.\n\nThe L2 space supports element-wise variable order, which can be set\nfor ELEMENT-nodes.\n\nPer default, all dofs are local dofs and are condensed if static\ncondensation is performed. The lowest order can be kept in the\nWIRE_BASKET via the flag 'lowest_order_wb=True'.\n\nAll dofs can be hidden. Then the basis functions don't show up in the\nglobal system.\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\nall_dofs_together: bool = False\n  Change ordering of dofs. If this flag ist set,\n  all dofs of an element are ordered successively.\n  Otherwise, the lowest order dofs (the constants)\n  of all elements are ordered first.\nlowest_order_wb: bool = False\n  Keep lowest order dof in WIRE_BASKET\nhide_all_dofs: bool = False\n  Set all used dofs to HIDDEN_DOFs\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE412918>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE412978>, '__flags_doc__': <staticmethod object at 0x0000026DEE414108>})
    pass
class LinearForm(NGS_Object):
    """
    Used to store the left hand side of a PDE. Add integrators
    (ngsolve.LFI) to it to implement your PDE.

    Parameters:

    space : ngsolve.FESpace
      The space the linearform is defined on. Can be a compound
      FESpace for a mixed formulation.

    flags : dict
      Additional options for the linearform, for example:

        print : bool
          Write additional debug information to testout file. This
          file must be set by ngsolve.SetTestoutFile. Use
          ngsolve.SetNumThreads(1) for serial output.


     Keyword arguments can be:
    print: bool
      Write additional debug information to testout file.
      This file must be set by ngsolve.SetTestoutFile. Use
      ngsolve.SetNumThreads(1) for serial output.
    printelvec: bool
      print element vectors to testout file
    """
    def Add(self, integrator: ngsolve.fem.LFI) -> LinearForm: 
        """
        Add integrator to linear form.

        Parameters:

        integrator : ngsolve.fem.LFI
          input linear form integrator
        """
    def Assemble(self) -> None: 
        """
        Assemble linear form
        """
    def __call__(self, gf: GridFunction) -> float: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    @overload
    def __iadd__(self, lfi: ngsolve.fem.LFI) -> LinearForm: ...
    @overload
    def __iadd__(self, arg0: SumOfIntegrals) -> LinearForm: ...
    def __init__(self, space: FESpace, **kwargs) -> None: ...
    def __str__(self) -> str: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> list:
        """
        list of components for linearforms on compound-space

        :type: list
        """
    @property
    def integrators(self) -> tuple:
        """
        returns tuple of integrators of the linear form

        :type: tuple
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def space(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def vec(self) -> ngsolve.la.BaseVector:
        """
        vector of the assembled linear form

        :type: ngsolve.la.BaseVector
        """
    pass
class Mesh():
    """
    NGSolve interface to the Netgen mesh. Provides access and functionality
    to use the mesh for finite element calculations.

    Parameters:

    mesh (netgen.Mesh): a mesh generated from Netgen
    """
    def BBBoundaries(self, pattern: str) -> Region: 
        """
        Return co dim 3 boundary mesh-region matching the given regex pattern
        """
    def BBoundaries(self, pattern: str) -> Region: 
        """
        Return co dim 2 boundary mesh-region matching the given regex pattern
        """
    @overload
    def Boundaries(self, pattern: str) -> Region: 
        """
        Return boundary mesh-region matching the given regex pattern

        Generate boundary mesh-region by boundary condition numbers
        """
    @overload
    def Boundaries(self, bnds: List[int]) -> Region: ...
    def Contains(self, x: float = 0.0, y: float = 0.0, z: float = 0.0) -> bool: 
        """
        Check if the point (x,y,z) is in the meshed domain (is inside a volume element)
        """
    def Curve(self, order: int) -> None: 
        """
        Curve the mesh elements for geometry approximation of given order
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> ngcomp::ElementRange: 
        """
        Return an iterator over elements on VOL/BND
        """
    def GetBBBoundaries(self) -> tuple: 
        """
        Return list of boundary conditions for co dimension 3
        """
    def GetBBoundaries(self) -> tuple: 
        """
        Return list of boundary conditions for co dimension 2
        """
    def GetBoundaries(self) -> tuple: 
        """
        Return list of boundary condition names
        """
    def GetCurveOrder(self) -> int: ...
    def GetHPElementLevel(self, ei: ElementId) -> int: 
        """
        THIS FUNCTION IS WIP!
         Return HP-refinement level of element
        """
    def GetMaterials(self) -> tuple: 
        """
        Return list of material names
        """
    def GetNE(self, arg0: VorB) -> int: 
        """
        Return number of elements of codimension VorB.
        """
    def GetPMLTrafo(self, dom: int = 1) -> pml.PML: 
        """
        Return pml transformation on domain dom
        """
    def GetPMLTrafos(self) -> list: 
        """
        Return list of pml transformations
        """
    def GetParentElement(self, ei: ElementId) -> ElementId: 
        """
        Return parent element id on refined mesh
        """
    def GetParentVertices(self, vnum: int) -> tuple: 
        """
        Return parent vertex numbers on refined mesh
        """
    def GetPeriodicNodePairs(self, arg0: ngsolve.fem.NODE_TYPE) -> list: 
        """
        returns list of periodic nodes with their identification number as [((master_nr, minion_nr),idnr),...]
        """
    def GetTrafo(self, eid: ElementId) -> ngsolve.fem.ElementTransformation: 
        """
        returns element transformation of given element id
        """
    @overload
    def MapToAllElements(self, arg0: ngsolve.fem.IntegrationRule, arg1: VorB) -> numpy.ndarray[ngsolve.fem.MeshPoint]: ...
    @overload
    def MapToAllElements(self, arg0: Dict[ngsolve.fem.ET, ngsolve.fem.IntegrationRule], arg1: VorB) -> numpy.ndarray[ngsolve.fem.MeshPoint]: ...
    @overload
    def Materials(self, pattern: str) -> Region: 
        """
        Return mesh-region matching the given regex pattern

        Generate mesh-region by domain numbers
        """
    @overload
    def Materials(self, domains: List[int]) -> Region: ...
    def Refine(self, mark_surface_elements: bool = False) -> None: 
        """
        Local mesh refinement based on marked elements, uses element-bisection algorithm
        """
    def RefineHP(self, levels: int, factor: float = 0.125) -> None: 
        """
        Geometric mesh refinement towards marked vertices and edges, uses factor for placement of new points
        """
    def SetDeformation(self, gf: ngcomp::GridFunction) -> None: 
        """
        Deform the mesh with the given GridFunction
        """
    def SetElementOrder(self, eid: ElementId, order: int) -> None: 
        """
        For backward compatibility, not recommended to use
        """
    def SetPML(self, pmltrafo: pml.PML, definedon: object) -> None: 
        """
        Set PML transformation on domain
        """
    def SetRefinementFlag(self, ei: ElementId, refine: bool) -> None: 
        """
        Set refinementflag for mesh-refinement
        """
    def SetRefinementFlags(self, refine: List[bool]) -> None: 
        """
        Set refinementflags for mesh-refinement
        """
    def UnSetPML(self, definedon: object) -> None: 
        """
        Unset PML transformation on domain
        """
    def UnsetDeformation(self) -> None: 
        """
        Unset the deformation
        """
    def __call__(self, x: numpy.ndarray[float64] = 0.0, y: numpy.ndarray[float64] = 0.0, z: numpy.ndarray[float64] = 0.0, VOL_or_BND: VorB = VorB.VOL) -> object: 
        """
        Get a MappedIntegrationPoint in the point (x,y,z) on the matching volume (VorB=VOL, default) or surface (VorB=BND) element. BBND elements aren't supported
        """
    def __eq__(self, mesh: Mesh) -> bool: ...
    @overload
    def __getitem__(self, arg0: ElementId) -> Ngs_Element: 
        """
        Return Ngs_Element from given ElementId

        Return MeshNode from given NodeId
        """
    @overload
    def __getitem__(self, arg0: NodeId) -> MeshNode: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __init__(self, ngmesh: netgen.libngpy._meshing.Mesh) -> None: 
        """
        Make an NGSolve-mesh from a Netgen-mesh

        Load a mesh from file.
        In MPI-parallel mode the mesh is distributed over the MPI-group given by the communicator (WIP!)
        """
    @overload
    def __init__(self, filename: str, comm: netgen.libngpy._meshing.MPI_Comm = <netgen.libngpy._meshing.MPI_Comm object at 0x0000026DEE18EDB0>) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def _updateBuffers(self) -> None: 
        """
        Update NGSolve mesh information, needs to be called if Netgen mesh changes
        """
    def nnodes(self, arg0: ngsolve.fem.NODE_TYPE) -> int: 
        """
        number of nodes given type
        """
    def nodes(self, node_type: ngsolve.fem.NODE_TYPE) -> MeshNodeRange: 
        """
        iterable of mesh nodes of type node_type
        """
    @property
    def comm(self) -> netgen.libngpy._meshing.MPI_Comm:
        """
        MPI-communicator the Mesh lives in

        :type: netgen.libngpy._meshing.MPI_Comm
        """
    @property
    def dim(self) -> int:
        """
        mesh dimension

        :type: int
        """
    @property
    def edges(self) -> MeshNodeRange:
        """
        iterable of mesh edges

        :type: MeshNodeRange
        """
    @property
    def faces(self) -> MeshNodeRange:
        """
        iterable of mesh faces

        :type: MeshNodeRange
        """
    @property
    def facets(self) -> MeshNodeRange:
        """
        iterable of mesh facets

        :type: MeshNodeRange
        """
    @property
    def ne(self) -> int:
        """
        number of volume elements

        :type: int
        """
    @property
    def nedge(self) -> int:
        """
        number of edges

        :type: int
        """
    @property
    def nface(self) -> int:
        """
        number of faces

        :type: int
        """
    @property
    def nfacet(self) -> int:
        """
        number of facets

        :type: int
        """
    @property
    def ngmesh(self) -> netgen.libngpy._meshing.Mesh:
        """
        the Netgen mesh

        :type: netgen.libngpy._meshing.Mesh
        """
    @property
    def nv(self) -> int:
        """
        number of vertices

        :type: int
        """
    @property
    def vertices(self) -> MeshNodeRange:
        """
        iterable of mesh vertices

        :type: MeshNodeRange
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE402BB8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.Mesh' objects>, '__doc__': '\nNGSolve interface to the Netgen mesh. Provides access and functionality\nto use the mesh for finite element calculations.\n\nParameters:\n\nmesh (netgen.Mesh): a mesh generated from Netgen\n\n\n', '__module__': 'ngsolve.comp', '__eq__': <instancemethod __eq__ at 0x0000026DEE402BE8>, 'comm': <property object at 0x0000026DEE405228>, '__getstate__': <instancemethod __getstate__ at 0x0000026DEE402C78>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE402CD8>, 'Elements': <instancemethod Elements at 0x0000026DEE402D38>, '__getitem__': <instancemethod __getitem__ at 0x0000026DEE402DC8>, 'GetNE': <instancemethod GetNE at 0x0000026DEE402DF8>, 'nv': <property object at 0x0000026DEE4054A8>, 'ne': <property object at 0x0000026DEE405598>, 'nedge': <property object at 0x0000026DEE405638>, 'nface': <property object at 0x0000026DEE4056D8>, 'nfacet': <property object at 0x0000026DEE4057C8>, 'nnodes': <instancemethod nnodes at 0x0000026DEE402F48>, 'dim': <property object at 0x0000026DEE4058B8>, 'ngmesh': <property object at 0x0000026DEE405958>, 'vertices': <property object at 0x0000026DEE405A48>, 'edges': <property object at 0x0000026DEE405B38>, 'faces': <property object at 0x0000026DEE405C28>, 'facets': <property object at 0x0000026DEE405D18>, 'nodes': <instancemethod nodes at 0x0000026DEE406108>, 'GetPeriodicNodePairs': <instancemethod GetPeriodicNodePairs at 0x0000026DEE406168>, 'GetTrafo': <instancemethod GetTrafo at 0x0000026DEE4061C8>, 'SetDeformation': <instancemethod SetDeformation at 0x0000026DEE406228>, 'UnsetDeformation': <instancemethod UnsetDeformation at 0x0000026DEE406288>, 'SetPML': <instancemethod SetPML at 0x0000026DEE4062E8>, 'UnSetPML': <instancemethod UnSetPML at 0x0000026DEE406348>, 'GetPMLTrafos': <instancemethod GetPMLTrafos at 0x0000026DEE4063A8>, 'GetPMLTrafo': <instancemethod GetPMLTrafo at 0x0000026DEE406408>, 'GetMaterials': <instancemethod GetMaterials at 0x0000026DEE406468>, 'Materials': <instancemethod Materials at 0x0000026DEE4064F8>, 'GetBoundaries': <instancemethod GetBoundaries at 0x0000026DEE406528>, 'Boundaries': <instancemethod Boundaries at 0x0000026DEE4065B8>, 'GetBBoundaries': <instancemethod GetBBoundaries at 0x0000026DEE4065E8>, 'BBoundaries': <instancemethod BBoundaries at 0x0000026DEE406648>, 'GetBBBoundaries': <instancemethod GetBBBoundaries at 0x0000026DEE4066A8>, 'BBBoundaries': <instancemethod BBBoundaries at 0x0000026DEE406708>, 'Refine': <instancemethod Refine at 0x0000026DEE406768>, 'RefineHP': <instancemethod RefineHP at 0x0000026DEE4067C8>, '_updateBuffers': <instancemethod _updateBuffers at 0x0000026DEE406828>, 'SetRefinementFlag': <instancemethod SetRefinementFlag at 0x0000026DEE406888>, 'SetRefinementFlags': <instancemethod SetRefinementFlags at 0x0000026DEE4068E8>, 'GetParentElement': <instancemethod GetParentElement at 0x0000026DEE406948>, 'GetParentVertices': <instancemethod GetParentVertices at 0x0000026DEE4069A8>, 'GetHPElementLevel': <instancemethod GetHPElementLevel at 0x0000026DEE406A08>, 'SetElementOrder': <instancemethod SetElementOrder at 0x0000026DEE406A68>, 'Curve': <instancemethod Curve at 0x0000026DEE406AC8>, 'GetCurveOrder': <instancemethod GetCurveOrder at 0x0000026DEE406B28>, 'Contains': <instancemethod Contains at 0x0000026DEE406B88>, 'MapToAllElements': <instancemethod MapToAllElements at 0x0000026DEE406C18>, '__call__': <instancemethod __call__ at 0x0000026DEE406C48>})
    pass
class NodeId():
    """
    an node identifier containing node type and node nr
    """
    def __eq__(self, arg0: NodeId) -> bool: ...
    def __hash__(self) -> int: ...
    def __init__(self, type: ngsolve.fem.NODE_TYPE, nr: int) -> None: ...
    def __ne__(self, arg0: NodeId) -> bool: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
    @property
    def nr(self) -> int:
        """
        the node number

        :type: int
        """
    @property
    def type(self) -> ngsolve.fem.NODE_TYPE:
        """
        the node type

        :type: ngsolve.fem.NODE_TYPE
        """
    pass
class MeshNodeRange():
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: ...
    pass
class Preconditioner(ngsolve.la.BaseMatrix, NGS_Object):
    """
     Keyword arguments can be:
    inverse: 
      Inverse type used in Preconditioner.
    test: bool = False
      Computes condition number for preconditioner, if testout file
      is set, prints eigenvalues to file.
    """
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
    def Mult(self, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
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
    def Test(self) -> None: ...
    def Update(self) -> None: 
        """
        Update preconditioner
        """
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __iadd__(self, mat: ngsolve.la.BaseMatrix) -> None: ...
    def __init__(self, bf: BilinearForm, type: str, **kwargs) -> None: ...
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
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
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
    def mat(self) -> ngsolve.la.BaseMatrix:
        """
        matrix of the preconditioner

        :type: ngsolve.la.BaseMatrix
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
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
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE421D38>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.Preconditioner' objects>, '__doc__': '\n Keyword arguments can be:\ninverse: \n  Inverse type used in Preconditioner.\ntest: bool = False\n  Computes condition number for preconditioner, if testout file\n  is set, prints eigenvalues to file.\n', '__module__': 'ngsolve.comp', '__flags_doc__': <staticmethod object at 0x0000026DEE41CC88>, 'Test': <instancemethod Test at 0x0000026DEE421DC8>, 'Update': <instancemethod Update at 0x0000026DEE421E28>, 'mat': <property object at 0x0000026DEE422F48>})
    pass
class BilinearForm(NGS_Object):
    """
    Used to store the left hand side of a PDE. integrators (ngsolve.BFI)
    to it to implement your PDE. If the left hand side is linear
    you can use BilinearForm.Assemble to assemble it after adding
    your integrators. For nonlinear usage use BilinearForm.Apply or
    BilinearForm.AssembleLinearization instead of Bilinearform.Assemble.

    Parameters:

    space : ngsolve.FESpace
      The finite element space the bilinearform is defined on. This
      can be a compound FESpace for a mixed formulation.


     Keyword arguments can be:
    condense: bool = False
      (formerly known as 'eliminate_internal')
      Set up BilinearForm for static condensation of internal
      bubbles. Static condensation has to be done by user,
      this enables only the use of the members harmonic_extension,
      harmonic_extension_trans and inner_solve. Have a look at the
      documentation for further information.
    eliminate_internal: bool = False
      deprecated for static condensation, replaced by 'condense'

    eliminate_hidden: bool = False
      Set up BilinearForm for static condensation of hidden
      dofs. May be overruled by eliminate_internal.
    print: bool = False
      Write additional information to testout file. 
      This file must be set by ngsolve.SetTestoutFile. Use 
      ngsolve.SetNumThreads(1) for serial output
    printelmat: bool = False
      Write element matrices to testout file
    symmetric: bool = False
      If set true, only half the matrix is stored
    nonassemble: bool = False
      BilinearForm will not allocate memory for assembling.
      optimization feature for (nonlinear) problems where the
      form is only applied but never assembled.
    project: bool = False
      When calling bf.Assemble, all saved coarse matrices from
      mesh refinements are updated as well using a Galerkin projection
      of the matrix on the finest grid. This is needed to use the multigrid
      preconditioner with a changing bilinearform.
    nonsym_storage: bool = False
      The full matrix is stored, even if the symmetric flag is set.
    diagonal: bool = False
      Stores only the diagonal of the matrix.
    geom_free: bool = False
      when element matrices are independent of geometry, we store them 
      only for the referecne elements
    check_unused: bool = True
      If set prints warnings if not UNUSED_DOFS are not used.
    """
    def Add(self, integrator: ngsolve.fem.BFI) -> BilinearForm: 
        """
                 Add integrator to bilinear form.

        Parameters:

        integrator : ngsolve.fem.BFI
          input bilinear form integrator
        """
    def Apply(self, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: 
        """
        Applies a (non-)linear variational formulation to x and stores the result in y.

        Parameters:

        x : ngsolve.BaseVector
          input vector

        y : ngsolve.BaseVector
          output vector
        """
    def Assemble(self, reallocate: bool = False) -> None: 
        """
        Assemble the bilinear form.

        Parameters:

        reallocate : bool
          input reallocate
        """
    def AssembleLinearization(self, ulin: ngsolve.la.BaseVector, reallocate: bool = False) -> None: 
        """
        Computes linearization of the bilinear form at given vecor.

        Parameters:

        ulin : ngsolve.la.BaseVector
          input vector
        """
    def ComputeInternal(self, u: ngsolve.la.BaseVector, f: ngsolve.la.BaseVector) -> None: 
        """
        Parameters:

        u : ngsolve.la.BaseVector
          input vector

        f : ngsolve.la.BaseVector
          input right hand side
        """
    def DeleteMatrix(self) -> None: ...
    def DeleteSpecialElements(self) -> None: ...
    def Energy(self, x: ngsolve.la.BaseVector) -> float: 
        """
        Computes the energy of EnergyIntegrators like SymbolicEnergy for given input vector.

        Parameters:

        x : ngsolve.la.BaseVector
          input vector
        """
    def Flux(self, gf: GridFunction) -> ngsolve.fem.CoefficientFunction: 
        """
        Parameters:

        gf : ngsolve.comp.GridFunction
          input GridFunction
        """
    def SetPreconditioner(self, arg0: ngcomp::Preconditioner) -> None: ...
    def UnsetPreconditioner(self, arg0: ngcomp::Preconditioner) -> None: ...
    def __call__(self, gfu: GridFunction, gfv: GridFunction) -> float: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    @overload
    def __iadd__(self, arg0: SumOfIntegrals) -> BilinearForm: ...
    @overload
    def __iadd__(self, other: ngsolve.fem.BFI) -> BilinearForm: ...
    @overload
    def __iadd__(self, arg0: Variation) -> BilinearForm: ...
    @overload
    def __init__(self, trialspace: FESpace, testspace: FESpace, **kwargs) -> None: ...
    @overload
    def __init__(self, space: FESpace, **kwargs) -> None: ...
    def __str__(self) -> str: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> list:
        """
        list of components for bilinearforms on compound-space

        :type: list
        """
    @property
    def condense(self) -> bool:
        """
        use static condensation ?

        :type: bool
        """
    @property
    def harmonic_extension(self) -> ngsolve.la.BaseMatrix:
        """
        harmonic_extension used for static condensaition

        :type: ngsolve.la.BaseMatrix
        """
    @property
    def harmonic_extension_trans(self) -> ngsolve.la.BaseMatrix:
        """
        harmonic_extension_trans used for static condensation

        :type: ngsolve.la.BaseMatrix
        """
    @property
    def inner_matrix(self) -> ngsolve.la.BaseMatrix:
        """
        inner_matrix of the bilinear form

        :type: ngsolve.la.BaseMatrix
        """
    @property
    def inner_solve(self) -> ngsolve.la.BaseMatrix:
        """
        inner_solve used for static condensation

        :type: ngsolve.la.BaseMatrix
        """
    @property
    def integrators(self) -> tuple:
        """
        integrators of the bilinear form

        :type: tuple
        """
    @property
    def loform(self) -> BilinearForm:
        """
        :type: BilinearForm
        """
    @property
    def mat(self) -> ngsolve.la.BaseMatrix:
        """
        matrix of the assembled bilinear form

        :type: ngsolve.la.BaseMatrix
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def space(self) -> FESpace:
        """
        fespace on which the bilinear form is defined on

        :type: FESpace
        """
    pass
class FESpaceElement(Ngs_Element):
    def GetFE(self) -> ngsolve.fem.FiniteElement: 
        """
        the finite element containing shape functions
        """
    def GetTrafo(self) -> ngsolve.fem.ElementTransformation: 
        """
        the transformation from reference element to physical element
        """
    def VB(self) -> VorB: 
        """
        VorB of element
        """
    @property
    def dofs(self) -> list:
        """
        degrees of freedom of element

        :type: list
        """
    @property
    def edges(self) -> tuple:
        """
        tuple of global edge numbers

        :type: tuple
        """
    @property
    def elementnode(self) -> NodeId:
        """
        inner node, i.e. cell, face or edge node for 3D/2D/1D

        :type: NodeId
        """
    @property
    def faces(self) -> tuple:
        """
        tuple of global face numbers

        :type: tuple
        """
    @property
    def facets(self) -> tuple:
        """
        tuple of global face, edge or vertex numbers

        :type: tuple
        """
    @property
    def index(self) -> int:
        """
        material or boundary condition index

        :type: int
        """
    @property
    def mat(self) -> str:
        """
        material or boundary condition label

        :type: str
        """
    @property
    def nr(self) -> int:
        """
        the element number

        :type: int
        """
    @property
    def type(self) -> ngsolve.fem.ET:
        """
        geometric shape of element

        :type: ngsolve.fem.ET
        """
    @property
    def valid(self) -> bool:
        """
        is element valid

        :type: bool
        """
    @property
    def vertices(self) -> tuple:
        """
        tuple of global vertex numbers

        :type: tuple
        """
    pass
class NodalFESpace(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416A68>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.NodalFESpace' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE416AC8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE416B28>, '__flags_doc__': <staticmethod object at 0x0000026DEE414C88>})
    pass
class MeshNode(NodeId):
    """
    a node within a mesh
    """
    def __eq__(self, arg0: NodeId) -> bool: ...
    def __hash__(self) -> int: ...
    def __ne__(self, arg0: NodeId) -> bool: ...
    def __repr__(self) -> str: ...
    def __str__(self) -> str: ...
    @property
    def edges(self) -> tuple:
        """
        tuple of global edge numbers

        :type: tuple
        """
    @property
    def elements(self) -> tuple:
        """
        tuple of global element-ids

        :type: tuple
        """
    @property
    def faces(self) -> tuple:
        """
        tuple of global face numbers

        :type: tuple
        """
    @property
    def nr(self) -> int:
        """
        the node number

        :type: int
        """
    @property
    def point(self) -> tuple:
        """
        vertex coordinates

        :type: tuple
        """
    @property
    def type(self) -> ngsolve.fem.NODE_TYPE:
        """
        the node type

        :type: ngsolve.fem.NODE_TYPE
        """
    @property
    def vertices(self) -> tuple:
        """
        tuple of global vertex numbers

        :type: tuple
        """
    pass
class NodeRange():
    def __iter__(self) -> iterator: ...
    def __len__(self) -> int: ...
    pass
class NormalFacetFESpace(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    highest_order_dc: bool = False
      Splits highest order facet functions into two which are associated with
      the corresponding neighbors and are local dofs on the corresponding element
     (used to realize projected jumps)
    hide_highest_order_dc: bool = False
      if highest_order_dc is used this flag marks the corresponding local dofs
      as hidden dofs (reduces number of non-zero entries in a matrix). These dofs
      can also be compressed.
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE4160D8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.NormalFacetFESpace' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\nhighest_order_dc: bool = False\n  Splits highest order facet functions into two which are associated with\n  the corresponding neighbors and are local dofs on the corresponding element\n (used to realize projected jumps)\nhide_highest_order_dc: bool = False\n  if highest_order_dc is used this flag marks the corresponding local dofs\n  as hidden dofs (reduces number of non-zero entries in a matrix). These dofs\n  can also be compressed.\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE416138>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE416198>, '__flags_doc__': <staticmethod object at 0x0000026DEE414648>})
    pass
class NumProc(NGS_Object):
    def Do(self) -> None: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    pass
class NumberSpace(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412768>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.NumberSpace' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE4127C8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE412828>, '__flags_doc__': <staticmethod object at 0x0000026DEE40EF08>})
    pass
class ORDER_POLICY():
    """
    Enumeration of all supported order policies

    Members:

      CONSTANT

      NODETYPE

      VARIABLE

      OLDSTYLE
    """
    def __init__(self, arg0: int) -> None: ...
    def __int__(self) -> int: ...
    @property
    def name(self) -> str:
        """
        (self: handle) -> str

        :type: str
        """
    CONSTANT: ngsolve.comp.ORDER_POLICY # value = ORDER_POLICY.CONSTANT
    NODETYPE: ngsolve.comp.ORDER_POLICY # value = ORDER_POLICY.NODETYPE
    OLDSTYLE: ngsolve.comp.ORDER_POLICY # value = ORDER_POLICY.OLDSTYLE
    VARIABLE: ngsolve.comp.ORDER_POLICY # value = ORDER_POLICY.VARIABLE
    __entries: dict # value = {'CONSTANT': (ORDER_POLICY.CONSTANT, None), 'NODETYPE': (ORDER_POLICY.NODETYPE, None), 'VARIABLE': (ORDER_POLICY.VARIABLE, None), 'OLDSTYLE': (ORDER_POLICY.OLDSTYLE, None)}
    __members__: dict # value = {'CONSTANT': ORDER_POLICY.CONSTANT, 'NODETYPE': ORDER_POLICY.NODETYPE, 'VARIABLE': ORDER_POLICY.VARIABLE, 'OLDSTYLE': ORDER_POLICY.OLDSTYLE}
    pass
class PDE():
    @overload
    def Add(self, np: NumProc) -> None: ...
    @overload
    def Add(self, cf: ngsolve.fem.CoefficientFunction, name: str) -> None: ...
    @overload
    def Add(self, space: FESpace) -> None: ...
    @overload
    def Add(self, lf: LinearForm) -> None: ...
    @overload
    def Add(self, list: list) -> None: ...
    @overload
    def Add(self, name: str, value: float) -> None: ...
    @overload
    def Add(self, gf: GridFunction) -> None: ...
    @overload
    def Add(self, bf: BilinearForm) -> None: ...
    @overload
    def Add(self, pre: Preconditioner) -> None: ...
    @overload
    def Add(self, mesh: Mesh) -> None: ...
    def AddPDE_File(self, filename: str) -> None: 
        """
        Adds definitions of other PDE file into existing one
        """
    def LoadSolution(self, filename: str, ascii: bool = False) -> None: ...
    def Mesh(self, meshnr: int = 0) -> Mesh: ...
    def SetCurveIntegrator(self, filename: str, lfi: ngsolve.fem.LFI) -> None: ...
    def Solve(self) -> None: ...
    @overload
    def __init__(self, filename: str) -> None: ...
    @overload
    def __init__(self) -> None: ...
    def __str__(self) -> str: ...
    @property
    def bilinearforms(self) -> object:
        """
        :type: object
        """
    @property
    def coefficients(self) -> object:
        """
        :type: object
        """
    @property
    def constants(self) -> object:
        """
        :type: object
        """
    @property
    def gridfunctions(self) -> object:
        """
        :type: object
        """
    @property
    def linearforms(self) -> object:
        """
        :type: object
        """
    @property
    def numprocs(self) -> object:
        """
        :type: object
        """
    @property
    def preconditioners(self) -> object:
        """
        :type: object
        """
    @property
    def spaces(self) -> object:
        """
        :type: object
        """
    @property
    def variables(self) -> object:
        """
        :type: object
        """
    pass
class Periodic(FESpace, NGS_Object):
    """
    Periodic or quasi-periodic Finite Element Spaces.
    The periodic fespace is a wrapper around a standard fespace with an 
    additional dof mapping for the periodic degrees of freedom. All dofs 
    on minion boundaries are mapped to their master dofs. Because of this, 
    the mesh needs to be periodic. Low order fespaces are currently not
    supported, so methods using them will not work.

    Parameters:

    fespace : ngsolve.comp.FESpace
        finite element space

    phase : list of Complex = None
        phase shift for quasi-periodic finite element space. The basis
        functions on the minion boundary are multiplied by the factor
        given in this list. If None (default) is given, a periodic
        fespace is created. The order of the list must match the order
        of the definition of the periodic boundaries in the mesh.

    used_idnrs : list of int = None
        identification numbers to be made periodic if you don't want to
        use all periodic identifications defined in the mesh, if None
        (default) all available periodic identifications are used.


     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, fespace: FESpace, phase: Optional[list] = None, use_idnrs: object = []) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416D08>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.Periodic' objects>, '__doc__': "Periodic or quasi-periodic Finite Element Spaces.\nThe periodic fespace is a wrapper around a standard fespace with an \nadditional dof mapping for the periodic degrees of freedom. All dofs \non minion boundaries are mapped to their master dofs. Because of this, \nthe mesh needs to be periodic. Low order fespaces are currently not\nsupported, so methods using them will not work.\n\nParameters:\n\nfespace : ngsolve.comp.FESpace\n    finite element space\n\nphase : list of Complex = None\n    phase shift for quasi-periodic finite element space. The basis\n    functions on the minion boundary are multiplied by the factor\n    given in this list. If None (default) is given, a periodic\n    fespace is created. The order of the list must match the order\n    of the definition of the periodic boundaries in the mesh.\n\nused_idnrs : list of int = None\n    identification numbers to be made periodic if you don't want to\n    use all periodic identifications defined in the mesh, if None\n    (default) all available periodic identifications are used.\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE416D68>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE416DC8>})
    pass
class MultiGridPreconditioner(Preconditioner, ngsolve.la.BaseMatrix, NGS_Object):
    """
     Keyword arguments can be:
    inverse: 
      Inverse type used in Preconditioner.
    test: bool = False
      Computes condition number for preconditioner, if testout file
      is set, prints eigenvalues to file.
    updateall: bool = False
      Update all smoothing levels when calling Update
    smoother: string = 'point'
      Smoother between multigrid levels, available options are:
        'point': Gauss-Seidel-Smoother
        'line':  Anisotropic smoother
        'block': Block smoother
    coarsetype: string = direct
      How to solve coarse problem.
    coarsesmoothingsteps: int = 1
      If coarsetype is smoothing, then how many smoothingsteps will be done.
    updatealways: bool = False
    """
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
    def Mult(self, x: ngsolve.la.BaseVector, y: ngsolve.la.BaseVector) -> None: ...
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
    def Test(self) -> None: ...
    def Update(self) -> None: 
        """
        Update preconditioner
        """
    def __add__(self, mat: BaseMatrix) -> BaseMatrix: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __iadd__(self, mat: ngsolve.la.BaseMatrix) -> None: ...
    def __init__(self, bf: BilinearForm, **kwargs) -> None: ...
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
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
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
    def mat(self) -> ngsolve.la.BaseMatrix:
        """
        matrix of the preconditioner

        :type: ngsolve.la.BaseMatrix
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
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
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE421EB8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.MultiGridPreconditioner' objects>, '__doc__': "\n Keyword arguments can be:\ninverse: \n  Inverse type used in Preconditioner.\ntest: bool = False\n  Computes condition number for preconditioner, if testout file\n  is set, prints eigenvalues to file.\nupdateall: bool = False\n  Update all smoothing levels when calling Update\nsmoother: string = 'point'\n  Smoother between multigrid levels, available options are:\n    'point': Gauss-Seidel-Smoother\n    'line':  Anisotropic smoother\n    'block': Block smoother\ncoarsetype: string = direct\n  How to solve coarse problem.\ncoarsesmoothingsteps: int = 1\n  If coarsetype is smoothing, then how many smoothingsteps will be done.\nupdatealways: bool = False\n\n", '__module__': 'ngsolve.comp', '__flags_doc__': <staticmethod object at 0x0000026DEE41CA48>})
    pass
class Prolongation():
    def Prolongate(self, finelevel: int, vec: ngsolve.la.BaseVector) -> None: ...
    def Restrict(self, finelevel: int, vec: ngsolve.la.BaseVector) -> None: ...
    pass
class DualProxyFunction(ProxyFunction, ngsolve.fem.CoefficientFunction):
    def Compile(self, realcompile: bool = False, maxderiv: int = 2, wait: bool = False) -> ngsolve.fem.CoefficientFunction: 
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
    def Deriv(self) -> ProxyFunction: 
        """
        take canonical derivative (grad, curl, div)
        """
    def Derive(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, variable: ngsolve.fem.CoefficientFunction, direction: ngsolve.fem.CoefficientFunction = None) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute directional derivative with respect to variable
        """
    def DiffShape(self, direction: ngsolve.fem.CoefficientFunction = 1.0) -> ngsolve.fem.CoefficientFunction: 
        """
        Compute shape derivative in direction
        """
    def Eig(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns eigenvectors and eigenvalues of matrix-valued CF
        """
    def Freeze(self) -> ngsolve.fem.CoefficientFunction: 
        """
        don't differentiate this expression
        """
    def InnerProduct(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: 
        """
         
        Returns InnerProduct with another CoefficientFunction.

        Parameters:

        cf : ngsolve.CoefficientFunction
          input CoefficientFunction

         
        """
    def Norm(self) -> ngsolve.fem.CoefficientFunction: 
        """
        Returns Norm of the CF
        """
    def Operator(self, name: str) -> ProxyFunction: 
        """
        Use an additional operator of the finite element space
        """
    def Operators(self) -> list: 
        """
        returns list of available differential operators
        """
    def Other(self, bnd: object = <ngsolve.ngstd.DummyArgument>) -> ProxyFunction: 
        """
        take value from neighbour element (DG)
        """
    def Trace(self) -> ProxyFunction: 
        """
        take canonical boundary trace
        """
    @overload
    def __add__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __add__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    def __call__(self, arg0: ngsolve.fem.CoefficientFunction) -> ngfem::SumOfIntegrals: ...
    def __diffop__(self) -> ngsolve.fem.DifferentialOperator: ...
    @overload
    def __getitem__(self, components: tuple) -> ngsolve.fem.CoefficientFunction: 
        """
        returns component comp of vectorial CF
        """
    @overload
    def __getitem__(self, comp: int) -> ngsolve.fem.CoefficientFunction: ...
    def __getstate__(self) -> tuple: ...
    @overload
    def __mul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __mul__(self, arg0: ngfem::DifferentialSymbol) -> ngfem::SumOfIntegrals: ...
    def __neg__(self) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: float) -> object: ...
    @overload
    def __pow__(self, exponent: int) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __pow__(self, arg0: ngsolve.fem.CoefficientFunction) -> object: ...
    @overload
    def __radd__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __radd__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rmul__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    def __rsub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __rtruediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    def __str__(self) -> str: ...
    @overload
    def __sub__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __sub__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: float) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, value: complex) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __truediv__(self, cf: ngsolve.fem.CoefficientFunction) -> ngsolve.fem.CoefficientFunction: ...
    @property
    def data(self) -> dict:
        """
        :type: dict
        """
    @property
    def derivname(self) -> str:
        """
        name of the canonical derivative

        :type: str
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
    def imag(self) -> ngsolve.fem.CoefficientFunction:
        """
        imaginary part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def is_complex(self) -> bool:
        """
        is CoefficientFunction complex-valued ?

        :type: bool
        """
    @property
    def real(self) -> ngsolve.fem.CoefficientFunction:
        """
        real part of CF

        :type: ngsolve.fem.CoefficientFunction
        """
    @property
    def trans(self) -> ngsolve.fem.CoefficientFunction:
        """
        transpose of matrix-valued CF

        :type: ngsolve.fem.CoefficientFunction
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <slot wrapper '__init__' of 'ngsolve.comp.DualProxyFunction' objects>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.DualProxyFunction' objects>, '__doc__': None, '__module__': 'ngsolve.comp', '__call__': <instancemethod __call__ at 0x0000026DEE40AB58>})
    pass
class PyNumProc(NumProc, NGS_Object):
    def Do(self, lh: ngsolve.ngstd.LocalHeap) -> None: ...
    def __init__(self, pde: ngcomp::PDE, flags: pyngcore.Flags) -> None: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def pde(self) -> ngcomp::PDE:
        """
        PDE of the NumProc

        :type: ngcomp::PDE
        """
    pass
class Region():
    """
    a subset of volume or boundary elements
    """
    def Mask(self) -> pyngcore.BitArray: 
        """
        BitArray mask of the region
        """
    def VB(self) -> VorB: 
        """
        VorB of the region
        """
    @overload
    def __add__(self, arg0: Region) -> Region: ...
    @overload
    def __add__(self, arg0: str) -> Region: ...
    @overload
    def __init__(self, mesh: ngcomp::MeshAccess, vb: VorB, name: str) -> None: ...
    @overload
    def __init__(self, mesh: ngcomp::MeshAccess, vb: VorB, mask: pyngcore.BitArray) -> None: ...
    def __invert__(self) -> Region: ...
    @overload
    def __sub__(self, arg0: Region) -> Region: ...
    @overload
    def __sub__(self, arg0: str) -> Region: ...
    pass
class Reorder(FESpace, NGS_Object):
    """
    Reordered Finite Element Spaces.
    ...

     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, fespace: FESpace) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416EB8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.Reorder' objects>, '__doc__': "Reordered Finite Element Spaces.\n...\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp'})
    pass
class SumOfIntegrals():
    def Compile(self, realcompile: bool = False, wait: bool = False) -> SumOfIntegrals: ...
    def Derive(self, arg0: ngsolve.fem.CoefficientFunction, arg1: ngsolve.fem.CoefficientFunction) -> SumOfIntegrals: 
        """
        depricated: use 'Diff' instead
        """
    def Diff(self, arg0: ngsolve.fem.CoefficientFunction, arg1: ngsolve.fem.CoefficientFunction) -> SumOfIntegrals: ...
    def DiffShape(self, arg0: ngsolve.fem.CoefficientFunction) -> SumOfIntegrals: ...
    def __add__(self, arg0: SumOfIntegrals) -> SumOfIntegrals: ...
    def __getitem__(self, arg0: int) -> Integral: ...
    def __len__(self) -> int: ...
    def __rmul__(self, arg0: float) -> SumOfIntegrals: ...
    def __str__(self) -> str: ...
    def __sub__(self, arg0: SumOfIntegrals) -> SumOfIntegrals: ...
    pass
class SurfaceL2(FESpace, NGS_Object):
    """
    An L2-conforming finite element space.

    The L2 finite element space on surfaces consists of element-wise polynomials,
    which are discontinuous from element to element. It uses an
    L2-orthogonal hierarchical basis which leads to orthogonal
    mass-matrices on non-curved elements.

    The L2 space supports element-wise variable order, which can be set
    for ELEMENT-nodes.

    Per default, all dofs are local dofs and are condensed if static
    condensation is performed. The lowest order can be kept in the
    WIRE_BASKET via the flag 'lowest_order_wb=True'.


     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    lowest_order_wb: bool = False
      Keep lowest order dof in WIRE_BASKET and make other dofs LOCAL
    discontinuous: bool = False
      Make all dofs LOCAL
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412618>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.SurfaceL2' objects>, '__doc__': "An L2-conforming finite element space.\n\nThe L2 finite element space on surfaces consists of element-wise polynomials,\nwhich are discontinuous from element to element. It uses an\nL2-orthogonal hierarchical basis which leads to orthogonal\nmass-matrices on non-curved elements.\n\nThe L2 space supports element-wise variable order, which can be set\nfor ELEMENT-nodes.\n\nPer default, all dofs are local dofs and are condensed if static\ncondensation is performed. The lowest order can be kept in the\nWIRE_BASKET via the flag 'lowest_order_wb=True'.\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\nlowest_order_wb: bool = False\n  Keep lowest order dof in WIRE_BASKET and make other dofs LOCAL\ndiscontinuous: bool = False\n  Make all dofs LOCAL\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE412678>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE4126D8>, '__flags_doc__': <staticmethod object at 0x0000026DEE40EB48>})
    pass
class SymbolTable_D():
    def GetName(self, arg0: int) -> str: ...
    def __contains__(self, arg0: str) -> bool: ...
    @overload
    def __getitem__(self, pos: int) -> float: ...
    @overload
    def __getitem__(self, name: str) -> float: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    pass
class SymbolTable_sp_D():
    def GetName(self, pos: int) -> str: ...
    def __contains__(self, arg0: str) -> bool: ...
    @overload
    def __getitem__(self, name: str) -> float: ...
    @overload
    def __getitem__(self, pos: int) -> float: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    pass
class SymbolTable_sp_class ngcomp::BilinearForm():
    def GetName(self, arg0: int) -> str: ...
    def __contains__(self, arg0: str) -> bool: ...
    @overload
    def __getitem__(self, name: str) -> BilinearForm: ...
    @overload
    def __getitem__(self, pos: int) -> BilinearForm: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    pass
class SymbolTable_sp_class ngcomp::FESpace():
    def GetName(self, arg0: int) -> str: ...
    def __contains__(self, arg0: str) -> bool: ...
    @overload
    def __getitem__(self, pos: int) -> FESpace: ...
    @overload
    def __getitem__(self, name: str) -> FESpace: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    pass
class SymbolTable_sp_class ngcomp::GridFunction():
    def GetName(self, arg0: int) -> str: ...
    def __contains__(self, arg0: str) -> bool: ...
    @overload
    def __getitem__(self, pos: int) -> GridFunction: ...
    @overload
    def __getitem__(self, name: str) -> GridFunction: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    pass
class SymbolTable_sp_class ngcomp::LinearForm():
    def GetName(self, arg0: int) -> str: ...
    def __contains__(self, arg0: str) -> bool: ...
    @overload
    def __getitem__(self, pos: int) -> LinearForm: ...
    @overload
    def __getitem__(self, name: str) -> LinearForm: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    pass
class SymbolTable_sp_class ngcomp::NumProc():
    def GetName(self, arg0: int) -> str: ...
    def __contains__(self, arg0: str) -> bool: ...
    @overload
    def __getitem__(self, pos: int) -> NumProc: ...
    @overload
    def __getitem__(self, name: str) -> NumProc: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    pass
class SymbolTable_sp_class ngcomp::Preconditioner():
    def GetName(self, arg0: int) -> str: ...
    def __contains__(self, arg0: str) -> bool: ...
    @overload
    def __getitem__(self, name: str) -> Preconditioner: ...
    @overload
    def __getitem__(self, pos: int) -> Preconditioner: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    pass
class SymbolTable_sp_class ngfem::CoefficientFunction():
    def GetName(self, arg0: int) -> str: ...
    def __contains__(self, arg0: str) -> bool: ...
    @overload
    def __getitem__(self, name: str) -> ngsolve.fem.CoefficientFunction: ...
    @overload
    def __getitem__(self, pos: int) -> ngsolve.fem.CoefficientFunction: ...
    def __len__(self) -> int: ...
    def __str__(self) -> str: ...
    pass
class TangentialFacetFESpace(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    highest_order_dc: bool = False
      Splits highest order facet functions into two which are associated with
      the corresponding neighbors and are local dofs on the corresponding element
     (used to realize projected jumps)
    hide_highest_order_dc: bool = False
      if highest_order_dc is used this flag marks the corresponding local dofs
      as hidden dofs (reduces number of non-zero entries in a matrix). These dofs
      can also be compressed.
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412F48>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.TangentialFacetFESpace' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\nhighest_order_dc: bool = False\n  Splits highest order facet functions into two which are associated with\n  the corresponding neighbors and are local dofs on the corresponding element\n (used to realize projected jumps)\nhide_highest_order_dc: bool = False\n  if highest_order_dc is used this flag marks the corresponding local dofs\n  as hidden dofs (reduces number of non-zero entries in a matrix). These dofs\n  can also be compressed.\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE412FA8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE416048>, '__flags_doc__': <staticmethod object at 0x0000026DEE414588>})
    pass
class VTKOutput():
    @overload
    def Do(self, vb: VorB = VorB.VOL) -> None: ...
    @overload
    def Do(self, vb: VorB = VorB.VOL, drawelems: pyngcore.BitArray) -> None: ...
    def __init__(self, ma: Mesh, coefs: list = [], names: list = [], filename: str = 'vtkout', subdivision: int = 0, only_element: int = -1) -> None: ...
    pass
class Variation():
    def Compile(self, realcompile: bool = False, wait: bool = False) -> Variation: ...
    def __init__(self, arg0: SumOfIntegrals) -> None: ...
    pass
class VectorFacetFESpace(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE4167C8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.VectorFacetFESpace' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE416828>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE416888>, '__flags_doc__': <staticmethod object at 0x0000026DEE414AC8>})
    pass
class VectorFacetSurface(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416918>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.VectorFacetSurface' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE416978>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE4169D8>, '__flags_doc__': <staticmethod object at 0x0000026DEE414B88>})
    pass
class VectorH1(CompoundFESpace, FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    dirichletx: regexpr
      Regular expression string defining the dirichlet boundary
      on the first component of VectorH1.
      More than one boundary can be combined by the | operator,
      i.e.: dirichletx = 'top|right'
    dirichlety: regexpr
      Dirichlet boundary for the second component
    dirichletz: regexpr
      Dirichlet boundary for the third component
    dirichletx_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D, on the first component.
      More than one bboundary can be combined by the | operator,
      i.e.: dirichletx_bbnd = 'top|right'
    dirichlety_bbnd: regexpr
      Dirichlet bboundary for the second component
    dirichletz_bbnd: regexpr
      Dirichlet bboundary for the third component
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE412378>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.VectorH1' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\ndirichletx: regexpr\n  Regular expression string defining the dirichlet boundary\n  on the first component of VectorH1.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichletx = 'top|right'\ndirichlety: regexpr\n  Dirichlet boundary for the second component\ndirichletz: regexpr\n  Dirichlet boundary for the third component\ndirichletx_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D, on the first component.\n  More than one bboundary can be combined by the | operator,\n  i.e.: dirichletx_bbnd = 'top|right'\ndirichlety_bbnd: regexpr\n  Dirichlet bboundary for the second component\ndirichletz_bbnd: regexpr\n  Dirichlet bboundary for the third component\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE4123D8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE412438>, '__flags_doc__': <staticmethod object at 0x0000026DEE40E648>})
    pass
class VectorL2(CompoundFESpace, FESpace, NGS_Object):
    """
    A vector-valued L2-conforming finite element space.

    The Vector-L2 finite element space is a product-space of L2 spaces,
    where the number of components coincides with the mesh dimension.

    It is implemented by means of a CompoundFESpace, as one could do it at the
    user-level. Additionally, some operators are added for convenience and performance:
    One can evaluate the vector-valued function, and one can take the gradient.

     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    piola: bool = False
      Use Piola transform to map to physical element
      allows to use the div-differential operator.
    covariant: bool = False
      Use the covariant transform to map to physical element
      allows to use the curl-differential operator.
    all_dofs_together: bool = True
      dofs within one scalar component are together.
    hide_all_dofs: bool = False
      all dofs are condensed without a global dofnr
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE4124C8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.VectorL2' objects>, '__doc__': "A vector-valued L2-conforming finite element space.\n\nThe Vector-L2 finite element space is a product-space of L2 spaces,\nwhere the number of components coincides with the mesh dimension.\n\nIt is implemented by means of a CompoundFESpace, as one could do it at the\nuser-level. Additionally, some operators are added for convenience and performance:\nOne can evaluate the vector-valued function, and one can take the gradient.\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\npiola: bool = False\n  Use Piola transform to map to physical element\n  allows to use the div-differential operator.\ncovariant: bool = False\n  Use the covariant transform to map to physical element\n  allows to use the curl-differential operator.\nall_dofs_together: bool = True\n  dofs within one scalar component are together.\nhide_all_dofs: bool = False\n  all dofs are condensed without a global dofnr\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE412528>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE412588>, '__flags_doc__': <staticmethod object at 0x0000026DEE40EC88>})
    pass
class VectorNodalFESpace(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416BB8>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.VectorNodalFESpace' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE416C18>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE416C78>, '__flags_doc__': <staticmethod object at 0x0000026DEE414D48>})
    pass
class VectorSurfaceL2(FESpace, NGS_Object):
    """
     Keyword arguments can be:
    order: int = 1
      order of finite element space
    complex: bool = False
      Set if FESpace should be complex
    dirichlet: regexpr
      Regular expression string defining the dirichlet boundary.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet = 'top|right'
    dirichlet_bbnd: regexpr
      Regular expression string defining the dirichlet bboundary,
      i.e. points in 2D and edges in 3D.
      More than one boundary can be combined by the | operator,
      i.e.: dirichlet_bbnd = 'top|right'
    definedon: Region or regexpr
      FESpace is only defined on specific Region, created with mesh.Materials('regexpr')
      or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be
      mesh.Materials('regexpr').
    dim: int = 1
      Create multi dimensional FESpace (i.e. [H1]^3)
    dgjumps: bool = False
      Enable discontinuous space for DG methods, this flag is needed for DG methods,
      since the dofs have a different coupling then and this changes the sparsity
      pattern of matrices.
    low_order_space: bool = True
      Generate a lowest order space together with the high-order space,
      needed for some preconditioners.
    order_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE
      CONSTANT .. use the same fixed order for all elements,
      NODAL ..... use the same order for nodes of same shape,
      VARIBLE ... use an individual order for each edge, face and cell,
      OLDSTYLE .. as it used to be for the last decade
    """
    def ApplyM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
        Apply mass-matrix. Available only for L2-like spaces
        """
    def ConvertL2Operator(self, l2space: FESpace) -> ngsolve.la.BaseMatrix: ...
    def CouplingType(self, dofnr: int) -> COUPLING_TYPE: 
        """
                 Get coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number
        """
    def Elements(self, VOL_or_BND: VorB = VorB.VOL) -> FESpaceElementRange: 
        """
        Returns an iterable range of elements.

        Parameters:

        VOL_or_BND : ngsolve.comp.VorB
          input VOL, BND, BBND,...
        """
    def FinalizeUpdate(self) -> None: 
        """
        finalize update
        """
    def FreeDofs(self, coupling: bool = False) -> pyngcore.BitArray: 
        """
        Return BitArray of free (non-Dirichlet) dofs\n
        coupling=False ... all free dofs including local dofs\n
        coupling=True .... only element-boundary free dofs

        Parameters:

        coupling : bool
          input coupling
        """
    @overload
    def GetDofNrs(self, ni: NodeId) -> tuple: 
        """
        Parameters:

        ei : ngsolve.comp.ElementId
          input element id





        Parameters:

        ni : ngsolve.comp.NodeId
          input node id
        """
    @overload
    def GetDofNrs(self, ei: ElementId) -> tuple: ...
    def GetDofs(self, region: Region) -> pyngcore.BitArray: 
        """
        Returns all degrees of freedom in given region.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    def GetFE(self, ei: ElementId) -> object: 
        """
        Get the finite element to corresponding element id.

        Parameters:

        ei : ngsolve.comp.ElementId
           input element id
        """
    def GetOrder(self, nodeid: NodeId) -> int: 
        """
        return order of node.
        by now, only isotropic order is supported here
        """
    def GetTrace(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def GetTraceTrans(self, arg0: FESpace, arg1: ngsolve.la.BaseVector, arg2: ngsolve.la.BaseVector, arg3: bool) -> None: ...
    def HideAllDofs(self, component: object = <ngsolve.ngstd.DummyArgument>) -> None: 
        """
        set all visible coupling types to HIDDEN_DOFs (will be overwritten by any Update())
        """
    def InvM(self, rho: ngsolve.fem.CoefficientFunction = None) -> ngsolve.la.BaseMatrix: ...
    def Mass(self, rho: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None) -> ngsolve.la.BaseMatrix: ...
    def ParallelDofs(self) -> ngsolve.la.ParallelDofs: 
        """
        Return dof-identification for MPI-distributed meshes
        """
    def Prolongation(self) -> ngmg::Prolongation: 
        """
        Return prolongation operator for use in multi-grid
        """
    def Range(self, component: int) -> ngsolve.ngstd.IntRange: 
        """
                 Return interval of dofs of a component of a product space.

        Parameters:

        component : int
          input component
        """
    @overload
    def SetCouplingType(self, dofnrs: ngsolve.ngstd.IntRange, coupling_type: COUPLING_TYPE) -> None: 
        """
                 Set coupling type of a degree of freedom.

        Parameters:

        dofnr : int
          input dof number

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type




                 Set coupling type for interval of dofs.

        Parameters:

        dofnrs : Range
          range of dofs

        coupling_type : ngsolve.comp.COUPLING_TYPE
          input coupling type
        """
    @overload
    def SetCouplingType(self, dofnr: int, coupling_type: COUPLING_TYPE) -> None: ...
    def SetDefinedOn(self, region: Region) -> None: 
        """
        Set the regions on which the FESpace is defined.

        Parameters:

        region : ngsolve.comp.Region
          input region
        """
    @overload
    def SetOrder(self, nodeid: NodeId, order: int) -> None: 
        """
        Parameters:

        element_type : ngsolve.fem.ET
          input element type

        order : object
          input polynomial order



        Parameters:

        nodeid : ngsolve.comp.NodeId
          input node id

        order : int
          input polynomial order
        """
    @overload
    def SetOrder(self, element_type: ngsolve.fem.ET, order: int) -> None: ...
    def SolveM(self, vec: ngsolve.la.BaseVector, rho: ngsolve.fem.CoefficientFunction = None, definedon: Region = None) -> None: 
        """
                 Solve with the mass-matrix. Available only for L2-like spaces.

        Parameters:

        vec : ngsolve.la.BaseVector
          input right hand side vector

        rho : ngsolve.fem.CoefficientFunction
          input CF
        """
    def TestFunction(self) -> object: 
        """
        Return a proxy to be used as a testfunction for Symbolic Integrators
        """
    def TnT(self) -> Tuple[object, object]: 
        """
        Return a tuple of trial and testfunction
        """
    def TraceOperator(self, tracespace: FESpace, average: bool) -> ngsolve.la.BaseMatrix: ...
    def TrialFunction(self) -> object: 
        """
        Return a proxy to be used as a trialfunction in Symbolic Integrators
        """
    def Update(self) -> None: 
        """
        update space after mesh-refinement
        """
    def UpdateDofTables(self) -> None: 
        """
        update dof-tables after changing polynomial order distribution
        """
    def __eq__(self, space: FESpace) -> bool: ...
    @staticmethod
    def __flags_doc__() -> dict: ...
    def __getstate__(self) -> tuple: ...
    def __init__(self, mesh: Mesh, autoupdate: bool = False, **kwargs) -> None: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    @staticmethod
    def __special_treated_flags__() -> dict: ...
    def __str__(self) -> str: ...
    def __timing__(self) -> object: ...
    @property
    def __memory__(self) -> List[Tuple[str, int, int]]:
        """
        :type: List[Tuple[str, int, int]]
        """
    @property
    def components(self) -> tuple:
        """
        Return a list of the components of a product space

        :type: tuple
        """
    @property
    def couplingtype(self) -> FlatArray_enum ngcomp::COUPLING_TYPE_S:
        """
        :type: FlatArray_enum ngcomp::COUPLING_TYPE_S
        """
    @property
    def dim(self) -> int:
        """
        multi-dim of FESpace

        :type: int
        """
    @property
    def globalorder(self) -> int:
        """
        query global order of space

        :type: int
        """
    @property
    def is_complex(self) -> bool:
        """
        :type: bool
        """
    @property
    def loembedding(self) -> ngsolve.la.BaseMatrix:
        """
        :type: ngsolve.la.BaseMatrix
        """
    @property
    def lospace(self) -> FESpace:
        """
        :type: FESpace
        """
    @property
    def mesh(self) -> Mesh:
        """
        mesh on which the FESpace is created

        :type: Mesh
        """
    @property
    def name(self) -> str:
        """
        :type: str
        """
    @name.setter
    def name(self, arg1: str) -> None:
        pass
    @property
    def ndof(self) -> int:
        """
        number of degrees of freedom

        :type: int
        """
    @property
    def ndofglobal(self) -> int:
        """
        global number of dofs on MPI-distributed mesh

        :type: int
        """
    @property
    def type(self) -> str:
        """
        type of finite element space

        :type: str
        """
    __dict__: mappingproxy # value = mappingproxy({'__init__': <instancemethod __init__ at 0x0000026DEE416678>, '__dict__': <attribute '__dict__' of 'ngsolve.comp.VectorSurfaceL2' objects>, '__doc__': "\n\n\n Keyword arguments can be:\norder: int = 1\n  order of finite element space\ncomplex: bool = False\n  Set if FESpace should be complex\ndirichlet: regexpr\n  Regular expression string defining the dirichlet boundary.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet = 'top|right'\ndirichlet_bbnd: regexpr\n  Regular expression string defining the dirichlet bboundary,\n  i.e. points in 2D and edges in 3D.\n  More than one boundary can be combined by the | operator,\n  i.e.: dirichlet_bbnd = 'top|right'\ndefinedon: Region or regexpr\n  FESpace is only defined on specific Region, created with mesh.Materials('regexpr')\n  or mesh.Boundaries('regexpr'). If given a regexpr, the region is assumed to be\n  mesh.Materials('regexpr').\ndim: int = 1\n  Create multi dimensional FESpace (i.e. [H1]^3)\ndgjumps: bool = False\n  Enable discontinuous space for DG methods, this flag is needed for DG methods,\n  since the dofs have a different coupling then and this changes the sparsity\n  pattern of matrices.\nlow_order_space: bool = True\n  Generate a lowest order space together with the high-order space,\n  needed for some preconditioners.\norder_policy: ORDER_POLICY = ORDER_POLICY.OLDSTYLE\n  CONSTANT .. use the same fixed order for all elements,\n  NODAL ..... use the same order for nodes of same shape,\n  VARIBLE ... use an individual order for each edge, face and cell,\n  OLDSTYLE .. as it used to be for the last decade\n", '__module__': 'ngsolve.comp', '__getstate__': <instancemethod __getstate__ at 0x0000026DEE4166D8>, '__setstate__': <instancemethod __setstate__ at 0x0000026DEE416738>, '__flags_doc__': <staticmethod object at 0x0000026DEE4149C8>})
    pass
class VorB():
    """
    Enum specifying the codimension. VOL is volume, BND is boundary and BBND is codimension 2 (edges in 3D, points in 2D)

    Members:

      VOL

      BND

      BBND

      BBBND
    """
    def __init__(self, arg0: int) -> None: ...
    def __int__(self) -> int: ...
    @property
    def name(self) -> str:
        """
        (self: handle) -> str

        :type: str
        """
    BBBND: ngsolve.comp.VorB # value = VorB.BBBND
    BBND: ngsolve.comp.VorB # value = VorB.BBND
    BND: ngsolve.comp.VorB # value = VorB.BND
    VOL: ngsolve.comp.VorB # value = VorB.VOL
    __entries: dict # value = {'VOL': (VorB.VOL, None), 'BND': (VorB.BND, None), 'BBND': (VorB.BBND, None), 'BBBND': (VorB.BBBND, None)}
    __members__: dict # value = {'VOL': VorB.VOL, 'BND': VorB.BND, 'BBND': VorB.BBND, 'BBBND': VorB.BBBND}
    pass
def BndElementId(nr: int) -> ElementId:
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
def CompressCompound(fespace: FESpace, active_dofs: object = <ngsolve.ngstd.DummyArgument>) -> FESpace:
    pass
def ConvertOperator(spacea: FESpace, spaceb: FESpace, trial_proxy: ProxyFunction = None, trial_cf: ngsolve.fem.CoefficientFunction = None, definedon: Optional[Region] = None, vb: VorB = VorB.VOL, range_dofs: pyngcore.BitArray = None, localop: bool = False, parmat: bool = True, use_simd: bool = True, bonus_intorder_ab: int = 0, bonus_intorder_bb: int = 0) -> ngsolve.la.BaseMatrix:
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
@overload
def Integrate(igls: SumOfIntegrals, mesh: Mesh, element_wise: bool = False) -> object:
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
def Integrate(cf: ngsolve.fem.CoefficientFunction, mesh: Union[Mesh, Region], VOL_or_BND: VorB = VorB.VOL, order: int = 5, definedon: Region = None, region_wise: bool = False, element_wise: bool = False) -> object:
    pass
def Interpolate(cf: ngsolve.fem.CoefficientFunction, space: FESpace, bonus_intorder: int = 0) -> ngsolve.fem.CoefficientFunction:
    """
    Interpolate a CoefficientFunction into the finite element space.
    The interpolation is canonical interpolation using dual shapes.
    The result is a CoefficientFunction.
    Interpolation is done on the fly for each element, no global GridFunction is allocated.
    """
def KSpaceCoeffs(arg0: GridFunction, arg1: GridFunction, arg2: float, arg3: float) -> None:
    pass
def MPI_Init() -> netgen.libngpy._meshing.MPI_Comm:
    pass
def PatchwiseSolve(bf: SumOfIntegrals, lf: SumOfIntegrals, gf: GridFunction) -> None:
    pass
def Prolongate(arg0: GridFunction, arg1: GridFunction) -> None:
    pass
def ProlongateCoefficientFunction(arg0: ngsolve.fem.CoefficientFunction, arg1: int, arg2: FESpace) -> ngsolve.fem.CoefficientFunction:
    pass
def SetHeapSize(size: int) -> None:
    """
    Set a new heapsize.

    Parameters:

    size : int
      input heap size
    """
def SetTestoutFile(file: str) -> None:
    """
    Enable some logging into file with given filename

    Parameters:

    file : string
      input file name
    """
def SymbolicBFI(form: ngsolve.fem.CoefficientFunction, VOL_or_BND: VorB = VorB.VOL, element_boundary: bool = False, skeleton: bool = False, definedon: Optional[Union[Region, list]] = None, intrule: ngsolve.fem.IntegrationRule = <ngsolve.fem.IntegrationRule object at 0x0000026DEE4255F0>, bonus_intorder: int = 0, definedonelements: pyngcore.BitArray = None, simd_evaluate: bool = True, element_vb: VorB = VorB.VOL, geom_free: bool = False, deformation: GridFunction = None) -> ngsolve.fem.BFI:
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
def SymbolicEnergy(form: ngsolve.fem.CoefficientFunction, VOL_or_BND: VorB = VorB.VOL, definedon: object = <ngsolve.ngstd.DummyArgument>, element_boundary: bool = False, bonus_intorder: int = 0, definedonelements: object = <ngsolve.ngstd.DummyArgument>, simd_evaluate: bool = True, element_vb: VorB = VorB.VOL, deformation: GridFunction = None) -> ngsolve.fem.BFI:
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
def SymbolicLFI(form: ngsolve.fem.CoefficientFunction, VOL_or_BND: VorB = VorB.VOL, element_boundary: bool = False, skeleton: bool = False, definedon: Optional[Union[Region, list]] = None, intrule: ngsolve.fem.IntegrationRule = <ngsolve.fem.IntegrationRule object at 0x0000026DEE4256F0>, bonus_intorder: int = 0, definedonelements: pyngcore.BitArray = None, simd_evaluate: bool = True, element_vb: VorB = VorB.VOL, deformation: GridFunction = None) -> ngsolve.fem.LFI:
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
def SymbolicTPBFI(form: ngsolve.fem.CoefficientFunction, VOL_or_BND: VorB = VorB.VOL, element_boundary: bool = False, skeleton: bool = False, definedon: object = <ngsolve.ngstd.DummyArgument>) -> ngsolve.fem.BFI:
    pass
def TensorProductFESpace(spaces: list, flags: pyngcore.Flags = <pyngcore.Flags object at 0x0000026DEE425BB0>) -> FESpace:
    pass
@overload
def TensorProductIntegrate(arg0: GridFunction, arg1: list, arg2: ngsolve.fem.CoefficientFunction) -> float:
    pass
@overload
def TensorProductIntegrate(gftp: GridFunction, gfx: GridFunction, weight: ngsolve.fem.CoefficientFunction = None) -> None:
    pass
@overload
def Transfer2StdMesh(arg0: ngsolve.fem.CoefficientFunction, arg1: GridFunction) -> None:
    pass
@overload
def Transfer2StdMesh(gftp: GridFunction, gfstd: GridFunction) -> None:
    pass
BBBND: ngsolve.comp.VorB # value = VorB.BBBND
BBND: ngsolve.comp.VorB # value = VorB.BBND
BND: ngsolve.comp.VorB # value = VorB.BND
VOL: ngsolve.comp.VorB # value = VorB.VOL
ngsglobals: ngsolve.comp.GlobalVariables # value = <ngsolve.comp.GlobalVariables object at 0x0000026DEE409B70>
