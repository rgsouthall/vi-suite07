"""pybind solve"""
import ngsolve.solve
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import ngsolve.comp
import ngsolve.fem
import ngsolve.ngstd
__all__  = [
"BVP",
"CalcFlux",
"Draw",
"DrawFlux",
"SetVisualization",
"Tcl_Eval",
"_GetFacetValues",
"_GetValues",
"_GetVisualizationData",
"_SetLocale"
]
def BVP(bf: ngsolve.comp.BilinearForm, lf: ngsolve.comp.LinearForm, gf: ngsolve.comp.GridFunction, pre: ngsolve.comp.Preconditioner, maxsteps: int = 100, prec: float = 1e-08) -> ngsolve.comp.NumProc:
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
def CalcFlux(pde: ngsolve.comp.PDE, bf: ngsolve.comp.BilinearForm, gf: ngsolve.comp.GridFunction, flux: ngsolve.comp.GridFunction, applyd: bool = False) -> ngsolve.comp.NumProc:
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
@overload
def Draw(gf: ngsolve.comp.GridFunction, sd: int = 2, autoscale: bool = True, min: float = 0.0, max: float = 1.0, **kwargs) -> None:
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
def Draw(mesh: ngsolve.comp.Mesh, **kwargs) -> None:
    pass
@overload
def Draw(cf: ngsolve.fem.CoefficientFunction, mesh: ngsolve.comp.Mesh, name: str, sd: int = 2, autoscale: bool = True, min: float = 0.0, max: float = 1.0, draw_vol: bool = True, draw_surf: bool = True, reset: bool = False, **kwargs) -> None:
    pass
@overload
def Draw(arg0: object) -> None:
    pass
def DrawFlux(bf: ngsolve.comp.BilinearForm, gf: ngsolve.comp.GridFunction, label: str = 'flux', applyd: bool = False, useall: bool = False) -> ngsolve.comp.NumProc:
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
def SetVisualization(deformation: object = <ngsolve.ngstd.DummyArgument>, min: object = <ngsolve.ngstd.DummyArgument>, max: object = <ngsolve.ngstd.DummyArgument>, clipnormal: object = <ngsolve.ngstd.DummyArgument>, clipping: object = <ngsolve.ngstd.DummyArgument>) -> None:
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
def Tcl_Eval(arg0: str) -> None:
    pass
def _GetFacetValues(arg0: ngsolve.fem.CoefficientFunction, arg1: ngsolve.comp.Mesh, arg2: Dict[ngsolve.fem.ET, ngsolve.fem.IntegrationRule]) -> dict:
    pass
def _GetValues(arg0: ngsolve.fem.CoefficientFunction, arg1: ngsolve.comp.Mesh, arg2: ngsolve.comp.VorB, arg3: Dict[ngsolve.fem.ET, ngsolve.fem.IntegrationRule], arg4: bool) -> dict:
    pass
def _GetVisualizationData(arg0: ngsolve.comp.Mesh, arg1: Dict[ngsolve.fem.ET, ngsolve.fem.IntegrationRule]) -> dict:
    pass
def _SetLocale() -> None:
    pass
