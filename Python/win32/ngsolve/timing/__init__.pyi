import ngsolve.timing
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import i
import k
import c
import p
import e
import l
class TaskManager():
    def __enter__(self) -> None: ...
    def __exit__(self, arg0: object, arg1: object, arg2: object) -> None: ...
    @overload
    def __init__(self) -> None: 
        """
        Run paje-tracer, specify buffersize in bytes
        """
    @overload
    def __init__(self, pajetrace: int) -> None: ...
    @staticmethod
    def __timing__() -> List[Tuple[str, float]]: ...
    pass
class Timing():
    """
    Class for timing analysis of performance critical functions. Some 
    classes export a C++ function as __timing__, which returns a map 
    of performance critical parts with their timings. The class can save 
    these maps, load them and compare them. It can be saved as a benchmark 
    to be compared against.

    2 overloaded __init__ functions:

    1. __init__(name,obj,parallel=True,serial=True)
    2. __init__(filename)

    Parameters
    ----------

    name (str): Name for the timed class (for output formatting and 
        saving/loading of results)
    obj (NGSolve object): Some NGSolve class which has the __timing__ 
        functionality implemented. Currently supported classes:
            FESpace
    filename (str): Filename to load a previously saved Timing
    parallel (bool=True): Time in parallel (using TaskManager)
    serial (bool=True): Time not in parallel (not using TaskManager)
    """
    class type():
        pass
    __dict__: mappingproxy # value = mappingproxy({'__module__': 'ngsolve.timing', '__doc__': '\nClass for timing analysis of performance critical functions. Some \nclasses export a C++ function as __timing__, which returns a map \nof performance critical parts with their timings. The class can save \nthese maps, load them and compare them. It can be saved as a benchmark \nto be compared against.\n\n2 overloaded __init__ functions:\n\n1. __init__(name,obj,parallel=True,serial=True)\n2. __init__(filename)\n\nParameters\n----------\n\nname (str): Name for the timed class (for output formatting and \n    saving/loading of results)\nobj (NGSolve object): Some NGSolve class which has the __timing__ \n    functionality implemented. Currently supported classes:\n        FESpace\nfilename (str): Filename to load a previously saved Timing\nparallel (bool=True): Time in parallel (using TaskManager)\nserial (bool=True): Time not in parallel (not using TaskManager)\n\n', '__init__': <function Timing.__init__ at 0x0000026DF1E7D9D8>, '__str__': <function Timing.__str__ at 0x0000026DF1E7DA68>, 'Save': <function Timing.Save at 0x0000026DF1E7DAF8>, 'CompareTo': <function Timing.CompareTo at 0x0000026DF1E7DB88>, 'CompareToBenchmark': <function Timing.CompareToBenchmark at 0x0000026DF1E7DC18>, 'SaveBenchmark': <function Timing.SaveBenchmark at 0x0000026DF1E7DCA8>, '__dict__': <attribute '__dict__' of 'Timing' objects>, '__weakref__': <attribute '__weakref__' of 'Timing' objects>})
    __weakref__: getset_descriptor # value = <attribute '__weakref__' of 'Timing' objects>
    pass
__all__ = ['Timing']
