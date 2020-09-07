"""pybind ngstd"""
import ngsolve.ngstd
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
import pyngcore
__all__  = [
"Archive",
"DummyArgument",
"HeapReset",
"IntRange",
"LocalHeap",
"Timer",
"Tracer",
"_MemoryView",
"TestFlagsConversion",
"Timers",
"_PickleMemory",
"_UnpickleMemory"
]
class Archive():
    def __and__(self, array: pyngcore.Array_I_S) -> Archive: ...
    def __init__(self, filename: str, write: bool, binary: bool) -> None: ...
    pass
class DummyArgument():
    def __bool__(self) -> bool: ...
    def __repr__(self) -> str: ...
    pass
class HeapReset():
    """
    stores heap-pointer on init, and resets it on exit
    """
    def __init__(self, lh: LocalHeap) -> None: ...
    pass
class IntRange():
    def __init__(self, arg0: int, arg1: int) -> None: ...
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
class LocalHeap():
    """
    A heap for fast memory allocation
    """
    def __init__(self, size: int = 1000000, name: str = 'PyLocalHeap') -> None: ...
    pass
class Timer():
    def Start(self) -> None: 
        """
        start timer
        """
    def Stop(self) -> None: 
        """
        stop timer
        """
    def __enter__(self) -> None: ...
    def __exit__(self, arg0: object, arg1: object, arg2: object) -> None: ...
    def __init__(self, arg0: str) -> None: ...
    @property
    def time(self) -> float:
        """
        returns time

        :type: float
        """
    pass
class Tracer():
    def SetMaxTracefileSize(self) -> None: ...
    def SetTraceThreadCounter(self) -> None: ...
    def SetTraceThreads(self) -> None: ...
    pass
class _MemoryView():
    def __getstate__(self) -> tuple: ...
    def __setstate__(self, arg0: tuple) -> None: ...
    pass
def TestFlagsConversion(flags: pyngcore.Flags) -> None:
    pass
def Timers() -> list:
    """
    Returns list of timers
    """
def _PickleMemory(pickler: object, view: ngstd::MemoryView) -> None:
    pass
def _UnpickleMemory(unpickler: object) -> None:
    pass
