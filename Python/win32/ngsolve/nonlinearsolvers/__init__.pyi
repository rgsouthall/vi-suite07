import ngsolve.nonlinearsolvers
from typing import *
from typing import Iterable as iterable
from typing import Iterator as iterator
from numpy import float64
_Shape = Tuple[int, ...]
__all__  = [
"NewtonSolver",
"Projector",
"InnerProduct",
"Newton",
"NewtonMinimization",
"Norm",
"TimeFunction",
"sqrt"
]
class NewtonSolver():
    class type():
        """
        type(object_or_name, bases, dict)
        type(object) -> the object's type
        type(name, bases, dict) -> a new type
        """
        class object():
            """
            The most base type
            """
            class type():
                pass
            pass
        class type():
            pass
        @staticmethod
        def __prepare__() -> dict: ...
        __abstractmethods__: getset_descriptor # value = <attribute '__abstractmethods__' of 'type' objects>
        __bases__: tuple # value = (<class 'object'>,)
        __basicsize__ = 864
        __dict__: mappingproxy # value = mappingproxy({'__repr__': <slot wrapper '__repr__' of 'type' objects>, '__call__': <slot wrapper '__call__' of 'type' objects>, '__getattribute__': <slot wrapper '__getattribute__' of 'type' objects>, '__setattr__': <slot wrapper '__setattr__' of 'type' objects>, '__delattr__': <slot wrapper '__delattr__' of 'type' objects>, '__init__': <slot wrapper '__init__' of 'type' objects>, '__new__': <built-in method __new__ of type object at 0x00007FFB865159A0>, 'mro': <method 'mro' of 'type' objects>, '__subclasses__': <method '__subclasses__' of 'type' objects>, '__prepare__': <method '__prepare__' of 'type' objects>, '__instancecheck__': <method '__instancecheck__' of 'type' objects>, '__subclasscheck__': <method '__subclasscheck__' of 'type' objects>, '__dir__': <method '__dir__' of 'type' objects>, '__sizeof__': <method '__sizeof__' of 'type' objects>, '__basicsize__': <member '__basicsize__' of 'type' objects>, '__itemsize__': <member '__itemsize__' of 'type' objects>, '__flags__': <member '__flags__' of 'type' objects>, '__weakrefoffset__': <member '__weakrefoffset__' of 'type' objects>, '__base__': <member '__base__' of 'type' objects>, '__dictoffset__': <member '__dictoffset__' of 'type' objects>, '__mro__': <member '__mro__' of 'type' objects>, '__name__': <attribute '__name__' of 'type' objects>, '__qualname__': <attribute '__qualname__' of 'type' objects>, '__bases__': <attribute '__bases__' of 'type' objects>, '__module__': <attribute '__module__' of 'type' objects>, '__abstractmethods__': <attribute '__abstractmethods__' of 'type' objects>, '__dict__': <attribute '__dict__' of 'type' objects>, '__doc__': <attribute '__doc__' of 'type' objects>, '__text_signature__': <attribute '__text_signature__' of 'type' objects>})
        __dictoffset__ = 264
        __flags__ = 2148291584
        __itemsize__ = 40
        __mro__: tuple # value = (<class 'type'>, <class 'object'>)
        __name__ = 'type'
        __text_signature__: NoneType # value = None
        __weakrefoffset__ = 368
        pass
    __dict__: mappingproxy # value = mappingproxy({'__module__': 'ngsolve.nonlinearsolvers', '__init__': <function NewtonSolver.__init__ at 0x0000026DF1CFEE58>, 'Solve': <function TimeFunction.<locals>.retfunc at 0x0000026DF1CFEF78>, 'SetDirichlet': <function NewtonSolver.SetDirichlet at 0x0000026DF1E7D048>, '_UpdateInverse': <function NewtonSolver._UpdateInverse at 0x0000026DF1E7D0D8>, '__dict__': <attribute '__dict__' of 'NewtonSolver' objects>, '__weakref__': <attribute '__weakref__' of 'NewtonSolver' objects>, '__doc__': None})
    __weakref__: getset_descriptor # value = <attribute '__weakref__' of 'NewtonSolver' objects>
    pass
class Projector():
    pass
def InnerProduct(x: object, y: object, **kwargs) -> object:
    """
    Computes InnerProduct of given objects
    """
def Norm(x: object) -> object:
    """
    Compute Norm
    """
