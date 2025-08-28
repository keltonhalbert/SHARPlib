import sys, inspect
from . import _calc

try:
    from ._version import __version__
except ImportError:
    __version__ = "0.0.0" 

__doc__ = _calc.__doc__

__all__ = [s for s in dir(_calc) if not s.startswith('_')]

def _patch_submodules():
    """
    A helper function to patch pybind11 submodules for Sphinx.
    """
    
    def _recursive_patch(module_obj, public_path_prefix):
        """Recursively update the __module__ attribute of nested objects."""
        for name, obj in inspect.getmembers(module_obj):
            if name.startswith("__"):
                continue

            if hasattr(obj, "__module__") and isinstance(obj.__module__, str) and obj.__module__.startswith(_calc.__name__):
                obj.__module__ = public_path_prefix


    for name in __all__:
        if hasattr(_calc, name):
            submodule = getattr(_calc, name)
            public_path = f"{__name__}.{name}"
            submodule.__name__ = public_path
            sys.modules[public_path] = submodule
            _recursive_patch(submodule, public_path)

_patch_submodules()
del _patch_submodules 
del sys


