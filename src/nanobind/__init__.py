import sys
from . import _calc

try:
    from ._version import __version__
except ImportError:
    __version__ = "0.0.0" 

__doc__ = _calc.__doc__

__all__ = [s for s in dir(_calc) if not s.startswith('_')]

def _patch_submodules():
    """
    A helper function to patch nanobind submodules for Sphinx.
    """
    for name in __all__:
        if hasattr(_calc, name):
            submodule = getattr(_calc, name)
            
            public_path = f"{__name__}.{name}"
            
            submodule.__name__ = public_path
            
            sys.modules[public_path] = submodule

_patch_submodules()
del _patch_submodules 
del sys


