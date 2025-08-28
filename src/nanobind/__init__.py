from . import _calc

try:
    from ._version import __version__
except ImportError:
    __version__ = "0.0.0" 

__doc__ = _calc.__doc__

__all__ = [s for s in dir(_calc) if not s.startswith('_')]

