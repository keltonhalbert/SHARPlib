from ._calc import (
    constants as constants,
    interp as interp,
    layer as layer,
    params as params,
    parcel as parcel,
    thermo as thermo,
    winds as winds,
)

try:
    from ._version import __version__
except ImportError:
    __version__ = "0.0.0" 

__all__ = [
    "constants",
    "interp",
    "layer",
    "params",
    "parcel",
    "thermo",
    "winds",
]

__doc__ = ._calc.__doc__

