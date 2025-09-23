from typing import Annotated, overload

import numpy
from numpy.typing import NDArray


class HeightLayer:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, bottom: float, top: float, delta: float = 100.0) -> None: ...

    @property
    def bottom(self) -> float:
        """The bottom of the HeightLayer (meters)"""

    @bottom.setter
    def bottom(self, arg: float, /) -> None: ...

    @property
    def top(self) -> float:
        """The top of the HeightLayer (meters)"""

    @top.setter
    def top(self, arg: float, /) -> None: ...

    @property
    def delta(self) -> float:
        """The HeightLayer delta (increment) to use if iterating (meters)."""

    @delta.setter
    def delta(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

class PressureLayer:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, bottom: float, top: float, delta: float = -1000.0) -> None: ...

    @property
    def bottom(self) -> float:
        """The bottom of the PressureLayer (Pa)"""

    @bottom.setter
    def bottom(self, arg: float, /) -> None: ...

    @property
    def top(self) -> float:
        """The top of the PressureLayer (Pa)"""

    @top.setter
    def top(self, arg: float, /) -> None: ...

    @property
    def delta(self) -> float:
        """The PressureLayer delta (increment) to use if iterating (Pa)"""

    @delta.setter
    def delta(self, arg: float, /) -> None: ...

    def __repr__(self) -> str: ...

class LayerIndex:
    def __init__(self, kbot: int, ktop: int) -> None: ...

    @property
    def kbot(self) -> int:
        """
        The bottom index of a layer on a coordinate array (pressure or height).
        """

    @kbot.setter
    def kbot(self, arg: int, /) -> None: ...

    @property
    def ktop(self) -> int:
        """The top index of a layer on a coordinate array (pressure or height)."""

    @ktop.setter
    def ktop(self, arg: int, /) -> None: ...

    def __repr__(self) -> str: ...

@overload
def get_layer_index(layer: PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> LayerIndex:
    """
    Finds the interior array indices encapsulated by the given PressureLayer. 
    More specifically, the returned LayerIndex excludes the exact top and bottom 
    values (in mathematical notation, [bottom, top]). This behavior is due to the 
    fact many algorithms use interpolation to get the exact values of arbitrary
    top/bottom locations within a profile. 

    If layer.bottom or layer.top are out of bounds, this function 
    will truncate the layer to the coordinate range of data provided 
    by coord[] in an attempt to gracefully continue and produce a result. 
    This will modify the value of the given PressureLayer.

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressures

    Returns
    -------
    nwsspc.sharp.calc.layer.LayerIndex
        A LayerIndex with {kbot, ktop}.
    """

@overload
def get_layer_index(layer: HeightLayer, height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> LayerIndex:
    """
    Finds the interior array indices encapsulated by the given HeightLayer. 
    More specifically, the returned LayerIndex excludes the exact top and bottom 
    values (in mathematical notation, [bottom, top]). This behavior is due to the 
    fact many algorithms use interpolation to get the exact values of arbitrary
    top/bottom locations within a profile. 

    If layer.bottom or layer.top are out of bounds, this function 
    will truncate the layer to the coordinate range of data provided 
    by coord[] in an attempt to gracefully continue and produce a result. 
    This will modify the value of the given HeightLayer.

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.HeightLayer
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of heights

    Returns
    -------
    nwsspc.sharp.calc.layer.LayerIndex
        A LayerIndex with {kbot, ktop}.
    """

def height_layer_to_pressure(layer: HeightLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], isAGL: bool = False) -> PressureLayer:
    """
    Converts a HeightLayer to a PressureLayer via interpolation, with the
    optional argument to convert the HeightLayer from meters AGL to MSL by adding
    the station height to the HeightLayer. If for some strange reason 
    you provide a HeightLayer that is out of the bounds of height[], then
    the bottom and top of the output layer will be set to MISSING.

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.HeightLayer
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of heights
    isAGL : bool 
        Whether or not the station height needs to be added for interpolation (default: False) 

    Returns
    -------
    nwsspc.sharp.calc.layer.PressureLayer
    """

def pressure_layer_to_height(layer: PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], toAGL: bool = False) -> HeightLayer:
    """
    Converts a PressureLayer to a HeightLayer via interpolation, with the 
    optional argument to convert the HeightLayer to meters AGL by subtracting off 
    the station height from the returned HeightLayer. If for some strange reason
    you provide a PressureLayer that is out of the bounds of pressure[], then 
    the bottom and top of the output layer will be set to MISSING. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure 
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of heights
    toAGL : bool
        Whether or not to subtract the station height from the HeightLayer (default: False)

    Returns 
    -------
    nwsspc.sharp.calc.layer.HeightLayer
    """

@overload
def layer_min(layer: HeightLayer, height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[float, float]:
    """
    Returns the minimum value of the data array within the given HeightLayer. The 
    function bounds checks the layer by calling get_layer_index. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.HeightLayer 
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of heights 
    data : numpy.ndarray[dtype=float32]
        1D array of data 

    Returns
    -------
    tuple[float, float]
        (min_value, level_of_min)
    """

@overload
def layer_min(layer: PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[float, float]:
    """
    Returns the minimum value of the data array within the given PressureLayer. The 
    function bounds checks the layer by calling get_layer_index. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer 
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of pressure 
    data : numpy.ndarray[dtype=float32] 
        1D array of data 

    Returns
    -------
    tuple[float, float]
        (min_value, level_of_min)
    """

@overload
def layer_max(layer: HeightLayer, height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[float, float]:
    """
    Returns the maximum value observed within the given data array over
    the given HeightLayer. The function bounds checks the layer by calling 
    get_layer_index. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.HeightLayer 
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values 
    data : numpy.ndarray[dtype=float32] 
        1D NumPy array of data values

    Returns
    -------
    tuple[float, float]
        (max_value, level_of_max)
    """

@overload
def layer_max(layer: PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[float, float]:
    """
    Returns the maximum value observed within the given data array over the given 
    PressureLayer. The function bounds checks the layer by calling get_layer_index. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer 
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values 
    data : numpy.ndarray[dtype=float32] 
        1D NumPy array of data values

    Returns
    -------
    tuple[float, float]
        (max_value, level_of_max)
    """

@overload
def layer_mean(PressureLayer: PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Computes the pressure-weighted mean value of a field over 
    a given PressureLayer. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer 
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of pressure values 
    data : numpy.ndarray[dtype=float32] 
        1D NumPy array of data values

    Returns
    -------
    float
    """

@overload
def layer_mean(HeightLayer: HeightLayer, height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], isAGL: bool = False) -> float:
    """
    Computes the pressure-weighted mean value of a field over 
    a given HeightLayer. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.HeightLayer 
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values 
    data : numpy.ndarray[dtype=float32]
        1D NumPy array of data values
    isAGL : bool
        Whether or not the surface station height should be added to the HeightLayer (default: False)

    Returns
    -------
    float
    """

@overload
def integrate_layer_trapz(layer: HeightLayer, data: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], integ_sign: int = 0, weighted: bool = False) -> float:
    """
    Returns a trapezoidal integration of the given data array over 
    the given HeightLayer. There is an additional argument that 
    determines whether this is a weighted average or not. The sign 
    of the integration may be passed as well, i.e. integrating only 
    positive or negative area, by passing a 1, 0, or -1 to integ_sign. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.HeightLayer 
    data : numpy.ndarray[dtype=float32] 
        1D NumPy array of data values 
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height values 
    integ_sign : int 
        The sign of the area to integrate (-1: negative; 1: positive; 0: both; default: 0)
    weighted : bool 
        Boolean determining whether or not the integration is weighted by the coordinate array 

    Returns
    -------
    float
    """

@overload
def integrate_layer_trapz(layer: PressureLayer, data: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], integ_sign: int = 0, weighted: bool = False) -> float:
    """
    Returns the trapezoidal integration of the given data array over 
    the given PressureLayer. There is an additional argument that 
    determines whether this is a weighted average or not. The sign 
    of the integration may be passed as well, i.e. integrating only 
    positive or negative area, by passing a 1, 0, or -1 to integ_sign. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer
    data : numpy.ndarray[dtype=float32] 
        1D NumPy array of data values
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of pressure values 
    integ_sign : int 
        The sign of the area to integrate (-1: negative; 1: positive; 0: both; default: 0)
    weighted : bool 
        Boolean determining whether or not the integration is weighted by the coordinate array

    Returns 
    -------
    float
    """
