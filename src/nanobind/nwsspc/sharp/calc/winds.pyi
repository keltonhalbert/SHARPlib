from typing import Annotated, overload

import numpy
from numpy.typing import NDArray

import nwsspc.sharp.calc.layer


class WindVector:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, speed: float, direction: float) -> None: ...

    @property
    def speed(self) -> float:
        """Wind Speed (m/s)"""

    @speed.setter
    def speed(self, arg: float, /) -> None: ...

    @property
    def direction(self) -> float:
        """Wind Direction (degrees from North)"""

    @direction.setter
    def direction(self, arg: float, /) -> None: ...

class WindComponents:
    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, u_comp: float, v_comp: float) -> None: ...

    @property
    def u(self) -> float:
        """U wind component (m/s)"""

    @u.setter
    def u(self, arg: float, /) -> None: ...

    @property
    def v(self) -> float:
        """V wind component (m/s)"""

    @v.setter
    def v(self, arg: float, /) -> None: ...

@overload
def helicity(layer: nwsspc.sharp.calc.layer.HeightLayer, storm_motion: WindComponents, height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Computes the Storm Relative Helicity (SRH) over a given layer using 
    storm motion vector components stored in a WindComponents object.

    This integration occurs over the given arrays, using interpolation 
    for the top and bottom of the integration layer, and native data levels
    in between. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.HeightLayer 
        A HeightLayer to integrate over (meters AGL)
    storm_motion : nwsspc.sharp.calc.winds.WindComponents 
        A WindComponents object with the storm motion in m/s 
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values (meters)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of environment U-wind components (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of environment V-wind components (m/s)

    Returns
    -------
    float
        Storm Relative Helicity (m^2/s^2)
    """

@overload
def helicity(layer: nwsspc.sharp.calc.layer.PressureLayer, storm_motion: WindComponents, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Computes the Storm Relative Helicity (SRH) over a given layer using 
    storm motion vector components stored in a WindComponents object.

    This integration occurs over the given arrays, using interpolation 
    for the top and bottom of the integration layer, and native data levels
    in between. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer 
        A PressureLayer to integrate over (Pa)
    storm_motion : nwsspc.sharp.calc.winds.WindComponents 
        A WindComponents object with the storm motion in m/s 
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values (meters)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of environment U-wind components (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of environment V-wind components (m/s)

    Returns
    -------
    float
        Storm Relative Helicity (m^2/s^2)
    """

@overload
def wind_shear(layer: nwsspc.sharp.calc.layer.HeightLayer, height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> WindComponents:
    """
    Computes the U and V components of the wind shear over a 
    layer given the vertical sounding arrays of height, u_wind,
    and v_wind.

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.HeightLayer 
        HeightLayer for which to compute wind shear 
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values (meters)
    u_wind : numpy.ndarray[dtype=float32] 
        1D NumPy array of U-wind components (m/s)
    v_wind : numpy.ndarray[dtype=float32] 
        1D NumPy array of V-wind components (m/s)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindComponents
        WindComponents of U and V wind shear components (m/s)
    """

@overload
def wind_shear(layer: nwsspc.sharp.calc.layer.PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> WindComponents:
    """
    Computes the U and V components of the wind shear over a 
    layer given the vertical sounding arrays of height, u_wind,
    and v_wind.

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer 
        PressureLayer for which to compute wind shear 
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values (meters)
    u_wind : numpy.ndarray[dtype=float32] 
        1D NumPy array of U-wind components (m/s)
    v_wind : numpy.ndarray[dtype=float32] 
        1D NumPy array of V-wind components (m/s)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindComponents
        WindComponents of U and V wind shear components (m/s)
    """

def mean_wind(layer: nwsspc.sharp.calc.layer.PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu')], weighted: bool = False) -> WindComponents:
    """
    Computes the mean wind over the given PressureLayer and input profile
    arrays of pressure, U-wind, and V-wind components. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer 
        PressureLayer over which to compute mean
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure coordinate values (Pa)
    u_wind : numpy.ndarray[dtype=float32] 
        1D NumPy array of U-wind component values (m/s)
    v_wind : numpy.ndarray[dtype=float32] 
        1D NumPy array of V-wind component values (m/s)
    weighted : bool 
        Boolean flag to compute pressure-weighted mean wind (default: False)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindComponents
        WindComponents of U anf V mean wind components (m/s)
    """

@overload
def vector_magnitude(u_comp: float, v_comp: float) -> float:
    """
    Given the zonal (U) and meridional (V) wind components of a vector,
    compute and return the magnitude (m/s) of the vector.

    Parameters
    ----------
    u_comp : float
        U-wind component (m/s)
    v_comp : float
        V-wind component (m/s)

    Returns
    -------
    float
        Wind speed (m/s)
    """

@overload
def vector_magnitude(u_comp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_comp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Given the zonal (U) and meridional (V) components of a vector, 
    compute and return the magnitude (m/s) of the vector. 

    Parameters
    ----------
    u_comp : numpy.ndarray[dtype=float32] 
        1D NumPy array of U-wind component (m/s)
    v_comp : numpy.ndarray[dtype=float32] 
        1D NumPy array of V-wind component (m/s)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of wind speed (m/s)
    """

@overload
def vector_angle(u_comp: float, v_comp: float) -> float:
    """
    Given the zonal (U) and meridional (V) components of a vector,
    compute and return the angle (from North) of the vector. 

    Parameters
    ----------
    u_comp : float 
        The U-wind component
    v_comp : float 
        The V-wind component

    Returns
    -------
    float
        Wind direction (degrees from North)
    """

@overload
def vector_angle(u_comp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_comp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Given the zonal (U) and meridional (V) components of a vector, 
    compute and return the angle (from North) of the vector. 

    Parameters
    ----------
    u_comp : numpy.ndarray[dtype=float32] 
        1D NumPy array of U-wind component (m/s)
    v_comp : numpy.ndarray[dtype=float32] 
        1D NumPy array of V-wind component (m/s)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of wind direction (degrees from North)
    """

@overload
def u_component(wind_speed: float, wind_direction: float) -> float:
    """
    Computes the zonal (U) wind component from a wind vector.

    Parameters
    ----------
    wind_speed : float 
        The vector speed (m/s)
    wind_direction : float 
        The vector direction (degrees from North)

    Returns
    -------
    float
        The U-wind component (m/s)
    """

@overload
def u_component(wind_speed: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], wind_direction: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Computes the zonal (U) wind component from a wind vector.

    Parameters
    ----------
    wind_speed : numpy.ndarray[dtype=float32] 
        1D NumPy array of wind speeds (m/s)
    wind_direction : numpy.ndarray[dtype=float32] 
        1D NumPy array of wind direction (degrees)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of U-wind component values (m/s)
    """

@overload
def v_component(wind_speed: float, wind_direction: float) -> float:
    """
    Computes the meridional (V) wind component from a wind vector.

    Parameters
    ----------
    wind_speed : float 
        Vector speed (m/s)
    wind_direction : float 
        Vector angle (degrees from North)

    Returns
    -------
    float
        The V-wind component
    """

@overload
def v_component(wind_speed: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], wind_direction: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Computes the meridional (V) wind component from a wind vector.

    Parameters
    ----------
    wind_speed : numpy.ndarray[dtype=float32] 
        1D NumPy array of wind speeds (m/s)
    wind_direction : numpy.ndarray[dtype=float32] 
        1D NumPy array of wind direction (degrees)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of V-wind component values (m/s)
    """

@overload
def components_to_vector(u_comp: float, v_comp: float) -> WindVector:
    """
    Given the zonal (U) and meridional (V) components of a vector, 
    compute and return the wind speed (m/s) and direction 
    (degrees from North) as a WindVector type.

    Parameters
    ----------
    u_comp : float 
        The U-wind component (m/s)
    v_comp : float 
        The V-wind component (m/s)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindVector
        WindVector containing wind speed (m/s) and direction (degrees from North)
    """

@overload
def components_to_vector(wind_comp: WindComponents) -> WindVector:
    """
    Given the components of a vector via a WindComponents object (m/s),
    compute and return the wind speed (m/s) and wind direction (degrees from North)
    as a WindVector type.

    Parameters
    ----------
    wind_comp : nwsspc.sharp.calc.winds.WindComponents 
        WindComponents (m/s)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindVector
        WindVector containing wind speed (m/s) and direction (degrees from North)
    """

@overload
def components_to_vector(u_comp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_comp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')], Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]]:
    """
    Given 1D NumPy arrays of zonal (U) and meridional (V) wind components,
    compute the wind speed (m/s) and direction (degrees from North).

    Parameters
    ----------
    u_comp : numpy.ndarray[dtype=float32] 
        1D NumPy array of U-wind component values (m/s)
    v_comp : numpy.ndarray[dtype=float32]
        1D NumPy array of V-wind component values (m/s)

    Returns
    -------
    tuple[numpy.ndarray[dtype=float32], numpy.ndarray[dtype=float32]]
        wspd: 1D NumPy array of wind speeds (m/s)

        wdir: 1D NumPy array of wind directiond (degrees from North)
    """

@overload
def vector_to_components(wind_speed: float, wind_direction: float) -> WindComponents:
    """
    Given the wind speed (m/s) and wind direction (degrees from North),
    compute and return the zonal and meridional vector components as WindComponents.

    Parameters
    ----------
    wind_speed : float 
        The magnitude/speed of the vector (m/s)
    wind_direction : float 
        The direction of the vector (degrees from North)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindComponents
        The U-wind and V-wind components (m/s) as WindComponents
    """

@overload
def vector_to_components(wind_vector: WindVector) -> WindComponents:
    """
    Given the wind speed (m/s) and direction (degrees from North),
    compute and return the zonal (U) and meridional (V) vector 
    components as WindComponents.

    Parameters
    ----------
    wind_vector : nwsspc.sharp.calc.winds.WindVector 
        WindVector containing wind wpeed (m/s) and direction (degrees from North)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindComponents
        U and V wind components (m/s) in a WindComponents object
    """

@overload
def vector_to_components(wspd: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], wdir: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')], Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]]:
    """
    Given 1D NumPy arrays of the wind speed (m/s) and direction (degrees from North),
    compute and return 1D NumPy arays of the zonal (U) and meridional (V) vector 
    componenWindComponents.

    Parameters
    ----------
    wspd : numpy.ndarray[dtype=float32] 
        1D NumPy array of wind speeds (m/s)
    wdir : numpy.ndarray[dtype=float32]
        1D NumPy array of wind directions (degrees from North)

    Returns
    -------
    tuple[numpy.ndarray[dtype=float32], numpy.ndarray[dtype=float32]]
        uwin: 1D NumPy array of U-wind components

        vwin: 1D NumPy array of V-wind components
    """
