from typing import Annotated

import numpy
from numpy.typing import NDArray


def interp_height(hght_val: float, hght_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Interpolate a value from an array in height coordinates (meters).
    The coordinate array (hght_arr) is assumed to be sorted (ascending).

    Parameters
    ----------
    hght_val : float 
        The coordinate height value to interpolate to (meters)
    hght_arr : numpy.ndarray[dtype=float32]
        1D numpy array of height values to interpolate from (meters)
    data_arr : numpy.ndarray[dtype=float32]
        float 1D numpy array of data values to interpolate from 

    Returns
    -------
    float 
        An interpolated data value from data_arr corresponding to hght_val
    """

def interp_pressure(pres_val: float, pres_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Interpolate a value from an array in pressure coordinates (Pa).
    All pressure interpolation happens in log10 space.
    The coordinate array (pres_arr) is assumed to be sorted (descending).

    Parameters
    ----------
    pres_val : float 
        The coordinate pressure value to interpolate to (Pa)
    pres_arr : nump.ndarray[dtype=float32] 
        1D numpy array of pressure values to interpolate from (Pa)
    data_arr : numpy.ndarray[dtype=float32]
        1D numpy array of data values to interpolate from 

    Returns
    -------
    float 
        An interpolated data value from data_arr corresponding to pres_val
    """

def find_first_pressure(data_val: float, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data_array: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Conducts a bottom-up search for the first occurrence of a given value,
    and interpolates in order to get the pressure level it occurs at.

    Parameters
    ----------
    data_val : float 
        The value to search for 
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    data_array : numpy.ndarray[dtype=float32]
        1D NumPy array of values to search over

    Returns
    -------
    float
        The pressure level of first occurrence (Pa)
    """

def find_first_height(data_val: float, height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], data_array: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Conducts a bottom-up search for the first occurrence of a given value,
    and interpolates in order to get the height level it occurs at.

    Parameters
    ----------
    data_val : float 
        The value to search for 
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height values (meters)
    data_array : numpy.ndarray[dtype=float32]
        1D NumPy array of values to search over

    Returns
    -------
    float
        The height level of first occurrence (meters)
    """
