import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray

import nwsspc.sharp.calc.layer
import nwsspc.sharp.calc.parcel


class adiabat(enum.Enum):
    pseudo_liq = 1
    """Use a pseudoadiabatic parcel with liquid only processes"""

    adiab_liq = 2
    """Use an adiabatic parcel with liquid only processes"""

    pseudo_ice = 3
    """Use a pseudoadiabatic parcel with liquid and ice processes"""

    adiab_ice = 4
    """Use an adiabatic parcel with liquid and ice processes"""

@overload
def wobf(temperature: float) -> float:
    """
    Computes the difference between the wet-bulb potential temperatures
    (Kelvin) for saturated and dry air, given the temperature (Kelvin).

    The Wobus Function (wobf) is defined as the difference between the 
    wet-bulb potential temperature for saturated air (WBPTS) and the 
    wet-bulb potential temperature for dry air (WBPTD) given an iput
    air temperature. This function is used in the computation of moist 
    adiabatic ascent as part of the Wobus parcel lifter (lifter_wobus).

    WOBF(T) := WBPTS - WBPTD

    Although WBPTS and WBPTD are functions of both pressure and 
    temperature, it is assumed their difference is a function of 
    temperature only. The difference is also proportional to the 
    heat imparted to a parcel of air.

    This function uses a polynomial approximation to the wobus function, 
    fitted to value in Table 78 of PP.319-322 of the Smithsonian Meteorological
    Table by Rolan List (6th Revised Edition). Herman Wobus, a mathematician 
    for the Navy Weather Research Facility in Norfolk, VA computed these 
    coefficients a very long time ago, as he was retired as of the time of 
    the documentation found on this routine written in 1981.

    It was shown by Robert Davies-Jones (2007) that the Wobus function has
    a slight dependence on pressure, which results in errors of up to 1.2
    Kelvin in the temperature of a lifted parcel. 

    Parameters
    ----------
    temperature : float 
        The air temperature (Kelvin)

    Returns
    -------
    float
        The Wobus function temperature (Kelvin)
    """

@overload
def wobf(tmpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Computes the difference between the wet-bulb potential temperatures
    (Kelvin) for saturated and dry air, given the temperature (Kelvin).

    The Wobus Function (wobf) is defined as the difference between the 
    wet-bulb potential temperature for saturated air (WBPTS) and the 
    wet-bulb potential temperature for dry air (WBPTD) given an iput
    air temperature. This function is used in the computation of moist 
    adiabatic ascent as part of the Wobus parcel lifter (lifter_wobus).

    WOBF(T) := WBPTS - WBPTD

    Although WBPTS and WBPTD are functions of both pressure and 
    temperature, it is assumed their difference is a function of 
    temperature only. The difference is also proportional to the 
    heat imparted to a parcel of air.

    This function uses a polynomial approximation to the wobus function, 
    fitted to value in Table 78 of PP.319-322 of the Smithsonian Meteorological
    Table by Rolan List (6th Revised Edition). Herman Wobus, a mathematician 
    for the Navy Weather Research Facility in Norfolk, VA computed these 
    coefficients a very long time ago, as he was retired as of the time of 
    the documentation found on this routine written in 1981.

    It was shown by Robert Davies-Jones (2007) that the Wobus function has
    a slight dependence on pressure, which results in errors of up to 1.2
    Kelvin in the temperature of a lifted parcel. 

    Parameters
    ----------
    tmpk_arr : numpy.ndarray[dtype=float32] 
        The 1D array of air temperatures (Kelvin)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of Wobus function temperatures (Kelvin)
    """

@overload
def lcl_temperature(temperature: float, dewpoint: float) -> float:
    """
    Compute the Lifted Condensation Level (LCL) temperature (Kelvin) 
    given an air temperature (Kelvin) and dewpoint temperature (Kelvin).

    The LCL temperature is computed as in Bolton (1980) eq 15, and is
    considered to be within a 10th of a degree of the more exact 
    iterative formula.

    Parameters
    ----------
    temperature : float 
        The air temperature (Kelvin)
    dewpoint : float 
        The dewpoint temperature (Kelvin)

    Returns
    -------
    float
        The LCL temperature (Kelvin)
    """

@overload
def lcl_temperature(tmpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dwpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the Lifted Condensation Level (LCL) temperature (Kelvin) 
    given an array of air temperatures (Kelvin) and dewpoint 
    temperatures (Kelvin).

    The LCL temperatures are computed as in Bolton (1980) eq 15, and is
    considered to be within a 10th of a degree of the more exact 
    iterative formula.

    Parameters
    ----------
    tmpk_arr : numpy.ndarray[dtype=float32]
        The 1D air temperature array (Kelvin)
    dwpk_arr : numpy.ndarray[dtype=float32]
        The 1D dewpoint temperature array (Kelvin)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of LCL temperatures (Kelvin)
    """

@overload
def vapor_pressure(pressure: float, temperature: float) -> float:
    """
    Compute the vapor pressure with respect to liquid water.

    The vapor pressure (or saturation vapor pressure) is computed with 
    respect to liquid water when the dewpoint temperature (or air temperature)
    is passed (in Kelvin), along with the air pressure (Pa). The air pressure 
    is used as a minimum floor value for extremely cold temperatures at low
    pressures, and is consistent with how vapor pressure is treated in CM1. 

    This function uses the formulation by Bolton (1980), and is
    accurate to within 0.3% for the temperature range of -35C <= T <= 35C.

    Parameters
    ----------
    pressure : float 
        The air pressure (Pa)
    temperature : float 
        The air temperature (K) or dewpoint temperature (K)

    Returns
    -------
    float
        The vapor pressure (Pa) given dewpoint temperature, or the 
        saturation vapor pressure given the air temperature.
    """

@overload
def vapor_pressure(pres_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], tmpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the vapor pressure with respect to liquid water.

    The vapor pressure (or saturation vapor pressure) is computed with 
    respect to liquid water when the dewpoint temperature (or air temperature)
    is passed (in Kelvin), along with the air pressure (Pa). The air pressure 
    is used as a minimum floor value for extremely cold temperatures at low
    pressures, and is consistent with how vapor pressure is treated in CM1. 

    This function uses the formulation by Bolton (1980), and is
    accurate to within 0.3% for the temperature range of -35C <= T <= 35C.

    Parameters
    ----------
    pres_arr : numpy.ndarray[dtype=float32] 
        The 1D array of air pressure (Pa)
    tmpk_arr : numpy.ndarray[dtype=float32]
        The 1D array of air temperature (K) or dewpoint temperature (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of vapor pressure (Pa) given dewpoint temperature, or the 
        saturation vapor pressure given the air temperature.
    """

@overload
def vapor_pressure_ice(pressure: float, temperature: float) -> float:
    """
    Compute the vapor pressure with respect to ice.

    The vapor pressure (or saturation vapor pressure) is computed with 
    respect to ice when the dewpoint temperature (or air temperature)
    is passed (in Kelvin), along with the air pressure (Pa). The air pressure 
    is used as a minimum floor value for extremely cold temperatures at low
    pressures, and is consistent with how vapor pressure is treated in CM1. 

    This function uses the formulation by Bolton (1980), and is
    accurate to within 0.3% for the temperature range of -35C <= T <= 35C.

    Parameters
    ----------
    pressure : float 
        The air pressure (Pa)
    temperature : float 
        The air temperature (K) or dewpoint temperature (K)

    Returns
    -------
    float
        The vapor pressure (Pa) given dewpoint temperature, or the 
        saturation vapor pressure given the air temperature.
    """

@overload
def vapor_pressure_ice(pres_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], tmpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the vapor pressure with respect to ice.

    The vapor pressure (or saturation vapor pressure) is computed with 
    respect to ice when the dewpoint temperature (or air temperature)
    is passed (in Kelvin), along with the air pressure (Pa). The air pressure 
    is used as a minimum floor value for extremely cold temperatures at low
    pressures, and is consistent with how vapor pressure is treated in CM1. 

    This function uses the formulation by Bolton (1980), and is
    accurate to within 0.3% for the temperature range of -35C <= T <= 35C.

    Parameters
    ----------
    pres_arr : numpy.ndarray[dtype=float32] 
        The 1D array of air pressure (Pa)
    tmpk_arr : numpy.ndarray[dtype=float32]
        The 1D array of air temperature (K) or dewpoint temperature (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of vapor pressure (Pa) given dewpoint temperature, or the 
        saturation vapor pressure given the air temperature.
    """

@overload
def temperature_at_mixratio(wv_mixratio: float, pressure: float) -> float:
    """
    Computes the temperature (K) of air at the given water vapor mixing ratio
    (kg/kg) and air pressure (Pa). Can be used to compute the dewpoint temperature 
    from mixing ratio. 

    The routine is implemented as in Bolton (1980) eq 11, and is considered to be 
    accurate to 0.03 K for -35C <= T <= 35C. 

    Parameters
    ----------
    wv_mixratio : float 
        The water vapor mixing ratio (kg/kg)
    pressure : float 
        The air pressure (Pa)

    Returns
    -------
    float
        The temperature (K) of an air parcel at a given mixing ratio and pressure.
    """

@overload
def temperature_at_mixratio(mixr_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], pres_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Computes the temperature (K) of air at the given water vapor mixing ratio
    (kg/kg) and air pressure (Pa). Can be used to compute the dewpoint temperature 
    from mixing ratio. 

    The routine is implemented as in Bolton (1980) eq 11, and is considered to be 
    accurate to 0.03 K for -35C <= T <= 35C. 

    Parameters
    ----------
    mixr_arr : numpy.ndarray[dtype=float32] 
        The 1D array of water vapor mixing ratio (kg/kg)
    pres_arr : numpy.ndarray[dtype=float32]
        The 1D array of air pressure (Pa)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of temperature (K) of an air parcel at a given mixing ratio and pressure.
    """

@overload
def theta_level(potential_temperature: float, temperature: float) -> float:
    """
    Computes the pressure level (Pa) of a parcel given the potential temperature (K) and air 
    temperature (K).

    Parameters
    ----------
    potential_temperature : float 
        The potential temperature, or theta (K)
    temperature : float 
        The air temperature (K)

    Returns
    -------
    float
        The pressure level (Pa) corresponding to the potential temperature and air temperature.
    """

@overload
def theta_level(theta_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], tmpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Computes the pressure level (Pa) of a parcel given the potential temperature (K) and air 
    temperature (K).

    Parameters
    ----------
    theta_arr : numpy.ndarray[dtype=float32] 
        The 1D array of potential temperature, or theta (K)
    tmpk_arr : numpy.ndarray[dtype=float32]
        The 1D array of air temperature (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D pressure level (Pa) corresponding to the potential temperature and air temperature.
    """

@overload
def theta(pressure: float, temperature: float, ref_pressure: float = 100000.0) -> float:
    """
    Computes the potential temperature (K), or theta, given the air pressure (Pa), air temperature (K),
    and a reference pressure (default value is 100000 Pa).

    Parameters
    ----------
    pressure : float 
        The air pressure (Pa)
    temperature : float 
        The air temperature (K)
    ref_pressure : float, optional 
        The reference pressure (default 100000.0 Pa)

    Returns
    -------
    float
        The potential temperature (K), or theta
    """

@overload
def theta(pres_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], tmpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], ref_pressure: float = 100000.0) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Computes the potential temperature (K), or theta, given the air pressure (Pa), air temperature (K),
    and a reference pressure (default value is 100000 Pa).

    Parameters
    ----------
    pressure : numpy.ndarray[dtype=float32]
        The 1D array of air pressure (Pa)
    temperature : numpy.ndarray[dtype=float32]
        The 1D array of air temperature (K)
    ref_pressure : float, optional
        The reference pressure (default 100000.0 Pa)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of potential temperature (K), or theta
    """

@overload
def mixratio(q: float) -> float:
    """
    Compute the water vapor mixing ratio (kg/kg) from specific humidity (kg/kg).

    Parameters
    ----------
    q : float 
        The specific humidity (kg/kg)

    Returns
    -------
    float
        The water vapor mixing ratio (kg/kg)
    """

@overload
def mixratio(pressure: float, temperature: float) -> float:
    """
    Compute the water vapor mixing ratio (kg/kg) from the air pressure (Pa) and temperature (K).
    If given the air temperature, this is the saturation mixing ratio. If given the dewpoint 
    temperature, tis is the mixing ratio. 

    Parameters
    ----------
    pressure : float
        The air pressure (Pa)
    temperature : float 
        The air temperature (K)

    Returns
    -------
    float
        The water vapor mixing ratio (kg/kg)
    """

@overload
def mixratio(spfh_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the water vapor mixing ratio (kg/kg) from the specific humidity (kg/kg).

    Parameters
    ----------
    spfh_arr : numpy.ndarray[dtype=float32] 
        The 1D array of specific humidity (kg/kg)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of water vapor mixing ratio (kg/kg)
    """

@overload
def mixratio(pres_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], tmpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the water vapor mixing ratio (kg/kg) from the air pressure (Pa) and temperature (K).
    If given the air temperature, this is the saturation mixing ratio. If given the dewpoint 
    temperature, this is the mixing ratio. 

    Parameters
    ----------
    pres_arr : numpy.ndarray[dtype=float32] 
        The 1D array of air pressure (Pa)
    tmpk_arr : numpy.ndarray[dtype=float32]
        The 1D array of air temperature (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of water vapor mixing ratio (kg/kg)
    """

@overload
def mixratio_ice(pressure: float, temperature: float) -> float:
    """
    Compute the ice water mixing ratio (kg/kg) from the air pressure (Pa) and temperature (K).
    If given the air temperatuer, this is the saturation mixing ratio. If given the dewpoint 
    temperature, this is the mixing ratio. 

    Parameters
    ----------
    pressure : float 
        The air pressure (Pa)
    temperature : float 
        The air temperature (K)

    Returns
    -------
    float
        The ice water mixing ratio (kg/kg)
    """

@overload
def mixratio_ice(pres_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], tmpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the ice water mixing ratio (kg/kg) from the air pressure (Pa) and temperature (K).
    If given the air temperatuer, this is the saturation mixing ratio. If given the dewpoint 
    temperature, this is the mixing ratio. 

    Parameters
    ----------
    pres_arr : numpy.ndarray[dtype=float32]
        The 1D array of air pressure (Pa)
    tmpl_arr : numpy.ndarray[dtype=float32] 
        The 1D array of air temperature (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The ice water mixing ratio (kg/kg)
    """

@overload
def specific_humidity(rv: float) -> float:
    """
    Compute the specific humidity (kg/kg) from a mixing ratio (kg/kg).

    Parameters
    ----------
    rv : float 
        The water vapor mixing ratio (kg/kg)

    Returns
    -------
    float
        The specific humidity (kg/kg)
    """

@overload
def specific_humidity(mixr_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the specific humidity (kg/kg) from a mixing ratio (kg/kg).

    Parameters
    ----------
    mixr_arr : numpy.ndarray[dtype=float32]
        The 1D array of water vapor mixing ratios (kg/kg)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of specific humidity (unitless)
    """

@overload
def relative_humidity(pressure: float, temperature: float, dewpoint: float) -> float:
    """
    Compute the relative humidity of water vapor with respect to liquid water.
    NOTE: The pressure variable is only used as a sanity check when computing 
    vapor pressure at extremely low pressures and temperatures. If you do not 
    want or need this behavior, you can pass the THETA_REF_PRESSURE constant
    in place of an air pressure.

    Parameters
    ----------
    pressure : float 
        the ambient air pressure (Pa)
    temperature : float 
        the ambient air temperature (K)
    dewpoint : float 
        the ambient dewpoint temperature (K)

    Returns
    -------
    float
        Relative Humidity (fraction, unitless)
    """

@overload
def relative_humidity(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the relative humidity of water vapor with respect to liquid water.
    NOTE: The pressure variable is only used as a sanity check when computing 
    vapor pressure at extremely low pressures and temperatures. If you do not
    want or need this behavior, you canpass in an array of THETA_REF_PRESSURE
    in place of air pressure.

    Parameters
    ----------
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient air pressure Pa)
    temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of ambient air temperature (K)
    dewpoint : numpy.ndarray[dtype=float32]
        1D NumPy array of ambient dewpoint temperature (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of relative humidity values (fraction, unitless)
    """

@overload
def virtual_temperature(temperature: float, rv: float, rl: float = 0.0, ri: float = 0.0) -> float:
    """
    Returns the virtual temperature in Kelvin given the dry-bulb 
    temperature (Kelvin), the water vapor mixing ratio (kg/kg), the 
    liquid water mixing ratio (kg/kg), and the ice water mixing ratios
    (kg/kg). The liquid and ice water mixing ratios have default values 
    of zero, if unspecified. 

    Parameters
    ----------
    temperature : float 
        The dry-bulb temperature (K)
    rv : float 
        The water vapor mixing ratio (kg/kg)
    rl : float 
        The liquid water mixing ratio (kg/kg) 
    ri : float 
        The ice water mixing ratio (kg/kg)

    Returns
    -------
    float
        The virtual temperature (K)
    """

@overload
def virtual_temperature(tmpk_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], rv_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], rl_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)] | None = None, ri_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)] | None = None) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Returns the virtual temperature in Kelvin given the dry-bulb 
    temperature (Kelvin), the water vapor mixing ratio (kg/kg), the 
    liquid water mixing ratio (kg/kg), and the ice water mixing ratios
    (kg/kg). The liquid and ice water mixing ratios have default values 
    of zero, if unspecified. 

    Parameters
    ----------
    tmpk_arr : numpy.ndarray[dtype=float32] 
        The 1D array of dry-bulb temperature (K)
    rv_arr : numpy.ndarray[dtype=float32] 
        The 1D array of water vapor mixing ratio (kg/kg)
    rl_arr : numpy.ndarray[dtype=float32], optional 
        The 1D array of liquid water mixing ratio (kg/kg) 
    ri_arr : numpy.ndarray[dtype=float32], optional 
        The 1D array of ice water mixing ratio (kg/kg)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of virtual temperature (K)
    """

@overload
def density(pressure: float, temperature: float, wv_mixratio: float = 0.0) -> float:
    """
    Compute the density of air.

    Computes the air density given pressure and temperature.

    Parameters 
    ----------
    pressure : float 
        The air pressure (Pa)
    temperature : float 
        The air temperature (K)
    wv_mixratio : float, optional 
        The water vapor mixing ratio (kg/kg). Default is 0.

    Returns 
    -------
    float 
        Air density (kg m^-3)
    """

@overload
def density(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixr_arr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)] | None = None) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]: ...

@overload
def wetlift(pressure: float, temperature: float, lifted_pressure: float, converge: float = 0.0010000000474974513) -> float:
    """
    Compute the temperature of a parcel lifted moist adiabatically to a new level. 

    With a given parcel defined by a pressure (Pa) and temperature (K), lift it 
    moist adiabatically to a new pressure level (Pa) and return the temperature of 
    the parcel at that level. 

    This function relies on the Wobus Function (thermo.wobf), and it was shown by 
    Robert Davies-Jones (2007) that the WObus function has a slight dependence on 
    pressure, which results in errors of up to 1.2 K in the temperature of a lifted 
    parcel. 

    Parameters
    ----------
    pressure : float 
        The air pressure (Pa)
    temperature : float 
        The saturated air temperature (K)
    lifted_pressure : float 
        The new pressure level to lift to (Pa)
    converge : float 
        The iterative convergence criteria (K; default = 0.001)

    Returns
    -------
    float
        The new temperature (K) when lifted moist adiabatically to the new pressure level
    """

@overload
def wetlift(arg0: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], arg1: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], arg2: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], /) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the temperature of a parcel lifted moist adiabatically to a new level. 

    With a given parcel defined by a pressure (Pa) and temperature (K), lift it 
    moist adiabatically to a new pressure level (Pa) and return the temperature of 
    the parcel at that level. 

    This function relies on the Wobus Function (thermo.wobf), and it was shown by 
    Robert Davies-Jones (2007) that the WObus function has a slight dependence on 
    pressure, which results in errors of up to 1.2 K in the temperature of a lifted 
    parcel. 

    Parameters
    ----------
    pressure : numpy.ndarray[dtype=float32] 
        The 1D array of air pressures (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        The 1D array of saturated air temperatures (K)
    lifted_pressure : numpy.ndarray[dtype=float32]
        The 1D array of new pressure levels to lift to (Pa)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        The 1D array of new temperatures (K) when lifted moist adiabatically to the new pressure levels
    """

def drylift(pressure: float, temperature: float, dewpoint: float) -> tuple[float, float]:
    """
    Given a parcel's initial pressure (Pa), temperature (K), and 
    dewpoint temperature (K), lift the parcel dry adiabaitically to its
    Lifted Condensation Level and return the resulting LCL pressure (Pa)
    and temperature (K).

    The LCL temperature is computed using an approximation. See the 
    lcl_temperature function documentation for more information.

    Parameters
    ----------
    pressure : float 
        Parcel starting pressure (Pa)
    temperature : float 
        Parcel starting temperature (K)
    dewpoint : float 
        Parcel starting dewpoint temperature (K)

    Returns
    -------
    tuple[float, float]
        A tuple of (lcl_pres, lcl_temperature)
    """

@overload
def wetbulb(lifter: nwsspc.sharp.calc.parcel.lifter_wobus, pressure: float, temperature: float, dewpoint: float) -> float:
    """
    Compute the wet bulb temperature (K) given the ambient pressure (Pa), 
    temperature (K), and dewpoint temperature (K).

    First, it lifts a parcel with the given pressure, temperature, and
    dewpoint temperature to its Lifted Condensation Level (LCL). To compute the 
    temperature and pressure of the LCL, an approximation is used. See the 
    lcl_temperature function for further detail. 

    After the parcel has reached the LCL, the lifter passed to the function 
    lowers the parcel to its initial pressure level along a moist adiabat or 
    pseudoadiabat. 

    Parameters
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_wobus 
        a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure : float 
        The ambient pressure (Pa)
    temperature : float 
        The ambient temperature (K)
    dewpoint : float 
        The ambient dewpoint temperature (K)

    Returns
    -------
    float
        The wetbulb temperature (K)
    """

@overload
def wetbulb(lifter: nwsspc.sharp.calc.parcel.lifter_cm1, pressure: float, temperature: float, dewpoint: float) -> float:
    """
    Compute the wet bulb temperature (K) given the ambient pressure (Pa), 
    temperature (K), and dewpoint temperature (K).

    First, it lifts a parcel with the given pressure, temperature, and
    dewpoint temperature to its Lifted Condensation Level (LCL). To compute the 
    temperature and pressure of the LCL, an approximation is used. See the 
    lcl_temperature function for further detail. 

    After the parcel has reached the LCL, the lifter passed to the function 
    lowers the parcel to its initial pressure level along a moist adiabat or 
    pseudoadiabat. 

    Parameters
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_cm1 
        a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure : float 
        The ambient pressure (Pa)
    temperature : float 
        The ambient temperature (K)
    dewpoint : float 
        The ambient dewpoint temperature (K)

    Returns
    -------
    float
        The wetbulb temperature (K)
    """

@overload
def wetbulb(lifter: nwsspc.sharp.calc.parcel.lifter_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the wet bulb temperature (K) given the ambient pressure (Pa), 
    temperature (K), and dewpoint temperature (K).

    First, it lifts a parcel with the given pressure, temperature, and
    dewpoint temperature to its Lifted Condensation Level (LCL). To compute the 
    temperature and pressure of the LCL, an approximation is used. See the 
    lcl_temperature function for further detail. 

    After the parcel has reached the LCL, the lifter passed to the function 
    lowers the parcel to its initial pressure level along a moist adiabat or 
    pseudoadiabat. 

    Parameters
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_wobus 
        a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient pressures (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient temperatureds (K)
    dewpoint : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient dewpoint temperatures (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of wetbulb temperatures (K)
    """

@overload
def wetbulb(lifter: nwsspc.sharp.calc.parcel.lifter_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the wet bulb temperature (K) given the ambient pressure (Pa), 
    temperature (K), and dewpoint temperature (K).

    First, it lifts a parcel with the given pressure, temperature, and
    dewpoint temperature to its Lifted Condensation Level (LCL). To compute the 
    temperature and pressure of the LCL, an approximation is used. See the 
    lcl_temperature function for further detail. 

    After the parcel has reached the LCL, the lifter passed to the function 
    lowers the parcel to its initial pressure level along a moist adiabat or 
    pseudoadiabat. 

    Parameters
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_cm1 
        a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient pressures (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient temperatureds (K)
    dewpoint : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient dewpoint temperatures (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of wetbulb temperatures (K)
    """

@overload
def theta_wetbulb(lifter: nwsspc.sharp.calc.parcel.lifter_wobus, pressure: float, temperature: float, dewpoint: float) -> float:
    """
    Compute the wet-bulb potential temperature (K) given the pressure (Pa), 
    temperature (K), and dewpoint temperature (K).

    First, it lifts a parcel with the given pressure, temperature, 
    and dewpoint temperature to its Lifted Condensation Level (LCL). 
    To compute the temperature and pressure of the LCL, an approximation 
    is used. See the lcl_temperature fuction for further detail. 

    After the parcel has reached the LCL, the lifted passed to the function 
    lowers the parcel to the standard parcel reference pressure level 
    (1000 hPa) along a moist adiabat.

    Parameters
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_wobus 
        a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure : float 
        the ambient air pressure (Pa)
    temperature : float 
        the ambient air temperature (K)
    dewpoint : float 
        the ambient dewpoint temperature

    Returns
    -------
    float
        The wet-bulb potential temperature (K)
    """

@overload
def theta_wetbulb(lifter: nwsspc.sharp.calc.parcel.lifter_cm1, pressure: float, temperature: float, dewpoint: float) -> float:
    """
    Compute the wet-bulb potential temperature (K) given the pressure (Pa), 
    temperature (K), and dewpoint temperature (K).

    First, it lifts a parcel with the given pressure, temperature, 
    and dewpoint temperature to its Lifted Condensation Level (LCL). 
    To compute the temperature and pressure of the LCL, an approximation 
    is used. See the lcl_temperature fuction for further detail. 

    After the parcel has reached the LCL, the lifted passed to the function 
    lowers the parcel to the standard parcel reference pressure level 
    (1000 hPa) along a moist adiabat.

    Parameters
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_cm1 
        a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure : float 
        the ambient air pressure (Pa)
    temperature : float 
        the ambient air temperature (K)
    dewpoint : float 
        the ambient dewpoint temperature

    Returns
    -------
    float
        The wet-bulb potential temperature (K)
    """

@overload
def theta_wetbulb(lifter: nwsspc.sharp.calc.parcel.lifter_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the wet-bulb potential temperature (K) given the pressure (Pa), 
    temperature (K), and dewpoint temperature (K).

    First, it lifts a parcel with the given pressure, temperature, 
    and dewpoint temperature to its Lifted Condensation Level (LCL). 
    To compute the temperature and pressure of the LCL, an approximation 
    is used. See the lcl_temperature fuction for further detail. 

    After the parcel has reached the LCL, the lifted passed to the function 
    lowers the parcel to the standard parcel reference pressure level 
    (1000 hPa) along a moist adiabat.

    Parameters
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_wobus 
        a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient air pressure (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient air temperature (K)
    dewpoint : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient dewpoint temperature

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of wet-bulb potential temperatures (K)
    """

@overload
def theta_wetbulb(lifter: nwsspc.sharp.calc.parcel.lifter_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the wet-bulb potential temperature (K) given the pressure (Pa), 
    temperature (K), and dewpoint temperature (K).

    First, it lifts a parcel with the given pressure, temperature, 
    and dewpoint temperature to its Lifted Condensation Level (LCL). 
    To compute the temperature and pressure of the LCL, an approximation 
    is used. See the lcl_temperature fuction for further detail. 

    After the parcel has reached the LCL, the lifted passed to the function 
    lowers the parcel to the standard parcel reference pressure level 
    (1000 hPa) along a moist adiabat.

    Parameters
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_cm1 
        a parcel lifter (e.g. lifter_cm1 or lifter_wobus)
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient air pressure (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient air temperature (K)
    dewpoint : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient dewpoint temperature

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of wet-bulb potential temperatures (K)
    """

@overload
def thetae(pressure: float, temperature: float, dewpoint: float) -> float:
    """
    Compute the equivalent potential temperature (K) given the 
    air pressure (Pa), air temperature (K), and dewpoint temperature (K).
    The equivalent potential temperature is computed as in Bolton 1980. 
    It was found in Davies-Jones 2009 to be the most accurate non-iterative 
    formulation of theta-e.

    References
    ----------
    Bolton 1980: https://doi.org/10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2

    Davies-Jones 2009: https://doi.org/10.1175/2009MWR2774.1

    Parameters
    ----------
    pressure : float 
        the ambient air pressure (Pa)
    temperature : float 
        the ambient air temperature (K)
    dewpoint : float 
        the dewpoint temperature (K)

    Returns
    -------
    float
        The equivalent potential temperature (K)
    """

@overload
def thetae(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the equivalent potential temperature (K) given the 
    air pressure (Pa), air temperature (K), and dewpoint temperature (K).
    The equivalent potential temperature is computed as in Bolton 1980. 
    It was found in Davies-Jones 2009 to be the most accurate non-iterative 
    formulation of theta-e.

    References 
    ----------
    Bolton 1980: https://doi.org/10.1175/1520-0493(1980)108%3C1046:TCOEPT%3E2.0.CO;2

    Davies-Jones 2009: https://doi.org/10.1175/2009MWR2774.1

    Parameters
    ----------
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient air pressure (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of ambient air temperature (K)
    dewpoint : numpy.ndarray[dtype=float32] 
        1D NumPy array of dewpoint temperature (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of equivalent potential temperature (K)
    """

@overload
def lapse_rate(layer_agl: nwsspc.sharp.calc.layer.HeightLayer, height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Computes the lapse rate over a given HeightLayer (meters AGL).
    This routine handles converting the height AGL to MSL by adding
    the surface height value to the HeightLayer.

    Parameters
    ----------
    layer_agl : nwsspc.sharp.calc.layer.HeightLayer 
        a HeightLayer (meters AGL) to compute the lapse rate over
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values (meters)
    temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of temperature values (K)

    Returns
    -------
    float
        The temperature lapse rate (K)
    """

@overload
def lapse_rate(layer: nwsspc.sharp.calc.layer.PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Computes the lapse rate over a given PressureLayer (Pa).

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer 
        a PressureLayer (Pa) to compute the lapse rate over
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values (meters)
    temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of temperature values (K)

    Returns
    -------
    float
        The temperature lapse rate (K)
    """

@overload
def lapse_rate_max(layer: nwsspc.sharp.calc.layer.HeightLayer, depth: float, height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[float, nwsspc.sharp.calc.layer.HeightLayer]:
    """
    Given a layer of the atmosphere (e.g. 2 - 6 km), find the maximum
    lapse rate over the provided depth (e.g. 2 km) within that given layer. 
    Returns the maximum lapse rate, as well as the layer it was found in. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.HeightLayer 
        a HeightLayer to search for the maximum lapse lapse_rate (meters AGL) 
    depth : float 
        the depth used to compute a lapse rate within a layer (meters)
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values (meters)
    temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of temperature values (K)

    Returns
    -------
    tuple[float, nwsspc.sharp.calc.layer.HeightLayer]
        A tuple containing the maximum lapse rate and the layer it was found in
    """

@overload
def lapse_rate_max(layer: nwsspc.sharp.calc.layer.PressureLayer, depth: float, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[float, nwsspc.sharp.calc.layer.PressureLayer]:
    """
    Given a layer of the atmosphere (e.g. 800 hPa - 500 hPa), find the maximum
    lapse rate over the provided depth (e.g. 100 hPa) within that given layer. 
    Returns the maximum lapse rate, as well as the layer it was found in. 

    Parameters
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer 
        a PressureLayer to search for the maximum lapse lapse_rate (Pa) 
    depth : float 
        the depth used to compute a lapse rate within a layer (Pa)
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values (meters)
    temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of temperature values (K)

    Returns
    -------
    tuple[float, nwsspc.sharp.calc.layer.PressureLayer]
        A tuple containing the maximum lapse rate and the layer it was found in
    """

@overload
def buoyancy(parcel_temperature: float, environment_temperature: float) -> float:
    """
    Compute buoyancy given parcel & environment temperatures.

    Parameters
    ----------
    parcel_temperature : float 
        parcel temperature (K)
    environment_temperature : float 
        environment temperature (K)

    Returns
    -------
    float
        Buoyancy (m/s^2)
    """

@overload
def buoyancy(parcel_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], environment_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute buoyancy given arrays of parcel & environment temperatures.

    Parameters
    ----------
    parcel_temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of parcel temperatures (K)
    environment_temperature : numpy.ndarray[dtype=float32] 
        1D NumPy array of environment temperatures (K)

    Returns
    -------
    numpy.ndarray[dtype=float32]
        1D NumPy array of buoyancy values (m/s^2)
    """
