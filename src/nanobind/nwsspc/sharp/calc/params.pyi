import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray

import nwsspc.sharp.calc.layer
import nwsspc.sharp.calc.parcel
import nwsspc.sharp.calc.winds


@overload
def effective_inflow_layer(lifter: nwsspc.sharp.calc.parcel.lifter_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtemp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], cape_thresh: float = 100.0, cinh_thresh: float = -250.0, mupcl: nwsspc.sharp.calc.parcel.Parcel | None = None) -> nwsspc.sharp.calc.layer.PressureLayer:
    """
    Computes the Effective Inflow Layer, or the layer of the atmosphere
    beliefed to be the primary source of inflow for supercell thunderstorms. 
    The Effective Inflow Layer, and its use in computing shear and storm 
    relative helicity, is described by Thompson et al. 2007: 
    https://www.spc.noaa.gov/publications/thompson/effective.pdf

    Standard/default values for cape_thresh and cinh_thresh have been 
    experimentally determined to be cape_thresh = 100 J/kg and 
    cinh_thresh = -250.0 J/kg. If an empty parcel object is passed via the 
    'mupcl' kwarg, the Most Unstable parcel found during the EIL search will 
    be returned. 

    References 
    ----------
    Thompson et al. 2007: https://www.spc.noaa.gov/publications/thompson/effective.pdf

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_wobus
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    dewpoint : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of dewpoint values (K)
    virtemp : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    cape_thresh : float, default = 100.0
        The CAPE threshold used to compute the Effective Inflow Layer 
    cinh_thresh : float, default = -250.0 
        The CINH threshold used to compute the Effective Inflow Layer
    mupcl : None or nwsspc.sharp.calc.parcel.Parcel, optional

    Returns 
    -------
    nwsspc.sharp.calc.layer.PressureLayer
        The Effective Inflow Layer
    """

@overload
def effective_inflow_layer(lifter: nwsspc.sharp.calc.parcel.lifter_lut_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtemp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], cape_thresh: float = 100.0, cinh_thresh: float = -250.0, mupcl: nwsspc.sharp.calc.parcel.Parcel | None = None) -> nwsspc.sharp.calc.layer.PressureLayer:
    """
    Computes the Effective Inflow Layer, or the layer of the atmosphere
    beliefed to be the primary source of inflow for supercell thunderstorms. 
    The Effective Inflow Layer, and its use in computing shear and storm 
    relative helicity, is described by Thompson et al. 2007: 
    https://www.spc.noaa.gov/publications/thompson/effective.pdf

    Standard/default values for cape_thresh and cinh_thresh have been 
    experimentally determined to be cape_thresh = 100 J/kg and 
    cinh_thresh = -250.0 J/kg. If an empty parcel object is passed via the 
    'mupcl' kwarg, the Most Unstable parcel found during the EIL search will 
    be returned. 

    References 
    ----------
    Thompson et al. 2007: https://www.spc.noaa.gov/publications/thompson/effective.pdf

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_lut_wobus
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    dewpoint : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of dewpoint values (K)
    virtemp : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    cape_thresh : float, default = 100.0
        The CAPE threshold used to compute the Effective Inflow Layer 
    cinh_thresh : float, default = -250.0 
        The CINH threshold used to compute the Effective Inflow Layer
    mupcl : None or nwsspc.sharp.calc.parcel.Parcel, optional

    Returns 
    -------
    nwsspc.sharp.calc.layer.PressureLayer
        The Effective Inflow Layer
    """

@overload
def effective_inflow_layer(lifter: nwsspc.sharp.calc.parcel.lifter_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtemp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], cape_thresh: float = 100.0, cinh_thresh: float = -250.0, mupcl: nwsspc.sharp.calc.parcel.Parcel | None = None) -> nwsspc.sharp.calc.layer.PressureLayer:
    """
    Computes the Effective Inflow Layer, or the layer of the atmosphere
    beliefed to be the primary source of inflow for supercell thunderstorms. 
    The Effective Inflow Layer, and its use in computing shear and storm 
    relative helicity, is described by Thompson et al. 2007: 
    https://www.spc.noaa.gov/publications/thompson/effective.pdf

    Standard/default values for cape_thresh and cinh_thresh have been 
    experimentally determined to be cape_thresh = 100 J/kg and 
    cinh_thresh = -250.0 J/kg. If an empty parcel object is passed via the 
    'mupcl' kwarg, the Most Unstable parcel found during the EIL search will 
    be returned. 

    References 
    ----------
    Thompson et al. 2007: https://www.spc.noaa.gov/publications/thompson/effective.pdf

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_cm1
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    dewpoint : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of dewpoint values (K)
    virtemp : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    cape_thresh : float, default = 100.0
        The CAPE threshold used to compute the Effective Inflow Layer 
    cinh_thresh : float, default = -250.0 
        The CINH threshold used to compute the Effective Inflow Layer
    mupcl : None or nwsspc.sharp.calc.parcel.Parcel, optional

    Returns 
    -------
    nwsspc.sharp.calc.layer.PressureLayer
        The Effective Inflow Layer
    """

@overload
def effective_inflow_layer(lifter: nwsspc.sharp.calc.parcel.lifter_lut_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtemp: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], cape_thresh: float = 100.0, cinh_thresh: float = -250.0, mupcl: nwsspc.sharp.calc.parcel.Parcel | None = None) -> nwsspc.sharp.calc.layer.PressureLayer:
    """
    Computes the Effective Inflow Layer, or the layer of the atmosphere
    beliefed to be the primary source of inflow for supercell thunderstorms. 
    The Effective Inflow Layer, and its use in computing shear and storm 
    relative helicity, is described by Thompson et al. 2007: 
    https://www.spc.noaa.gov/publications/thompson/effective.pdf

    Standard/default values for cape_thresh and cinh_thresh have been 
    experimentally determined to be cape_thresh = 100 J/kg and 
    cinh_thresh = -250.0 J/kg. If an empty parcel object is passed via the 
    'mupcl' kwarg, the Most Unstable parcel found during the EIL search will 
    be returned. 

    References 
    ----------
    Thompson et al. 2007: https://www.spc.noaa.gov/publications/thompson/effective.pdf

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_lut_cm1
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    dewpoint : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of dewpoint values (K)
    virtemp : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    cape_thresh : float, default = 100.0
        The CAPE threshold used to compute the Effective Inflow Layer 
    cinh_thresh : float, default = -250.0 
        The CINH threshold used to compute the Effective Inflow Layer
    mupcl : None or nwsspc.sharp.calc.parcel.Parcel, optional

    Returns 
    -------
    nwsspc.sharp.calc.layer.PressureLayer
        The Effective Inflow Layer
    """

@overload
def storm_motion_bunkers(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mean_wind_layer_agl: nwsspc.sharp.calc.layer.HeightLayer, wind_shear_layer_agl: nwsspc.sharp.calc.layer.HeightLayer, leftMover: bool = False, pressureWeighted: bool = False) -> nwsspc.sharp.calc.winds.WindComponents:
    """
    Estimates the supercell storm motion using the Bunkers et al. 2000 method 
    described in the following paper:
    https://doi.org/10.1175/1520-0434(2000)015%3C0061:PSMUAN%3E2.0.CO;2

    This does not use any of the updated methods described by Bunkers et al. 2014, 
    which uses Effective Inflow Layer metrics to get better estimates of storm 
    motion, especially when considering elevated convection. 

    In QC builds, returns MISSING components if either layer is MISSING, or if
    the mean wind layer or either 500 m deep end of the wind shear layer can't
    be converted to pressure (see
    nwsspc.sharp.calc.layer.height_layer_to_pressure), for example because it
    extends past the profile.

    References 
    ----------

    Buners et al. 2000: https://doi.org/10.1175/1520-0434(2000)015%3C0061:PSMUAN%3E2.0.CO;2

    Parameters 
    ----------
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        1D NumPy array of height values (meters)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of U wind component values (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of V wind compnent values (m/s)
    mean_wind_layer_agl : nwsspc.sharp.calc.layer.HeightLayer 
        HeightLayer (AGL) for computing the mean wind 
    wind_shear_layer_agl : nwsspc.sharp.calc.layer.HeightLayer 
        HeightLayer (AGL) for computing wind_shear
    leftMover : bool 
        Whether to compute left mover supercell motion (default: False)
    pressureWeighted : bool 
        Whether to use the pressure weighted mean wind (default: False)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindComponents
        U, V wind components of storm motion (m/s)
    """

@overload
def storm_motion_bunkers(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], eff_infl_lyr: nwsspc.sharp.calc.layer.PressureLayer, mupcl: nwsspc.sharp.calc.parcel.Parcel, leftMover: bool = False) -> nwsspc.sharp.calc.winds.WindComponents:
    """
    Estimates supercell storm motion using the effective-inflow mean wind of
    Bunkers et al. (2014).

    The mean wind is pressure weighted, with the height of the effective
    inflow base (see effective_inflow_layer) as its base and 65% of the
    height of the most-unstable parcel's Equilibrium Level as its top, both
    in meters AGL. The storm moves 7.5 m/s from that mean wind, perpendicular
    to the 0-6 km AGL shear between the 0-0.5 km and 5.5-6 km mean winds, as
    in the Bunkers 2000 method. Bunkers et al. found this does as well as the
    Bunkers 2000 method overall, and better for elevated supercells.

    This falls back to the Bunkers 2000 method with 0-6 km AGL layers when:

    - the effective inflow layer or the parcel's EL pressure is MISSING;
    - in QC builds, the effective inflow layer or the EL can't be converted
      to height (see nwsspc.sharp.calc.layer.pressure_layer_to_height), for
      example because it lies outside the profile; or
    - the mean wind layer would be less than 3 km deep, that is, 65% of the
      EL height AGL is less than 3 km above the inflow base.

    The inflow base and the EL are converted to meters AGL internally, so
    height may be AGL or MSL, and the result doesn't depend on the station
    elevation.

    The input parameters of eff_infl_lyr and mupcl (effective inflow layer 
    pressure bounds and the most unstable parcel, respectively) are required
    to be precomputed and passed to this routine. These are expensive 
    operations that are presumed to be computed at some other point 
    in the analysis pipeline. 

    References
    ----------
    Bunkers, M. J., D. A. Barber, R. L. Thompson, R. Edwards, and J. Garner,
    2014: Choosing a universal mean wind for supercell motion prediction.
    J. Operational Meteor., 2 (11), 115-129,
    https://doi.org/10.15191/nwajom.2014.0211

    Parameters 
    ----------
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height values (meters)
    u_wind : numpy.ndarray[dtype=float32] 
        1D NumPy array of U wind component values (m/s)
    v_wind : numpy.ndarray[dtype=float32] 
        1D NumPy array of V wind component values (m/s)
    eff_infl_lyr : nwsspc.sharp.calc.layer.PressureLayer 
        Effective Inflow Layer PressureLayer
    mupcl : nwsspc.sharp.calc.parcel.Parcel 
        Most Unstable Parcel 
    leftMover : bool 
        Whether or not to compute left moving supercell motion (default: False)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindComponents
        U, V wind components of storm motion (m/s)
    """

def mcs_motion_corfidi(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[nwsspc.sharp.calc.winds.WindComponents, nwsspc.sharp.calc.winds.WindComponents]:
    """
    Compute the Corfidi upshear and downshear MCS motion vectors.

    Estimates the mesoscale convective system (MCS) motion vectors for upshear 
    and downshear propagating convective systems as in Corfidi et al. 2003.
    The method is based on observations that MCS motion is a function of 
    1) the advection of existing cells by the mean wind and 
    2) the propagation of new convection relative to existing storms.

    In QC builds, returns MISSING components for both vectors if the 0-1.5 km
    AGL layer can't be converted to pressure (see
    nwsspc.sharp.calc.layer.height_layer_to_pressure), for example because the
    profile ends below 1.5 km AGL.

    References
    ----------
    Corfidi et al. 2003: https://www.spc.noaa.gov/publications/corfidi/mcs2003.pdf

    Parameters 
    ----------
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height values (meters)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of u-wind components (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of v-wind components (m/s)

    Returns 
    -------
    tuple[nwsspc.sharp.calc.winds.WindComponents, nwsspc.sharp.calc.winds.WindComponents]
        (upshear, downshear)
    """

def effective_bulk_wind_difference(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], effective_inflow_layer: nwsspc.sharp.calc.layer.PressureLayer, equilibrium_level_pressure: float) -> nwsspc.sharp.calc.winds.WindComponents:
    """
    Compute the Effective Bulk Wind Difference 

    The effective bulk wind difference is the wind shear over the layer
    from the base of the effective inflow layer up half the distance from
    that base to the equilibrium level, normally that of the most unstable
    parcel (Thompson et al. 2007). This is analogous to the usage 
    of 0-6 km wind shear, but allows more flexibility for elevated 
    convection. The effective inflow layer and the equilibrium level are
    converted to meters AGL internally, so height may be AGL or MSL, and the
    result doesn't depend on the station elevation.

    Returns MISSING if the effective inflow layer or 
    equilibrium level pressure are MISSING. In QC builds it also returns
    MISSING if either can't be converted to height (see
    nwsspc.sharp.calc.layer.pressure_layer_to_height), for example because it
    lies outside the profile.

    Parameters 
    ----------
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height values (meters)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of u-wind components (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of v-wind components (m/s)
    effective_inflow_layer : nwsspc.sharp.calc.layer.PressureLayer 
        The PressureLayer that defines the Effective Inflow layer 
    equilibrium_level_pressure : float 
        The pressure of the equilibrium level of the most unstable parcel (Pa)

    Returns
    -------
    nwsspc.sharp.calc.winds.WindComponents 
        The (u, v) wind components of the effective bulk wind difference (m/s)
    """

def energy_helicity_index(cape: float, helicity: float) -> float:
    """
    Computes the Energy Helicity Index.

    EHI is a composite parameter based on the premise that 
    storm rotation shoudl be maximized when CAPE is large 
    and SRH is large. Typically, the layers used for helicity 
    are either 0-1 km AGL or 0-3 km AGL.

    References 
    ----------
    https://doi.org/10.1175/1520-0434(2003)18%3C530:RSATFP%3E2.0.CO;2

    Parameters 
    ----------
    CAPE : float 
        Convective Available Potential Energy (J/kg)
    helicity : float 
        Storm Relative Helicity (m^2 / s^2 a.k.a J/kg)

    Returns
    -------
    float
        Energy Helicity Index (umitless)
    """

def significant_tornado_parameter(pcl: nwsspc.sharp.calc.parcel.Parcel, lcl_hght_agl: float, storm_relative_helicity: float, bulk_wind_difference: float) -> float:
    """
    Computes the Significant Tornado Parameter.


    The Significant Tornado Parameter is used to diagnose environments
    where tornadoes are favored. STP traditionally comes in two flavors:
    fixed-layer, and effective-layer. Fixed-layer STP expects surface-based
    CAPE, the surface-based LCL, 0-1 km storm-relative helicity, the
    0-6 km bulk wind difference, and the surface-based CINH. For the
    effective inflow layer based STP, use 100mb mixed-layer CAPE,
    100mb mixed-layer LCL height AGL, effective-layer srh, the
    effective layer bulk wind difference, and the 100mb mixed-layer
    CINH. NOTE: The effective bulk wind difference is the shear between
    the bottom of the effective inflow layer and 50% of the height of the
    equilibrium level of the most unstable parcel.

    References
    ----------
    Thompson et al 2012: https://www.spc.noaa.gov/publications/thompson/waf-env.pdf

    Parameters 
    ----------
    pcl : nwsspc.sharp.calc.parcel.Parcel 
        For effective-layer STP, a mixed-layer parcel, and for fixed-layer STP, a surface-based parcel 
    lcl_hght_agl : float 
        The parcel LCL height in meters
    storm_relative_helicity : float 
        For effective-layer STP, effective SRH, and for fixed-layer, 0-1 km SRH (m^2 / s^2)
    bulk_wind_difference : float 
        For effective-layer STP, effective BWD, and for fixed-layer STP, 0-6 km BWD (m/s)

    Returns
    -------
    float 
        The Significant Tornado Parameter
    """

def supercell_composite_parameter(mu_cape: float, eff_srh: float, eff_shear: float) -> float:
    """
    Computes the Supercell Composite Parameter. 

    The supercell composite parameter is used to diagnose environments
    where supercells are favored. Requires computing most unstable
    CAPE, effective layer storm relative helicity, and effective
    bulk shear. Effective bulk shear is the vector difference between
    the winds at the bottom of the effective inflow layer, and 50% of
    the equilibrium level height. It is similar to the 0-6 km shear
    vector, but allows for elevated supercell thunderstorms.

    The left-moving supercell composite parameter can be computed by
    providing effective SRH calculated using the bunkers left-moving
    storm motion, and will return negative values.

    References
    ----------
    Thompson et al 2003: https://www.spc.noaa.gov/publications/thompson/ruc_waf.pdf

    Thompson et al 2007: https://www.spc.noaa.gov/publications/thompson/effective.pdf

    Thompson et al 2012: https://www.spc.noaa.gov/publications/thompson/waf-env.pdf



    Parameters 
    ----------
    mu_cape : float 
        The CAPE of the Most Unstable Parcel (J/kg)
    eff_srh : float 
        Effective inflow layer Storm Relative Helicity (m^2/s^2) 
    eff_shear : float 
        Effective layer shear (m/s)

    Returns
    -------
    float
        Supercell Composite Parameter (unitless)
    """

def significant_hail_parameter(mu_pcl: nwsspc.sharp.calc.parcel.Parcel, lapse_rate_700_500mb: float, tmpk_500mb: float, freezing_level_agl: float, shear_0_6km: float) -> float:
    """
    Compute the significant hail parameter, given a precomputed most-unstable parcel, 
    the 700-500 mb lapse rate, the 500mb temperature, the height (AGL) of the 
    freezing level, and the 0-6 km shear magnitude.

    The Sig. Hail Parameter (SHIP) was developed using a large database of 
    surface-modified, observed severe hail proximity soundings. It is based on 
    parameters, and is meant to delineate between SIG (>=2" diameter) and NON-SIG
    (<2" diameter) hail environments.

    SHIP = [(MUCAPE j/kg) * (Mixing Ratio of MU PARCEL g/kg) *  
            (700-500mb LAPSE RATE c/km) * (-500mb TEMP C) *
            (0-6km Shear m/s) ] / 42,000,000

    0-6 km shear is confined to a range of 7-27 m s-1, mixing ratio is confined to 
    a range of 11-13.6 g kg-1, and the 500 mb temperature is set to -5.5 C for 
    any warmer values.

    Once the initial version of SHIP is calculated, the values are modified in 
    the following scenarios:

    1) If MUCAPE < 1300 J kg-1, SHIP = SHIP * (MUCAPE/1300); 2) if 700-500 mb 
    lapse rate < 5.8 C km-1, SHIP = SHIP * (lr75/5.8); 3) if freezing 
    level < 2400 m AGL, SHIP = SHIP * (fzl/2400)

    It is important to note that SHIP is NOT a forecast hail size.

    Since SHIP is based on the RAP depiction of MUCAPE - unrepresentative MUCAPE 
    "bullseyes" may cause a similar increase in SHIP values. This typically occurs 
    when bad surface observations get into the RAP model.

    Developed in the same vein as the STP and SCP parameters, values of SHIP 
    greater than 1.00 indicate a favorable environment for SIG hail. Values greater 
    than 4 are considered very high. In practice, maximum contour values of 1.5-2.0 
    or higher will typically be present when SIG hail is going to be reported. 

    Parameters 
    ----------
    mu_pcl : nwsspc.sharp.calc.parcel.Parcel 
        A precomputed Most Unstable parcel
    lapse_rate_700_500mb : float
        The 700-500 mb lapse rate (K/km)
    tmpk_500mb : float 
        The 500mb temperature (K)
    freezing_level_agl : float 
        The height of the freezing level (AGL, meters)
    shear_0_6km : float
        The 0-6 km shear vector magnitude (m/s)

    Returns
    -------
    float
        The significant hail parameter
    """

def derecho_composite_parameter(dcape: float, mucape: float, shear_0_6km: float, mean_wind_0_6km: float) -> float:
    """
    The Derecho Composite Parameter (DCP) is based on a dataset of 113 derecho 
    events compiled by the Evans and Boswell (2001) study. It is intended to 
    hightlight environments favorable for cold pool/outflow driven convective 
    events. The physical mechanisms behind this parameter focus on cold pool 
    production (DCAPE), ability to sustain strong convection (MUCAPE), 
    convective organization (0 - 6 km shear), and sufficient deep-layer flow 
    within the environment.

    References
    ----------
    Evans and Doswell 2001: https://doi.org/10.1175/1520-0434(2001)016%3C0329:EODEUP%3E2.0.CO;2

    Parameters
    ----------
    dcape : float 
        Downdraft Convective Available Potential Energy (J/kg)
    mucape : float 
        Most Unstable Parcel Convective Available Potential Energy (J/kg)
    shear_0_6km : float 
        Shear magnitude in the 0 - 6 km layer AGL (m/s)
    mean_wind_0_6km : float 
        Mean wind magnitude in the 0 - 6 km layer AGL (m/s)

    Returns
    -------
    float 
        The Derecho Composite Parameter
    """

def large_hail_parameter(mu_pcl: nwsspc.sharp.calc.parcel.Parcel, lapse_rate_700_500mb: float, hail_growth_zone: nwsspc.sharp.calc.layer.PressureLayer, storm_motion: nwsspc.sharp.calc.winds.WindComponents, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Computes the Large Hail Parameter (LHP). The LHP is a multi-ingredient, 
    composite index that includes thermodynamics and kinematics to attempt 
    to detect environments that support very large hail. LHP has shown skill 
    when differentiationg environments that support hail >= 3.5 in from those 
    with < 2.0 in.

    In QC builds, returns MISSING if the hail growth zone, the 1500 m below
    the equilibrium level, or the 0-1 km or 3-6 km AGL layer can't be
    converted between height and pressure (see
    nwsspc.sharp.calc.layer.pressure_layer_to_height and
    nwsspc.sharp.calc.layer.height_layer_to_pressure), for example because it
    extends past the profile or the hail growth zone is MISSING. It also
    returns MISSING if the equilibrium level can't be interpolated to height, for
    example because it lies above the top of the profile.

    References
    ----------
    Johnson and Sugden 2014: https://ejssm.org/archives/wp-content/uploads/2021/09/vol9-5.pdf

    Parameters 
    ----------
    mu_pcl : nwsspc.sharp.calc.parcel.Parcel 
        A previously computed most unstable parcel with CAPE and an equilibrium level
    lapse_rate_700_500mb : float 
        700 - 500 hPa Lapse Rate (K)
    hail_growth_zone : nwsspc.sharp.calc.layer.PressureLayer
        The layer of the atmosphere encompasing the Hail Growth Zone (Pa)
    storm_motion : nwsspc.sharp.calc.winds.WindComponents 
        The storm motion vector to be used (m/s)
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height values (m)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of u_wind values (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of v_wind values (m/s)

    Returns 
    -------
    float 
       Large Hail Parameter
    """

@overload
def convective_temperature(lifter: nwsspc.sharp.calc.parcel.lifter_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], cinh_thresh: float = -1.0) -> float:
    """
    Computes the convective temperature by iteratively lifting parcels from 
    the surface using a lowest 100 hPa mean mixing ratio and increasing 
    surface temperatures to find the first parcel that reaches the CINH 
    threshold. The first guess is the current surface temperature.

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_wobus
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    virtemp : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    mixratio : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of water vapor mixing ratio values (g/g)
    cinh_thresh : float, default = -1.0 
        The CINH threshold used to compute the convective temperature (J/kg)

    Returns
    -------
    float
        The convective temperature (K)
    """

@overload
def convective_temperature(lifter: nwsspc.sharp.calc.parcel.lifter_lut_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], cinh_thresh: float = -1.0) -> float:
    """
    Computes the convective temperature by iteratively lifting parcels from 
    the surface using a lowest 100 hPa mean mixing ratio and increasing 
    surface temperatures to find the first parcel that reaches the CINH 
    threshold. The first guess is the current surface temperature.

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_lut_wobus
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    virtemp : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    mixratio : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of water vapor mixing ratio values (g/g)
    cinh_thresh : float, default = -1.0 
        The CINH threshold used to compute the convective temperature (J/kg)

    Returns
    -------
    float
        The convective temperature (K)
    """

@overload
def convective_temperature(lifter: nwsspc.sharp.calc.parcel.lifter_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], cinh_thresh: float = -1.0) -> float:
    """
    Computes the convective temperature by iteratively lifting parcels from 
    the surface using a lowest 100 hPa mean mixing ratio and increasing 
    surface temperatures to find the first parcel that reaches the CINH 
    threshold. The first guess is the current surface temperature.

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_cm1
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    virtemp : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    mixratio : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of water vapor mixing ratio values (g/g)
    cinh_thresh : float, default = -1.0 
        The CINH threshold used to compute the convective temperature (J/kg)

    Returns
    -------
    float
        The convective temperature (K)
    """

@overload
def convective_temperature(lifter: nwsspc.sharp.calc.parcel.lifter_lut_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], cinh_thresh: float = -1.0) -> float:
    """
    Computes the convective temperature by iteratively lifting parcels from 
    the surface using a lowest 100 hPa mean mixing ratio and increasing 
    surface temperatures to find the first parcel that reaches the CINH 
    threshold. The first guess is the current surface temperature.

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_lut_cm1
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    virtemp : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    mixratio : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of water vapor mixing ratio values (g/g)
    cinh_thresh : float, default = -1.0 
        The CINH threshold used to compute the convective temperature (J/kg)

    Returns
    -------
    float
        The convective temperature (K)
    """

def precipitable_water(layer: nwsspc.sharp.calc.layer.PressureLayer, pres: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
    """
    Given a PressureLayer to integrate over, compute the precipitable water 
    from the given pressure and mixing ratio arrays.

    Parameters 
    ----------
    layer : nwsspc.sharp.calc.layer.PressureLayer 
        a PressureLayer over which to integrate (Pa)
    pres : numpy.ndarray[dtype=float32] 
        1D NumPy array of presssure values (Pa)
    mixr : numpy.ndarray[dtype=float32] 
        1D NumPy array of water vapor mixing ratio values (unitless)

    Returns
    -------
    float
        Precipitable water content (mm)
    """

def hail_growth_layer(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> nwsspc.sharp.calc.layer.PressureLayer:
    """
    Performs a top-down search for the hail growth zone, 
    returning the PressureLayer bounds. This zone is defined 
    as the pressure layer that contains the portion of the 
    profile between -30C and -10C.

    Parameters 
    ----------
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of temperature values (K)

    Returns 
    -------
    nwsspc.sharp.calc.layer.PressureLayer
        The PressureLayer containing the hail growth zone
    """

def dendritic_layer(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> nwsspc.sharp.calc.layer.PressureLayer:
    """
    Performs a top-down search for the dendritic growth zone, 
    returning the PressureLayer bounds. This zone is defined 
    as the pressure layer that contains the portion of the 
    profile between -17C and -12C.

    Parameters 
    ----------
    pressure : numpy.ndarray[dtype=float32] 
        1D NumPy array of pressure values (Pa)
    temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of temperature values (K)

    Returns 
    -------
    nwsspc.sharp.calc.layer.PressureLayer
        The PressureLayer containing the dendritic growth zone
    """

def snow_squall_parameter(wetbulb_2m: float, mean_relh_0_2km: float, delta_thetae_0_2km: float, mean_wind_0_2km: float) -> float:
    """
    The Snow Squall Parameter is a non-dimensional parameter that combines 
    several ingredients believed to be beneficial for identifying snow squall 
    environments by identifying the overlap of low-level potential instability, 
    sufficient moisture, and strong low-level winds.


    References
    ----------
    Banacos et al. 2014: https://www.weather.gov/media/btv/research/Snow%20Squalls%20Forecasting%20and%20Hazard%20Mitigation.pdf

    Parameters 
    ----------
    wetbulb_2m : float 
        The surface wetbulb temperature, used to mask the parameter (K)
    mean_relh_0_2km : float 
        The mean relative humidity between the surface and 2 km AGL (fraction)
    delta_thetae_0_2km : float 
        The difference in equivalent potential temperature between 2 km AGL and the surface (K)
    mean_wind_0_2km : float 
        The mean wind speed between the surface and 2 km AGL (m/s)

    Returns 
    -------
    float 
        The Snow Squall Parameter
    """

BOURGOUIN_PRESSURE_MIN: float = 25000.0

class BourgouinEnergy:
    """
    Wet-bulb melting and refreezing energies of the modified Bourgouin method.

    Energies are areas between the wet-bulb temperature profile and 0 C
    (Birk et al. 2021, Eq. 1 with the wet-bulb temperature), reported as
    positive values in J/kg. The near-surface cold layer is the lowest layer
    colder than 0 C that has a warmer-than-0 C layer above it.

    Every field defaults to MISSING.

    References
    ----------
    Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, melting_energy_total: float, melting_energy_aloft: float, refreezing_energy: float) -> None: ...

    @property
    def melting_energy_total(self) -> float:
        """
        Total wet-bulb melting energy of the column (J/kg). Used for snow and for freezing rain or rain.
        """

    @melting_energy_total.setter
    def melting_energy_total(self, arg: float, /) -> None: ...

    @property
    def melting_energy_aloft(self) -> float:
        """
        Wet-bulb melting energy above the near-surface cold layer (J/kg). Used for ice pellets. 0 when there is no near-surface cold layer.
        """

    @melting_energy_aloft.setter
    def melting_energy_aloft(self, arg: float, /) -> None: ...

    @property
    def refreezing_energy(self) -> float:
        """
        Wet-bulb refreezing energy of the near-surface cold layer (J/kg, positive). 0 when there is no near-surface cold layer.
        """

    @refreezing_energy.setter
    def refreezing_energy(self, arg: float, /) -> None: ...

class PrecipTypeProbabilities:
    """
    Probabilities of the four basic precipitation types.

    Each field is a fraction in [0, 1]. The probabilities are independent and
    do not sum to 1. Two types with high probabilities describe a mix of the
    two, not a contradiction.

    Every field defaults to MISSING.
    """

    def __init__(self) -> None: ...

    @property
    def rain(self) -> float:
        """Probability of rain (fraction)"""

    @rain.setter
    def rain(self, arg: float, /) -> None: ...

    @property
    def snow(self) -> float:
        """Probability of snow (fraction)"""

    @snow.setter
    def snow(self, arg: float, /) -> None: ...

    @property
    def freezing_rain(self) -> float:
        """Probability of freezing rain (fraction)"""

    @freezing_rain.setter
    def freezing_rain(self, arg: float, /) -> None: ...

    @property
    def ice_pellets(self) -> float:
        """Probability of ice pellets (fraction)"""

    @ice_pellets.setter
    def ice_pellets(self, arg: float, /) -> None: ...

def bourgouin_energy(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], wetbulb: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], min_energy: float = 0.0, pressure_min: float = 25000.0) -> BourgouinEnergy:
    """
    Computes the wet-bulb melting and refreezing energies of the modified
    Bourgouin method from a sounding.

    The energies are the areas of Birk et al. (2021, Eq. 1) with the wet-bulb
    temperature Tw: g (Tw - T0) / T0 integrated over height, with
    T0 = 273.15 K. Layers are bounded by linearly interpolated 0 C crossings,
    and the areas are trapezoids split exactly at those crossings. The
    function reports every energy as a positive value in J/kg.

    * melting_energy_total (ME_total): all energy above 0 C in the column.
    * refreezing_energy (RE): the energy of the near-surface cold layer, the
      lowest layer below 0 C that has a layer above 0 C over it. 0 if there
      is no such layer.
    * melting_energy_aloft (ME_aloft): all energy above 0 C over the top of
      the near-surface cold layer. 0 if there is no such layer.

    modified_bourgouin uses ME_total for snow and for freezing rain or rain,
    and ME_aloft for ice pellets. With the default min_energy = 0, ME_total
    equals ME_aloft whenever the lowest layer is cold. When the only warm layer
    is at the surface (Fig. 1b), ME_total holds all of the melting energy and
    ME_aloft is 0.

    With several warm layers, RE comes from the near-surface cold layer only,
    and ME_aloft adds up every warm layer above it. Other cold layers never
    enter the equations. Example, from the surface up: cold 100, warm 30,
    cold 80, and warm 20 J/kg give RE = 100, ME_aloft = 50, and ME_total = 50.

    A warm layer at the surface, below the near-surface cold layer (Fig. 1d),
    adds to ME_total only. As in the paper, it does not suppress ice pellets.
    Example, from the surface up: warm 150, cold 51, and warm 10 J/kg give
    ME_total = 160, ME_aloft = 10, and RE = 51. With ProbIce = 1,
    modified_bourgouin then gives rain 100 % and ice pellets about 20 %
    (19.6 %).

    min_energy sets the energy a layer needs to stand on its own. Weaker layers
    merge into their neighbors. Two warm or two cold layers absorb the weak
    layers between them. Weak layers between a warm and a cold layer all go to
    whichever kind makes up more of their depth, and a tie goes to the lower
    layer. Weak layers at the bottom or top of the profile join the layer next
    to them. Only a layer's own energy counts, never what it absorbs. The
    paper sets no minimum melting energy for the modified method, so the
    default of 0 is the paper as written and merges nothing. A positive value
    is an opt-in for noisy, high-resolution data, and a deviation from the
    paper. It changes more than the onset of ice pellets:

    * A warm or cold layer with less than min_energy no longer splits the
      layers around it, so ME_aloft is either 0 or at least min_energy.
    * A weak warm layer between two cold layers merges them into one
      near-surface cold layer, and RE includes both. Example: cold 100, warm
      1.99, cold 80, and warm 20 J/kg give RE = 180 with min_energy = 2, but
      RE = 100 with min_energy = 0. RE jumps as the weak layer crosses the
      threshold.
    * ME_total still counts every warm layer, including the merged ones, but
      ME_aloft leaves out a warm layer merged into the near-surface cold layer.
      So ME_total and ME_aloft can differ over a cold surface. Example: cold
      100, warm 1.99, cold 80, and warm 2.01 J/kg give ME_total = 4,
      ME_aloft = 2.01, and RE = 180 with min_energy = 2. The weak-melting
      taper of modified_bourgouin acts on ME_total.

    Consider min_energy for 1 Hz soundings and other high-resolution profiles
    whose wet-bulb temperature stays near 0 C over some depth, such as an
    isothermal melting layer or a surface layer close to 0 C. There, noise
    makes the profile cross 0 C many times. A warm sliver over a cold surface
    layer makes that layer a near-surface cold layer, ME_aloft and RE become
    positive, and the ice pellet probability jumps from 0 to about
    2.3 RE + 3 %. Start with 2 J/kg, the melting-layer minimum of the original
    Bourgouin method. On the three 1 Hz soundings in the SHARPlib test data,
    noise alone on a layer at exactly 0 C makes wet-bulb layers of which 99 %
    hold less than 1.2 J/kg, and the largest held 2.8 J/kg. With that noise
    added to a saturated 800 m layer at -0.15 C, a wet-bulb temperature of
    -0.28 C, over a cold surface, spurious ice pellets appeared in 25 % of
    trials with min_energy = 0 and in none with 0.5 J/kg or more. With the
    noise doubled, it took 2 J/kg. The cost is that warm layers weaker than
    min_energy no longer count, and noise can split a slightly stronger layer
    into pieces that each fall below it. A 3.1 J/kg melting layer was lost in
    7 % of trials at 2 J/kg and in 31 % at 3 J/kg.

    min_energy never changes ME_total, which counts every warm layer. On 1 Hz
    data, noise moved ME_total by up to about 1 J/kg in 95 % of trials. With
    ME_total below 5 J/kg, that alone moved rain or freezing rain by 0.1 or
    more in 9 to 32 % of trials, whatever the options.

    Only levels at pressures at or above pressure_min (by default 25000 Pa,
    250 hPa) are used, which keeps stratospheric temperatures out of the
    melting energy. The column ends at the highest such level, with nothing
    interpolated to pressure_min itself. A pressure_min of 0 uses every level.
    ME_total depends on how far up the data reach, up to that limit.

    The function never reads a wet-bulb temperature at a pressure below
    pressure_min. A caller can compute the wet-bulb temperature only up to
    pressure_min and fill the rest of the array with MISSING.

    The function skips levels whose wet-bulb temperature is MISSING or NaN,
    and joins the valid levels on either side with a straight line. With
    fewer than 2 valid levels at pressures at or above pressure_min, including
    empty arrays, every energy is MISSING, never 0, since zero energy would
    read as certain snow.

    Height must be strictly increasing, and pressure must be valid and
    strictly decreasing. This is not checked.

    References
    ----------
    Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1

    Bourgouin 2000: https://doi.org/10.1175/1520-0434(2000)015%3C0583:AMTDPT%3E2.0.CO;2

    Parameters
    ----------
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height values (m)
    wetbulb : numpy.ndarray[dtype=float32]
        1D NumPy array of wet-bulb temperature values (K)
    min_energy : float, default = 0.0
        The energy a layer needs to stand on its own (J/kg; 0 disables merging)
    pressure_min : float, default = 25000.0
        Levels at lower pressures are ignored (Pa)

    Returns
    -------
    nwsspc.sharp.calc.params.BourgouinEnergy
        The total melting energy, the melting energy aloft, and the refreezing
        energy (J/kg)
    """

def precipitation_generation_layer(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], min_depth: float = 0.0) -> nwsspc.sharp.calc.layer.HeightLayer:
    """
    Finds the precipitation generation layer of the modified Bourgouin
    method.

    Splits the profile into moist layers, where relative humidity is above
    75 %, and dry layers, where it is below 75 %. Returns the highest moist
    layer deeper than 1000 m that lies below the first dry layer deeper than
    1500 m (Birk et al. 2021, section 3e). The paper assumes that
    precipitation falling from above such a dry layer sublimates, so the dry
    layer eliminates every layer above it, even when it starts at the
    surface. Both depths are strict. A 1000 m moist layer is not a
    generation layer, and a 1500 m dry layer eliminates nothing. Layer
    boundaries are the linearly interpolated 75 % crossings. The function
    walks up the profile once and stops at the first eliminating dry layer.

    Relative humidity is over ice where the air temperature is below 0 C and
    over liquid water otherwise, from relative_humidity_ice and
    relative_humidity in nwsspc.sharp.calc.thermo. The paper uses relative
    humidity over ice at every temperature, so this deviates from it above
    0 C. The two agree at 0 C. For example, T = 283.15 K with Td = 280 K is
    0.808 over liquid, which is moist, but 0.733 over ice, which is dry.

    Levels at exactly 75 % continue the current layer and add their depth to
    it, so a layer ends only where the relative humidity crosses 75 %. For
    example, relative humidities of 0.8, 0.75, 0.75, and 0.8 at 0, 100, 2000,
    and 2100 m form one 2100 m moist layer, which is a generation layer. The
    paper defines moist as above 75 % and dry as below it, and doesn't say
    how levels at exactly 75 % count. Read strictly, they belong to neither
    kind of layer and end both. That reading turns the example into two
    100 m moist layers and no generation layer. The two readings differ only
    where the relative humidity is exactly 0.75.

    min_depth is an opt-in filter for noisy, high-resolution profiles. Moist
    and dry runs shallower than min_depth merge into the layers around them.
    Two layers of the same kind absorb the shallow runs between them.
    Shallow runs between a moist and a dry layer all go to whichever kind
    makes up more of their depth, and a tie goes to the lower layer. Shallow
    runs at the bottom or top of the profile join the layer next to them. So
    a dry layer's depth includes the shallow moist runs it absorbs. The
    default of 0 merges nothing, which is the paper as written. With
    min_depth above 0, a dry layer's depth is final once the moist run above
    it is min_depth deep, so the walk stops there.

    Consider min_depth for 1 Hz soundings and other high-resolution profiles
    where a moist layer is close to 1000 m deep, a dry layer is close to
    1500 m deep, or the relative humidity stays near 75 %. Noise there splits
    and joins layers, which can move the result between MISSING, a low cloud,
    and a higher, colder one. Start with 100 m. On the three 1 Hz soundings in
    the SHARPlib test data, noise alone on a layer at exactly 75 % makes
    humidity layers of which 90 % are under 110 m deep and 99 % under 200 m.
    With that noise added to profiles whose layers were 150 to 250 m from
    those depths, the generation layer changed in 19 to 29 % of trials with
    min_depth = 0, in 3 to 12 % with 50 m, and in at most 0.2 % with 100 m.
    Larger values merge real layers. One test sounding has a 324 m dry layer
    at the surface under a 684 m moist layer. With 300 m, noise often thinned
    the dry layer below 300 m, and the moist layer absorbed it and became a
    generation layer in 17 % of trials, against under 1 % with 200 m.

    A result of MISSING does not say why. The profile may have no moist layer
    deeper than 1000 m, a deep dry layer under the cloud may eliminate it
    (virga), or the moisture data may be missing. Callers can't tell these
    cases apart from the result.

    The walk skips levels with a MISSING or NaN temperature or dewpoint, and
    joins their valid neighbors with a straight line.

    Heights must be strictly increasing. This is not checked.

    References
    ----------
    Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1

    Parameters
    ----------
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height values (meters)
    temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of temperature values (K)
    dewpoint : numpy.ndarray[dtype=float32]
        1D NumPy array of dewpoint temperature values (K)
    min_depth : float, default = 0.0
        Moist and dry runs shallower than this merge into the layers around
        them (meters; 0 disables merging)

    Returns
    -------
    nwsspc.sharp.calc.layer.HeightLayer
        The precipitation generation layer (meters, AGL or MSL like height).
        Its bottom and top are MISSING if there is none.
    """

def probability_of_ice(temperature: float) -> float:
    """
    Computes the probability of heterogeneous ice nucleation (ProbIce) from
    the minimum air temperature T in the precipitation generation layer. It
    uses the piecewise relation of Baumgardt et al. (2017), as reproduced by
    Birk et al. (2021, Eq. 2 and Fig. 2):

    * 1 for T at or below -15 C
    * 0 for T at or above -7 C
    * otherwise (-0.065 T^4 - 3.1544 T^3 - 56.414 T^2 - 449.6 T - 1308) / 100,
      with T in C

    The polynomial does not meet the constant pieces. It gives 98.3 % at
    -15 C and 0.81 % at -7 C, so ProbIce jumps by about 0.017 at -15 C and
    about 0.008 at -7 C, as in the paper.

    MISSING or NaN input returns MISSING.

    Warnings
    --------
    The input is in Kelvin, and Celsius input gives wrong results. A value at
    or below 0 returns MISSING, and a positive value reads as a Kelvin
    temperature far below -15 C and returns 1.

    References
    ----------
    Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1

    Baumgardt et al. 2017: https://ams.confex.com/ams/97Annual/webprogram/Paper313165.html

    Parameters
    ----------
    temperature : float
        Minimum air temperature in the precipitation generation layer (K)

    Returns
    -------
    float
        Probability of ice (fraction)
    """

@overload
def modified_bourgouin(energy: BourgouinEnergy, prob_ice: float, surface_wetbulb: float) -> PrecipTypeProbabilities:
    """
    Computes the probabilities of rain, snow, freezing rain, and ice pellets
    from wet-bulb energies and the probability of ice, with the modified
    Bourgouin method. It follows the steps in the Appendix of Birk et al.
    (2021). The function works in percent and returns fractions. ProbIce is
    prob_ice in percent. Every clamp is to [0, 100] %.

    * Snow (Eq. 9): 1540 exp(-0.29 ME_total), clamped, then multiplied by
      ProbIce / 100 and clamped.
    * Ice pellets (Eq. 8): when both ME_aloft and RE are above 0,
      2.3 RE - 42 ln(ME_aloft + 1) + 3, clamped; otherwise 0. Then multiplied
      by ProbIce / 100 and clamped.
    * Freezing rain or rain (Eq. 7): -2.1 RE + 0.2 ME_total + 458, clamped
      first. When ME_total is below 5 J/kg, that clamped value is then
      multiplied by 0.2 ME_total. The result is
      (100 - ProbIce) + (ProbIce / 100) times that value, clamped.

    ME_total, ME_aloft, and RE are the melting_energy_total,
    melting_energy_aloft, and refreezing_energy fields of energy. Freezing
    rain and rain use the total melting energy, and ice pellets use the
    melting energy above the near-surface cold layer.

    The four probabilities are independent and do not sum to 1. Freezing rain
    and rain share one value. The function reports it as rain when
    surface_wetbulb is above 0 C (273.15 K) and as freezing rain otherwise,
    including at exactly 0 C, and sets the other to 0. A warm surface does not
    suppress ice pellets.

    As in the paper, the ice pellet probability jumps from 0 to its Eq. 8
    value as ME_aloft rises from 0. Just above 0, Eq. 8 gives 2.3 RE + 3 %
    before the clamp and the ProbIce scaling.

    If any input is MISSING or NaN, all four probabilities are MISSING.

    References
    ----------
    Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1

    Parameters
    ----------
    energy : nwsspc.sharp.calc.params.BourgouinEnergy
        Wet-bulb melting and refreezing energies (J/kg)
    prob_ice : float
        Probability of ice in the precipitation generation layer (fraction;
        see probability_of_ice)
    surface_wetbulb : float
        Surface wet-bulb temperature (K)

    Returns
    -------
    nwsspc.sharp.calc.params.PrecipTypeProbabilities
        The rain, snow, freezing rain, and ice pellet probabilities (fractions)
    """

@overload
def modified_bourgouin(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], wetbulb: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], min_depth: float = 0.0, min_energy: float = 0.0, pressure_min: float = 25000.0) -> PrecipTypeProbabilities:
    """
    Computes the probabilities of rain, snow, freezing rain, and ice pellets of
    Birk et al. (2021) from pressure, height, temperature, dewpoint, and
    wet-bulb temperature profiles. It runs these steps, and calling them
    yourself with the same options gives the same result:

    1. precipitation_generation_layer with min_depth finds the precipitation
       generation layer.
    2. layer_min in nwsspc.sharp.calc.layer finds the minimum air temperature
       in that layer, and probability_of_ice turns it into ProbIce.
    3. bourgouin_energy with min_energy and pressure_min computes the wet-bulb
       melting and refreezing energies of the whole column, not just the
       generation layer.
    4. The surface wet-bulb temperature is that of the lowest level. Above
       0 C, the liquid probability is rain, otherwise freezing rain.
    5. The modified_bourgouin overload that takes energies combines them into
       the four probabilities, which are independent and do not sum to 1.

    With the defaults (min_depth = 0, min_energy = 0, and pressure_min =
    25000 Pa), the function follows the paper except in three ways:

    * The energies use only levels at pressures at or above 250 hPa. This has
      no effect on realistic tropospheric profiles.
    * Relative humidity is over ice below 0 C and over liquid water otherwise,
      where the paper uses relative humidity over ice at every temperature.
      For example, T = 283.15 K with Td = 280 K is moist over liquid (0.808)
      but dry over ice (0.733).
    * Levels at exactly 75 % relative humidity continue the current layer,
      where the paper puts them in neither the moist nor the dry class.

    Positive min_depth and min_energy are opt-ins for noisy, high-resolution
    data and further deviations from the paper. min_depth applies only to the
    generation layer, and min_energy and pressure_min only to the energies.
    For 1 Hz soundings, start with min_depth = 100 m and min_energy = 2 J/kg.
    precipitation_generation_layer and bourgouin_energy describe their
    effects, the measurements behind these values, and what the options cost.

    Every probability is MISSING when there is no generation layer. The result
    does not say why. The profile may have no moist layer deeper than 1000 m,
    a deep dry layer under the cloud may eliminate it (virga), or the moisture
    data may be missing. In particular, a cloud 1 km deep or less, such as a
    drizzle cloud, gives MISSING. A caller who wants a result for such a cloud
    can call bourgouin_energy and then the modified_bourgouin overload that
    takes energies, with prob_ice = 0, which treats the cloud as having no
    ice.

    ME_total, and with it snow and freezing rain or rain, depends on how far
    up the data reach, up to pressure_min. The wet-bulb temperature at
    pressures below pressure_min does not affect the result, so a caller can
    compute it only up to pressure_min and fill the rest of the array with
    MISSING.

    The steps skip MISSING and NaN levels as their own documentation
    describes, and the surface wet-bulb temperature is that of the lowest
    level whose wet-bulb temperature is not MISSING or NaN. With fewer than 2
    levels, including empty arrays, every probability is MISSING.

    The profiles must start at the surface. Height must be strictly
    increasing, and pressure must be valid and strictly decreasing. This is
    not checked.

    References
    ----------
    Birk et al. 2021: https://doi.org/10.1175/WAF-D-20-0118.1

    Parameters
    ----------
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height values (meters)
    temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of temperature values (K)
    dewpoint : numpy.ndarray[dtype=float32]
        1D NumPy array of dewpoint temperature values (K)
    wetbulb : numpy.ndarray[dtype=float32]
        1D NumPy array of wet-bulb temperature values (K)
    min_depth : float, default = 0.0
        Moist and dry runs shallower than this merge into the layers around
        them when finding the generation layer (meters; 0 disables merging)
    min_energy : float, default = 0.0
        The energy a layer needs to stand on its own in the energies (J/kg;
        0 disables merging)
    pressure_min : float, default = 25000.0
        Levels at lower pressures are ignored in the energies (Pa)

    Returns
    -------
    nwsspc.sharp.calc.params.PrecipTypeProbabilities
        The rain, snow, freezing rain, and ice pellet probabilities (fractions)
    """

@overload
def equilibrium_moisture_content(temperature: float, rel_humidity: float) -> float:
    """
    Compute the equilibrium moisture content for fuel 
    as in Simard (1968).

    Parameters
    ----------
    temperature : float 
        The air temperature (K)
    rel_humidity : float 
        Relative Humidity (fraction)

    Returns 
    -------
    float 
        Equilibrium Moisture Content (fraction)
    """

@overload
def equilibrium_moisture_content(temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], rel_humidity: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the equilibrium moisture content for fuel 
    as in Simard (1968).

    Parameters
    ----------
    temperature : numpy.ndarray[dtype=float32] 
        The air temperature (K)
    rel_humidity : numpy.ndarray[dtype=float32]
        Relative Humidity (fraction)

    Returns 
    -------
    numpy.ndarray[dtype=float32]
        Equilibrium Moisture Content (fraction)
    """

@overload
def fosberg_fire_index(temperature: float, rel_humidity: float, wind_speed: float) -> float:
    """
    Compute the Fosberg Fire-Weather Index (FWWI) as in Fosberg (1978).

    Parameters 
    ----------
    temperature : float 
        The air temperature (K)
    rel_humidity : float 
        The relative humidity (fraction)
    wind_speed : float 
        Wind speed (m/s)

    Returns 
    -------
    float 
        Fosberg Fire-Weather Index
    """

@overload
def fosberg_fire_index(temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], rel_humidity: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], wind_speed: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
    """
    Compute the Fosberg Fire-Weather Index (FWWI) as in Fosberg (1978).

    Parameters 
    ----------
    temperature : numpy.ndarray[dtype=float32]
        The air temperature (K)
    rel_humidity : numpy.ndarray[dtype=float32]
        The relative humidity (fraction)
    wind_speed : numpy.ndarray[dtype=float32]
        Wind speed (m/s)

    Returns 
    -------
    numpy.ndarray[dtype=float32]
        Fosberg Fire-Weather Index
    """

@overload
def pyrocumulonimbus_firepower_threshold(lifter: nwsspc.sharp.calc.parcel.lifter_wobus, mix_layer: nwsspc.sharp.calc.layer.PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], potential_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], pcl: nwsspc.sharp.calc.parcel.Parcel | None = None, phi: float = 6.67e-05, beta_incr: float = 0.005) -> float:
    """
    Computes the Pyrocumulonimbus Firepower Threshold (PFT), or the minimum 
    amount of firepower required to generate pyrocumulonimbus clouds for a 
    given atmospheric profile. Requires a PressureLayer to define a mixing 
    layer used to average values of potential temperature, mixing ratio, 
    and wind speed. The beta increment determines how to vary the plume 
    buoyancy factor, with smaller values resulting in more iteration steps. 
    Phi is the fire moisture to potential temperature increment ratio.

    Default values for beta_incr and phi are 0.005 and 6.67e-5, respectively.
    If a parcel is passed, the values will be set with the PFT fire parcel.

    In QC builds, returns MISSING if the mix-layer mean potential
    temperature, mixing ratio, or wind speed is MISSING, or if the potential
    temperature is MISSING at the LFC or at the level where the formula
    evaluates air density.

    References 
    ----------
    Tory et al. 2018: https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
    Tory et al. 2021: https://journals.ametsoc.org/view/journals/wefo/36/2/WAF-D-20-0027.1.xml

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_wobus
    mix_layer : nwsspc.sharp.calc.layer.PressureLayer
        A mixing layer to take average values of the environment from
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    mixratio : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of water vapor mixing ratio values (g/g)
    virtual_temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of U wind component values (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of V wind compnent values (m/s)
    potential_temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of potential temperature values (K)
    pcl : None or nwsspc.sharp.calc.parcel.Parcel, optional
        If a parcel is provided, returns the parcel for the PFT
    phi : float 
        The fire moisture to potential temperature increment ratio (kg/(kg*K))
    beta_incr : float 
        The fire plume buoyancy factor (unitless)

    Returns 
    -------
    float 
        The PyroCB Firepower Threshold (Watts)
    """

@overload
def pyrocumulonimbus_firepower_threshold(lifter: nwsspc.sharp.calc.parcel.lifter_lut_wobus, mix_layer: nwsspc.sharp.calc.layer.PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], potential_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], pcl: nwsspc.sharp.calc.parcel.Parcel | None = None, phi: float = 6.67e-05, beta_incr: float = 0.005) -> float:
    """
    Computes the Pyrocumulonimbus Firepower Threshold (PFT), or the minimum 
    amount of firepower required to generate pyrocumulonimbus clouds for a 
    given atmospheric profile. Requires a PressureLayer to define a mixing 
    layer used to average values of potential temperature, mixing ratio, 
    and wind speed. The beta increment determines how to vary the plume 
    buoyancy factor, with smaller values resulting in more iteration steps. 
    Phi is the fire moisture to potential temperature increment ratio.

    Default values for beta_incr and phi are 0.005 and 6.67e-5, respectively.
    If a parcel is passed, the values will be set with the PFT fire parcel.

    In QC builds, returns MISSING if the mix-layer mean potential
    temperature, mixing ratio, or wind speed is MISSING, or if the potential
    temperature is MISSING at the LFC or at the level where the formula
    evaluates air density.

    References 
    ----------
    Tory et al. 2018: https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
    Tory et al. 2021: https://journals.ametsoc.org/view/journals/wefo/36/2/WAF-D-20-0027.1.xml

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_lut_wobus
    mix_layer : nwsspc.sharp.calc.layer.PressureLayer
        A mixing layer to take average values of the environment from
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    mixratio : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of water vapor mixing ratio values (g/g)
    virtual_temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of U wind component values (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of V wind compnent values (m/s)
    potential_temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of potential temperature values (K)
    pcl : None or nwsspc.sharp.calc.parcel.Parcel, optional
        If a parcel is provided, returns the parcel for the PFT
    phi : float 
        The fire moisture to potential temperature increment ratio (kg/(kg*K))
    beta_incr : float 
        The fire plume buoyancy factor (unitless)

    Returns 
    -------
    float 
        The PyroCB Firepower Threshold (Watts)
    """

@overload
def pyrocumulonimbus_firepower_threshold(lifter: nwsspc.sharp.calc.parcel.lifter_cm1, mix_layer: nwsspc.sharp.calc.layer.PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], potential_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], pcl: nwsspc.sharp.calc.parcel.Parcel | None = None, phi: float = 6.67e-05, beta_incr: float = 0.005) -> float:
    """
    Computes the Pyrocumulonimbus Firepower Threshold (PFT), or the minimum 
    amount of firepower required to generate pyrocumulonimbus clouds for a 
    given atmospheric profile. Requires a PressureLayer to define a mixing 
    layer used to average values of potential temperature, mixing ratio, 
    and wind speed. The beta increment determines how to vary the plume 
    buoyancy factor, with smaller values resulting in more iteration steps. 
    Phi is the fire moisture to potential temperature increment ratio.

    Default values for beta_incr and phi are 0.005 and 6.67e-5, respectively.
    If a parcel is passed, the values will be set with the PFT fire parcel.

    In QC builds, returns MISSING if the mix-layer mean potential
    temperature, mixing ratio, or wind speed is MISSING, or if the potential
    temperature is MISSING at the LFC or at the level where the formula
    evaluates air density.

    References 
    ----------
    Tory et al. 2018: https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
    Tory et al. 2021: https://journals.ametsoc.org/view/journals/wefo/36/2/WAF-D-20-0027.1.xml

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_cm1
    mix_layer : nwsspc.sharp.calc.layer.PressureLayer
        A mixing layer to take average values of the environment from
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    mixratio : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of water vapor mixing ratio values (g/g)
    virtual_temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of U wind component values (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of V wind compnent values (m/s)
    potential_temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of potential temperature values (K)
    pcl : None or nwsspc.sharp.calc.parcel.Parcel, optional
        If a parcel is provided, returns the parcel for the PFT
    phi : float 
        The fire moisture to potential temperature increment ratio (kg/(kg*K))
    beta_incr : float 
        The fire plume buoyancy factor (unitless)

    Returns 
    -------
    float 
        The PyroCB Firepower Threshold (Watts)
    """

@overload
def pyrocumulonimbus_firepower_threshold(lifter: nwsspc.sharp.calc.parcel.lifter_lut_cm1, mix_layer: nwsspc.sharp.calc.layer.PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], u_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], v_wind: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], potential_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], pcl: nwsspc.sharp.calc.parcel.Parcel | None = None, phi: float = 6.67e-05, beta_incr: float = 0.005) -> float:
    """
    Computes the Pyrocumulonimbus Firepower Threshold (PFT), or the minimum 
    amount of firepower required to generate pyrocumulonimbus clouds for a 
    given atmospheric profile. Requires a PressureLayer to define a mixing 
    layer used to average values of potential temperature, mixing ratio, 
    and wind speed. The beta increment determines how to vary the plume 
    buoyancy factor, with smaller values resulting in more iteration steps. 
    Phi is the fire moisture to potential temperature increment ratio.

    Default values for beta_incr and phi are 0.005 and 6.67e-5, respectively.
    If a parcel is passed, the values will be set with the PFT fire parcel.

    In QC builds, returns MISSING if the mix-layer mean potential
    temperature, mixing ratio, or wind speed is MISSING, or if the potential
    temperature is MISSING at the LFC or at the level where the formula
    evaluates air density.

    References 
    ----------
    Tory et al. 2018: https://journals.ametsoc.org/view/journals/mwre/146/8/mwr-d-17-0377.1.xml
    Tory et al. 2021: https://journals.ametsoc.org/view/journals/wefo/36/2/WAF-D-20-0027.1.xml

    Parameters 
    ----------
    lifter : nwsspc.sharp.calc.parcel.lifter_lut_cm1
    mix_layer : nwsspc.sharp.calc.layer.PressureLayer
        A mixing layer to take average values of the environment from
    pressure : numpy.ndarray[dtype=float32]
        A 1D NumPy array of pressure values (Pa)
    height : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of height values (Pa)
    temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of temperature values (K)
    mixratio : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of water vapor mixing ratio values (g/g)
    virtual_temperature : numpy.ndarray[dtype=float32] 
        A 1D NumPy array of virtual temperature values (K)
    u_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of U wind component values (m/s)
    v_wind : numpy.ndarray[dtype=float32]
        1D NumPy array of V wind compnent values (m/s)
    potential_temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of potential temperature values (K)
    pcl : None or nwsspc.sharp.calc.parcel.Parcel, optional
        If a parcel is provided, returns the parcel for the PFT
    phi : float 
        The fire moisture to potential temperature increment ratio (kg/(kg*K))
    beta_incr : float 
        The fire plume buoyancy factor (unitless)

    Returns 
    -------
    float 
        The PyroCB Firepower Threshold (Watts)
    """

SBC_MAX_BINS: int = 64

SBC_ICE_NUCLEATION_TEMPERATURE: float = 267.1499938964844

class PrecipType(enum.IntEnum):
    """
    Precipitation type from the spectral bin classifier.

    The seven categories of the 2023 version of the classifier, which adds
    rain mixed with ice pellets to the six categories of Reeves et al. (2016).
    The integer values are stable, for gridded output, and 0 is unused.
    PrecipType.missing is -9999, the value of MISSING.

    References
    ----------
    Reeves et al. 2016: https://doi.org/10.1175/JAMC-D-16-0044.1
    """

    missing = -9999
    """No classification"""

    rain = 1
    """Rain (RA)"""

    snow = 2
    """Snow (SN)"""

    rain_snow = 3
    """Rain and snow (RASN)"""

    freezing_rain = 4
    """Freezing rain (FZRA)"""

    ice_pellets = 5
    """Ice pellets (PL)"""

    freezing_rain_ice_pellets = 6
    """Freezing rain and ice pellets (FZRAPL)"""

    rain_ice_pellets = 7
    """Rain and ice pellets (RAPL)"""

class SpectralBinResult:
    """
    Result of the spectral bin classifier for one column.

    Every field defaults to missing: PrecipType.missing and MISSING. The fields
    are read-only.
    """

    def __init__(self) -> None: ...

    @property
    def precip_type(self) -> PrecipType:
        """Precipitation type at the surface (PrecipType)"""

    @property
    def liquid_fraction(self) -> float:
        """
        Liquid share of the precipitation mass reaching the surface (fraction). It is not a probability. 0.5 means a mix of liquid and ice, not even odds.
        """

    @property
    def supercooled_liquid_height(self) -> float:
        """
        Height of the lowest supercooled liquid water (m AGL). 0 when the surface type is freezing rain, alone or with ice pellets. MISSING when the classifier finds no supercooled liquid.
        """

class SpectralBinDSD:
    """
    A drop-size distribution for the spectral bin classifier.

    Build one with spectral_bin_dsd or spectral_bin_dsd_default, and reuse it
    for every column. It holds up to SBC_MAX_BINS bins, with the per-bin
    constants that the classifier reads. An invalid distribution has
    nbins == 0. The fields are read-only.
    """

    @property
    def nbins(self) -> int:
        """Number of bins, or 0 for an invalid distribution"""

    @property
    def rime_factor(self) -> float:
        """
        Degree of riming of the snow from a frozen cloud top (1 to 5). MISSING for an invalid distribution.
        """

    @property
    def diameter(self) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
        """
        Melted diameter of each bin, exactly as passed (mm). A new float32 array of length nbins.
        """

    @property
    def concentration(self) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
        """
        Number concentration of each bin, exactly as passed. A new float32 array of length nbins.
        """

def spectral_bin_dsd(diameter: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], concentration: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], rime_factor: float = 1.0) -> SpectralBinDSD:
    """
    Builds a drop-size distribution for the spectral bin classifier.

    Validates the bins and precomputes the per-bin constants that the
    classifier reads. The diameters are melted diameters in mm, as in the
    paper and the reference code. Only the ratios of the concentrations
    matter, so any unit works, and a bin may hold 0.

    The result is invalid, with nbins == 0, unless:

    * there are 1 to SBC_MAX_BINS bins
    * the diameters are finite, positive, and strictly increasing
    * the raindrop aspect-ratio fit of the reference is positive at every
      diameter, which holds below about 12.16 mm
    * the concentrations are finite and non-negative, and at least one is
      positive
    * rime_factor is in [1, 5]

    Empty arrays give an invalid distribution.

    rime_factor is the degree of riming of the snow that a frozen cloud top
    produces, from 1 (none) to 5 (graupel). It replaces the reference's fixed
    value of 1. Melting layers below a refreezing layer use 5 instead, as in
    the reference.

    References
    ----------
    Reeves et al. 2016: https://doi.org/10.1175/JAMC-D-16-0044.1

    Python reference (sbc_alg_2023Aug31.py, run_sbc.py): D. Tripp, 2023

    C++ MRMS code (sbcmodel_core.cc): A. Rosenow and D. Tripp

    Parameters
    ----------
    diameter : numpy.ndarray[dtype=float32]
        1D NumPy array of the melted diameter of each bin (mm)
    concentration : numpy.ndarray[dtype=float32]
        1D NumPy array of the number concentration of each bin (any unit)
    rime_factor : float, default = 1.0
        Degree of riming (1 to 5, unitless)

    Returns
    -------
    nwsspc.sharp.calc.params.SpectralBinDSD
        The drop-size distribution, or one with nbins == 0
    """

def spectral_bin_dsd_default() -> SpectralBinDSD:
    """
    The default drop-size distribution of the spectral bin classifier.

    The 4-bin distribution of the 2023 Python reference (run_sbc.py), with
    rime_factor = 1:

    ============== ======= ======= ======= =======
    Diameter (mm)  0.05    0.75    1.45    2.15
    Concentration  55.1843 146.647 11.6891 3.60886
    ============== ======= ======= ======= =======

    The reference interpolates its table of the DSD25 distribution of Reeves
    et al. (2016), 0.05 to 1.65 mm, every 0.7 mm. The 2.15 mm bin lies past
    the end of the table and takes its last value. The paper uses DSD25 as
    measured, with 18 bins 0.1 mm apart and a largest diameter of 1.84 mm.
    The C++ MRMS code (version 2.0.3) caps the bins at 1.85 mm: 0.05, 0.65,
    1.25, and 1.85 mm, with concentrations 55.1843, 206.606, 25.4924, and
    3.60886. Either can be passed to spectral_bin_dsd.

    References
    ----------
    Reeves et al. 2016: https://doi.org/10.1175/JAMC-D-16-0044.1

    Python reference (run_sbc.py): D. Tripp, 2023

    Returns
    -------
    nwsspc.sharp.calc.params.SpectralBinDSD
        The default drop-size distribution
    """
