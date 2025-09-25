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
    muplc : None or nwsspc.sharp.calc.parcel.Parcel, optional

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
    muplc : None or nwsspc.sharp.calc.parcel.Parcel, optional

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
    Estimates supercell storm motion using the Bunkers et al. 2014 
    method described in the following paper:
    http://dx.doi.org/10.15191/nwajom.2014.0211

    This method is parcel based, using a mean-wind vector defined as the 
    pressure-weighted mean wind between the Effective Inflow Layer surface 
    (see effective_inflow_layer routine) and 65% of the depth between that 
    surface and the most unstable parcel's Equilibrium Level. This method 
    produces the same storm motion estimate for surface based supercells, 
    and captures the motion of elevated supercells better than the 
    Bunkers 2000 method. 

    The input parameters of eff_infl_lyr and mupcl (effective inflow layer 
    pressure bounds and the most unstable parcel, respectively) are required
    to be precomputed and passed to this routine. These are expensive 
    operations that are presumed to be computed at some other point 
    in the analysis pipeline. 

    References
    ----------
    Bunkers et al. 2014: http://dx.doi.org/10.15191/nwajom.2014.0211

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

    The effective bulk wind difference is the wind shear between 
    the bottom height of the effective inflow layer, and 50% of 
    the equilibrium level depth. This is analogous to the usage 
    of 0-6 km wind shear, but allows more flexibility for elevated 
    convection. Returns MISSING if the effective inflow layer or 
    equilibrium level pressure are MISSING.
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
    Search for and return the PressureLayer of the lowest altitude 
    hail growth zone. If none is found, the top and bottom of the 
    PressureLayer are set to MISSING. 

    Parameters 
    ----------
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure values (Pa)
    temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of temperature values (K)

    Returns 
    nwsspc.sharp.calc.layer.PressureLayer
        The PressureLayer containing the hail growth zone
    """

def dendritic_layer(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> nwsspc.sharp.calc.layer.PressureLayer:
    """
    Search for and return the PressureLayer of the lowest altitude 
    dendritic growth zone. If none is found, the top and bottom of the 
    PressureLayer are set to MISSING. 

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

def rainfall_velocity(rho: float, rho_0: float, rainwater_mixratio: float) -> float:
    """
    Computes the estimated rainfall terminal velocity. 

    Computes the estimated rainfall terminal velocity as formulated 
    by Klemp and Wilhelmson 1978. 

    Parameters
    ----------
    rho : float 
        Air density (kg/m^3)
    rho_0 : float 
        Air density at the surface (kg/m^3)
    rainwater_mixratio : float
        Rainwater mixing ratio (kg/kg)

    Returns 
    -------
    float
        Rainfall velocity (m/s)
    """

def rainfall_efficiency(pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixr: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], pcl: nwsspc.sharp.calc.parcel.Parcel, rainwater_mixratio: float) -> float:
    """
    Computes an estimated rainfall efficiency.

    Computes an estimated rainfall efficiency, using formulations of 
    rainfall terminal velocity and evaporation rate as in Klemp and 
    Wilhelmson 1978.

    Starting from the parcel's lifted condensation level, it integrates 
    downward to the surface while evaporating water content and updating 
    the velocity between steps. Once it reaches the surface, it computes 
    the final fraction of rainwater mixing ratio relative to the starting 
    rainwater mixing ratio.

    Parameters
    ----------
    pressure : numpy.ndarray[dtype=float32]
        1D NumPy array of pressure (Pa)
    height : numpy.ndarray[dtype=float32]
        1D NumPy array of height (Pa)
    temperature : numpy.ndarray[dtype=float32]
        1D NumPy array of temperature (K)
    mixr : numpy.ndarray[dtype=float32]
        1D NumPy array of water vapor mixing ratio (kg/kg)
    pcl : nwsspc.sharp.calc.parcel.Parcel 
        A Parcel with a defined LCL pressure
    rainwater_mixratio : float
        Rainwater mixing ratio (kg/kg)

    Returns
    -------
    float
        The fraction of precipitation that reaches the surface.
    """
