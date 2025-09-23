import enum
from typing import Annotated, overload

import numpy
from numpy.typing import NDArray

import nwsspc.sharp.calc.layer
import nwsspc.sharp.calc.thermo


class lifter_wobus:
    """
    A functor that calls the Wobus Wetlift for computation of moist adiabats.

    This is used to wrap the Wobus Wetlift function for parcel lifting 
    routines. Functors -- classes with their operator() overloaded -- 
    are used so that functions can be passed to templates in a way that 
    the compiler can still optimize, rather than using function pointers
    or lambdas. 

    Specifically, this functor is designed to be passed as a template 
    argument to nwsspc.sharp.calc.parcel.Parcel.lift_parcel, so that 
    the method of computing moist adiabats can be changed without 
    changing the overall parcel lifting code. The reason this is awesome 
    is that the compiler can still optimize and inline this code, while 
    the user can configure the parcel lifting algorithm to their 
    specifications.
    """

    def __init__(self) -> None: ...

    lift_from_lcl: bool = ...
    """
    A static flag that helps the parcel lifting functions know where to lift from.
    """

    @property
    def converge(self) -> float:
        """The iterative convergence criteria (K)"""

    @converge.setter
    def converge(self, arg: float, /) -> None: ...

    def setup(self, lcl_pres: float, lcl_tmpk: float) -> None:
        """
        Some parcel lifters require setup in order to handle adiabatic ascent
        when tracking the water vapor, liquid, and ice mixing ratios. The Wobus 
        lifter does not require this, however, so this function does nothing.
        """

    def __call__(self, pres: float, tmpk: float, new_pres: float) -> float:
        """
        Overloads the call operator in order to facilitate parcel lifting. 

        Parameters
        ----------
        pres : float 
            Parcel pressure (Pa)
        tmpk : float 
            Parcel temperature (K)
        new_pres : float 
            Final level of parcel after lift (Pa)

        Returns
        -------
        float
            The temperature of the lifted parcel (K)
        """

    def parcel_virtual_temperature(self, pres: float, tmpk: float) -> float:
        """
        Computes the virtual temperature of the parcel (after saturation).

        Parameters
        ----------
        pres : float 
            Parcel pressure (Pa)
        tmpk : float 
            Parcel temperature (K)

        Returns
        -------
        float
            The virtual temperature of the parcel (K)
        """

class lifter_cm1:
    """
    Use the CM1 moist lift calculations for adiabatic and pseudoadiabatic
    parcel ascent.
    """

    def __init__(self) -> None: ...

    lift_from_lcl: bool = ...
    """
    A static flag that helps the parcel lifting functions know where to lift from.
    The lifter_cm1 lifts from the last lifted level, rather than the LCL, because
    it is an iterative solver. This results in major performance improvements while 
    maintaining accuracy.
    """

    @property
    def ma_type(self) -> nwsspc.sharp.calc.thermo.adiabat:
        """
        The type of moist adiabat to use, defined by nwsspc.sharp.calc.thermo.adiabat
        """

    @ma_type.setter
    def ma_type(self, arg: nwsspc.sharp.calc.thermo.adiabat, /) -> None: ...

    @property
    def pressure_incr(self) -> float:
        """The pressure increment (Pa) to use for the iterative solver."""

    @pressure_incr.setter
    def pressure_incr(self, arg: float, /) -> None: ...

    @property
    def converge(self) -> float:
        """The iterative convergence criteria (K)"""

    @converge.setter
    def converge(self, arg: float, /) -> None: ...

    def setup(self, lcl_pres: float, lcl_tmpk: float) -> None:
        """
        This function sets the total water mixing ratio for 
        adiabatic parcel ascent, and zeroes out the vapor,
        liquid, and ice mixing ratio from previous parcel 
        ascents. 

        Parameters
        ----------
        lcl_pres : float 
            The LCL pressure (Pa)
        lcl_tmpk : float 
            The LCL temperature (K)

        Returns
        -------
        None
        """

    def __call__(self, pres: float, tmpk: float, new_pres: float) -> float:
        """
        Lifts a parcel moist adiabatically/pseudoadiabatically using
        nwsspc.sharp.calc.thermo.moist_adiabat_cm1.

        Parameters
        ----------
        pres : float 
            Parcel pressure (Pa)
        tmpk : float 
            Parcel temperature (K)
        new_pres : float 
            Final level of parcel after lift (Pa)

        Returns
        -------
        float
            The temperature of the lifted parcel (K)
        """

    def parcel_virtual_temperature(self, pres: float, tmpk: float) -> float:
        """
        Computes the virtual temperature of the parcel (after saturation).

        Parameters
        ----------
        pres : float 
            Parcel pressure (Pa)
        tmpk : float 
            Parcel temperature (K)

        Returns
        -------
        float
            The virtual temperature of the parcel (K)
        """

class LPL(enum.Enum):
    SFC = 1
    """A Surface Based Parcel"""

    FCST = 2
    """A Forecast Surface Parcel"""

    MU = 3
    """A Most Unstable Parcel"""

    ML = 4
    """A Mixed-Layer Parcel"""

    USR = 5
    """A User Defined Parcel"""

class Parcel:
    """
    Contains information about a Parcel's starting level and
    thermodynamic attributes, as well as derived computations,
    methods for constructing a parcel, and parcel ascent routines.
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, pressure: float, temperature: float, dewpoint: float, lpl: LPL) -> None:
        """
        Constructor for a Parcel

        Parameters
        ----------
        pressure : float 
            Parcel initial pressure (Pa)
        temperature : float 
            Parcel initial temperature (K)
        dewpoint : float 
            Parcel initial dewpoint (K)
        lpl : nwsspc.sharp.calc.parcel.LPL 
            Parcel Lifted Parcel Level (LPL) definition
        """

    @property
    def pres(self) -> float:
        """Parcel starting pressure (Pa)"""

    @pres.setter
    def pres(self, arg: float, /) -> None: ...

    @property
    def tmpk(self) -> float:
        """Parcel starting temperature (K)"""

    @tmpk.setter
    def tmpk(self, arg: float, /) -> None: ...

    @property
    def dwpk(self) -> float:
        """Parcel starting dewpoint (K)"""

    @dwpk.setter
    def dwpk(self, arg: float, /) -> None: ...

    @property
    def lcl_pressure(self) -> float:
        """Pressure at the Lifted Condensation Level (Pa)"""

    @lcl_pressure.setter
    def lcl_pressure(self, arg: float, /) -> None: ...

    @property
    def lfc_pressure(self) -> float:
        """Pressure at the Level of Free Convection (Pa)"""

    @lfc_pressure.setter
    def lfc_pressure(self, arg: float, /) -> None: ...

    @property
    def eql_pressure(self) -> float:
        """Pressure at the parcel Equilibrium Level"""

    @eql_pressure.setter
    def eql_pressure(self, arg: float, /) -> None: ...

    @property
    def mpl_pressure(self) -> float:
        """Pressure at the Maximum Parcel Level"""

    @mpl_pressure.setter
    def mpl_pressure(self, arg: float, /) -> None: ...

    @property
    def cape(self) -> float:
        """
        Parcel Convective Available Potential Energy (J/kg) between the LFC and EL
        """

    @cape.setter
    def cape(self, arg: float, /) -> None: ...

    @property
    def cinh(self) -> float:
        """Parcel Convective Inhibition (J/kg) between the LFC and EL"""

    @cinh.setter
    def cinh(self, arg: float, /) -> None: ...

    @overload
    def lift_parcel(self, lifter: lifter_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
        """
        Lifts a Parcel dry adiabatically from its LPL to its LCL dry
        adiabatically, and then moist adiabatically from the LCL to 
        the top of the profile. The moist adiabat used is determined
        by the type of lifting functor passed to the function (i.e.
        lifter_wobus or lifter_cm1).

        Parameters
        ----------
        lifter : nwsspc.sharp.calc.parcel.lifter_wobus 
            An instantiated lifter_wobus functor
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of Pressure levels for lifting (Pa)

        Returns
        -------
        numpy.ndarray[dtype=float32]
            A 1D NumPy array of parcel virtual temperature values (K)
        """

    @overload
    def lift_parcel(self, lifter: lifter_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
        """
        Lifts a Parcel dry adiabatically from its LPL to its LCL dry
        adiabatically, and then moist adiabatically from the LCL to 
        the top of the profile. The moist adiabat used is determined
        by the type of lifting functor passed to the function (i.e.
        lifter_wobus or lifter_cm1).

        Parameters
        ----------
        lifter : nwsspc.sharp.calc.parcel.lifter_cm1 
            An instantiated lifter_cm1 functor
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of Pressure levels for lifting (Pa)

        Returns
        -------
        numpy.ndarray[dtype=float32]
            A 1D NumPy array of parcel virtual temperature values (K)
        """

    def find_lfc_el(self, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], buoyancy: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[float, float]:
        """
        Searches the buoyancy array for the LFC and EL combination that results
        in the maximum amount of CAPE in the given profile. The buoyancy array 
        is typically computed by calling nwsspc.sharp.calc.parcel.Parcel.lift_parcel. 
        Once the LFC and EL are found, the value are set in 
        nwsspc.sharp.calc.parcel.Parcel.lfc_pres and 
        nwsspc.sharp.calc.parcel.Parcel.eql_pres via the provided parcel.

        The value of eql_pres is MISSING if there is no qualifying level 
        found within the data bounds (e.g. incomplete data, or EL above 
        the available data). Any calls to nwsspc.sharp.calc.parcel.Parcel.cape_cinh 
        will still compute CAPE without the presence of an EL, using the best-available
        data. 

        Parameters
        ----------
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of pressure values (Pa)
        height : numpy.ndarray[dtype=float32]
            1D NumPy array of height values (meters)
        buoyancy : numpy.ndarray[dtype=float32] 
            1D NumPy array of buoyancy values (m/s^2)

        Returns
        -------
        tuple[float, float]
            (LFC_PRES, EL_PRES)
        """

    def maximum_parcel_level(self, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], buoyancy: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> float:
        """
        Find the pressure of the Maximum Parcel Level (MPL).

        The Maximum Parcel Level (MPL) is the level a parcel woud reach 
        if it expended all of its integrated positive buoyancy past the 
        Equilibrium Level. It is found by integrating negatively buoyant 
        area above the Equilibrium Level until the integrated negative 
        buoyancy is equal in magnitude to the Convective Available 
        Potential Energy between the Level of Free Convection and the 
        Equilibrium Level. 

        For valid calculations, nwsspc.sharp.calc.parcel.Parcel.cape_cinh 
        must be called first, or nwsspc.sharp.calc.parcel.Parcel.cape and 
        nwsspc.sharp.calc.parcel.Parcel.eql_pressure must be set. 

        A values of nwsspc.sharp.calc.constants.MISSING is returned if:
          * CAPE is 0 
          * nwsspc.sharp.calc.parce.Parcel.eql_pressure is MISSING
          * No valid MPL candidate is found within the profile
            In this scenario, it likely exceeds the top of the available data

        In addition to being returned, the result is stored inside of 
        nwsspc.sharp.calc.parcel.Parcel.mpl_pressure.

        Parameters
        ----------
        pres : numpy.ndarray[dtype=float32] 
            1D NumPy array of pressure values (Pa)
        hght : numpy.ndarray[dtype=float32] 
            1D NumPy array of height values (meters)
        buoyancy : numpy.ndarray[dtype=float32] 
            1D NumPy array of buoyancy values (m/s^2)

        Returns 
        -------
        float 
            The pressure of the Maximum Parcel Level (Pa)
        """

    def cape_cinh(self, pres: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], hght: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], buoy: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[float, float]:
        """
        Assuming that nwsspc.sharp.calc.parcel.Parcel.lift_parcel has 
        been called, cape_cinh will integrate the area between the LFC 
        and EL to compute CAPE, and integrate the area between the LPL 
        and LCL to compute CINH.

        If eql_pressure is MISSING, but lfc_pressure is defined, the 
        routine will compute CAPE with the available data despite the 
        lack of a defined equilibrium level. This is useful for 
        incomplete profile data, or pressure-level data where the EL 
        is above the top pressure value.

        The results are stored in nwsspc.sharp.calc.parcel.Parcel.cape 
        and nwsspc.sharp.calc.parcel.Parcel.cinh via the provided parcel.

        Parameters
        ----------
        pres : numpy.ndarray[dtype=float32] 
            1D NumPy array of pressure values (Pa)
        hght : numpy.ndarray[dtype=float32] 
            1D NumPy array of height values (meters)
        buoyancy : numpy.ndarray[dtype=float32] 
            1D NumPy array of buoyancy values (m/s^2)

        Returns
        -------
        tuple[float, float]
            (CAPE, CINH)
        """

    @staticmethod
    def surface_parcel(pressure: float, temperature: float, dewpoint: float) -> Parcel:
        """
        Given input values of surface pressure, temperature, and dewpoint 
        temperature, construct and return a Surface Based Parcel.

        Parameters
        ----------
        pressure : float 
            Surface pressure (Pa)
        temperature : float 
            Surface temperature (K)
        dewpoint : float 
            Surface dewpoint (K)

        Returns
        -------
        nwsspc.sharp.calc.parcel.Parcel
            Parcel with surface values
        """

    @overload
    @staticmethod
    def mixed_layer_parcel(mix_layer: nwsspc.sharp.calc.layer.PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], potential_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixing_ratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Parcel:
        """
        Given input arrays of pressure, potential temperature, 
        and water vapor mixing ratio, as well as a defined PressureLayer, 
        compute and return a mixed-layer Parcel.

        Parameters
        ----------
        mix_layer : nwsspc.sharp.calc.layer.PressureLayer 
            PressureLayer over which to compute a mixed-layer parcel 
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile pressure values (Pa)
        potential_temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile potential temperature (K)
        mixing_ratio : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile water vapor mixing ratio (unitless)

        Returns
        -------
        nwsspc.sharp.calc.parcel.Parcel
            Parcel with mixed layer values
        """

    @overload
    @staticmethod
    def mixed_layer_parcel(mix_layer: nwsspc.sharp.calc.layer.HeightLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], potential_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mixing_ratio: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Parcel:
        """
        Given input arrays of pressure, potential temperature, and water 
        vapor mixing ratio, as well as a defined PressureLayer, compute 
        and return a mixed-layer Parcel.

        Parameters
        ----------
        mix_layer : nwsspc.sharp.calc.layer.HeightLayer 
            HeightLayer over which to compute a mixed-layer parcel 
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile pressure values (Pa)
        height : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile height values (meters)
        potential_temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile potential temperature (K)
        mixing_ratio : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile water vapor mixing ratio (unitless)

        Returns
        -------
        nwsspc.sharp.calc.parcel.Parcel
            Parcel with mixed layer values
        """

    @overload
    @staticmethod
    def most_unstable_parcel(layer: nwsspc.sharp.calc.layer.PressureLayer, lifter: lifter_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Parcel:
        """
        Given input arrays of pressure, height, temperature, virtual temperature,
        and dewpoint temperature, as well as a defined PressureLayer/HeightLayer and 
        parcel lifter (lifter_wobus or lifter_cm1), find and return the most unstable parcel. 

        Parameters
        ----------
        layer : nwsspc.sharp.calc.layer.PressureLayer 
            PressureLayer for which to search for the Most Unstable Parcel
        lifter : nwsspc.sharp.calc.parcel.lifter_cm1 
            Parcel lifting routine to use for moist ascent
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile pressure values (Pa)
        height : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile height values (meters)
        temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile temperature values (K)
        virtual_temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile virtual temperature values (K)
        dewpoint : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile dewpoint values (K)

        Returns
        -------
        nwsspc.sharp.calc.parcel.Parcel
            Parcel with most-unstable values
        """

    @overload
    @staticmethod
    def most_unstable_parcel(layer: nwsspc.sharp.calc.layer.HeightLayer, lifter: lifter_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Parcel:
        """
        Given input arrays of pressure, height, temperature, virtual temperature,
        and dewpoint temperature, as well as a defined PressureLayer/HeightLayer and 
        parcel lifter (lifter_wobus or lifter_cm1), find and return the most unstable parcel. 

        Parameters
        ----------
        layer : nwsspc.sharp.calc.layer.HeightLayer 
            HeightLayer for which to search for the Most Unstable Parcel
        lifter : nwsspc.sharp.calc.parcel.lifter_cm1 
            Parcel lifting routine to use for moist ascent
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile pressure values (Pa)
        height : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile height values (meters)
        temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile temperature values (K)
        virtual_temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile virtual temperature values (K)
        dewpoint : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile dewpoint values (K)

        Returns
        -------
        nwsspc.sharp.calc.parcel.Parcel
            Parcel with most-unstable values
        """

    @overload
    @staticmethod
    def most_unstable_parcel(layer: nwsspc.sharp.calc.layer.PressureLayer, lifter: lifter_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Parcel:
        """
        Given input arrays of pressure, height, temperature, virtual temperature,
        and dewpoint temperature, as well as a defined PressureLayer/HeightLayer and 
        parcel lifter (lifter_wobus or lifter_cm1), find and return the most unstable parcel. 

        Parameters
        ----------
        layer : nwsspc.sharp.calc.layer.PressureLayer 
            PressureLayer for which to search for the Most Unstable Parcel
        lifter : nwsspc.sharp.calc.parcel.lifter_wobus 
            Parcel lifting routine to use for moist ascent
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile pressure values (Pa)
        height : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile height values (meters)
        temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile temperature values (K)
        virtual_temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile virtual temperature values (K)
        dewpoint : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile dewpoint values (K)

        Returns
        -------
        nwsspc.sharp.calc.parcel.Parcel
        """

    @overload
    @staticmethod
    def most_unstable_parcel(layer: nwsspc.sharp.calc.layer.HeightLayer, lifter: lifter_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], height: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], virtual_temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Parcel:
        """
        Given input arrays of pressure, height, temperature, virtual temperature,
        and dewpoint temperature, as well as a defined PressureLayer/HeightLayer and 
        parcel lifter (lifter_wobus or lifter_cm1), find and return the most unstable parcel. 

        Parameters
        ----------
        layer : nwsspc.sharp.calc.layer.HeightLayer 
            HeightLayer for which to search for the Most Unstable Parcel
        lifter : nwsspc.sharp.calc.parcel.lifter_wobus 
            Parcel lifting routine to use for moist ascent
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile pressure values (Pa)
        height : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile height values (meters)
        temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile temperature values (K)
        virtual_temperature : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile virtual temperature values (K)
        dewpoint : numpy.ndarray[dtype=float32] 
            1D NumPy array of profile dewpoint values (K)

        Returns
        -------
        nwsspc.sharp.calc.parcel.Parcel
            Parcel with most-unstable values
        """

class DowndraftParcel:
    """
    Contains information about a DowndraftParcel's starting level and
    thermodynamic attributes, as well as derived computations,
    methods for constructing a parcel, and parcel descent routines.
    """

    @overload
    def __init__(self) -> None: ...

    @overload
    def __init__(self, pressure: float, temperature: float, dewpoint: float) -> None:
        """
        Constructor for a DowndraftParcel

        Parameters
        ----------
        pressure : float 
            DowndraftParcel initial pressure (Pa)
        temperature : float 
            DowndraftParcel initial temperature (K)
        dewpoint : float 
            DowndraftParcel initial dewpoint (K)
        """

    @property
    def pres(self) -> float:
        """DowndraftParcel starting pressure (Pa)"""

    @pres.setter
    def pres(self, arg: float, /) -> None: ...

    @property
    def tmpk(self) -> float:
        """DowndraftParcel starting temperature (K)"""

    @tmpk.setter
    def tmpk(self, arg: float, /) -> None: ...

    @property
    def dwpk(self) -> float:
        """DowndraftParcel starting dewpoint (K)"""

    @dwpk.setter
    def dwpk(self, arg: float, /) -> None: ...

    @property
    def cape(self) -> float:
        """DowndraftParcel Convective Available Potential Energy (J/kg)"""

    @cape.setter
    def cape(self, arg: float, /) -> None: ...

    @property
    def cinh(self) -> float:
        """DowndraftParcel Convective Inhibition (J/kg)"""

    @cinh.setter
    def cinh(self, arg: float, /) -> None: ...

    @overload
    def lower_parcel(self, lifter: lifter_wobus, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
        """
        Lowers a saturated nwsspc.sharp.calc.parcel.DowndraftParcel moist 
        adiabatically from its LPL to the surface. The moist adiabat used 
        is determined by the type of lifting functor passed to the function 
        (i.e. lifter_wobus or lifter_cm1).

        Unlike nwsspc.sharp.calc.parcel.Parcel.lift_parcel, the virtual 
        temperature correction is not used for downdraft parcels.

        Parameters
        ----------
        lifter : nwsspc.sharp.calc.parcel.lifter_wobus 
            An instantiated lifter_wobus functor
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of Pressure levels for lifting (Pa)

        Returns
        -------
        numpy.ndarray[dtype=float32]
            A 1D NumPy array of parcel temperature values (K)
        """

    @overload
    def lower_parcel(self, lifter: lifter_cm1, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C')]:
        """
        Lowers a saturated nwsspc.sharp.calc.parcel.DowndraftParcel moist 
        adiabatically from its LPL to the surface. The moist adiabat used 
        is determined by the type of lifting functor passed to the function 
        (i.e. lifter_wobus or lifter_cm1).

        Unlike nwsspc.sharp.calc.parcel.Parcel.lift_parcel, the virtual 
        temperature correction is not used for downdraft parcels.

        Parameters
        ----------
        lifter : nwsspc.sharp.calc.parcel.lifter_cm1
            An instantiated lifter_cm1 functor
        pressure : numpy.ndarray[dtype=float32] 
            1D NumPy array of Pressure levels for lifting (Pa)

        Returns
        -------
        numpy.ndarray[dtype=float32]
            A 1D NumPy array of parcel temperature values (K)
        """

    def cape_cinh(self, pres: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], hght: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], buoy: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)]) -> tuple[float, float]:
        """
        Assuming that nwsspc.sharp.calc.parcel.DowndraftParcel.lower_parcel 
        has been called, cape_cinh will integrate the area between the LPL
        and the surface to compute downdraft CAPE and downdraft CINH.

        The results are stored in nwsspc.sharp.calc.parcel.DowndraftParcel.cape 
        and nwsspc.sharp.calc.parcel.DowndraftParcel.cinh via the provided parcel.

        Parameters
        ----------
        pres : numpy.ndarray[dtype=float32] 
            1D NumPy array of pressure values (Pa)
        hght : numpy.ndarray[dtype=float32] 
            1D NumPy array of height values (meters)
        buoyancy : numpy.ndarray[dtype=float32] 
            1D NumPy array of buoyancy values (m/s^2)

        Returns
        -------
        tuple[float, float]
            (DCAPE, DCINH)
        """

    @staticmethod
    def min_thetae(search_layer: nwsspc.sharp.calc.layer.PressureLayer, pressure: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], temperature: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], dewpoint: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], thetae: Annotated[NDArray[numpy.float32], dict(shape=(None,), order='C', device='cpu', writable=False)], mean_depth: float = 10000.0) -> DowndraftParcel:
        """
        Define a downdraft parcel. 

        Defines a downdraft parcel within a given search layer. 
        The downdraft parcel is defined as the minimum layer-mean 
        equivalent potential temperature (Theta-E) within the 
        search layer. Typical values are to search within the lowest
        400 hPa of the profile, and a mean depth of 100 hPa. 

        Parameters 
        ----------
        search_layer : nwsspc.sharp.calc.layer.PressureLayer 
            The layer over which to search for the downdraft parcel 
        pressure : numpy.ndarray[dtype=float32]
            1D NumPy array of pressure (Pa)
        temperature : numpy.ndarray[dtype=float32]
            1D NumPy array of temperature (K)
        dewpoint : numpy.ndarray[dtype=float32]
            1D NumPy array of dewpoint (K)
        thetae : numpy.ndarray[dtype=float32]
            1D NumPy array of thetae (K)
        mean_depth : float
            The layer depth for calculating mean thetae.

        Returns 
        -------
        nwsspc.sharp.calc.parcel.DowndraftParcel 
            Downdraft Parcel
        """
