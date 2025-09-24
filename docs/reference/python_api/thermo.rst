thermo
======

.. automodule:: nwsspc.sharp.calc.thermo

   Dry Variables
   -------------
   Dry thermodynamic variables.

   .. autofunction:: nwsspc.sharp.calc.thermo.theta
   .. autofunction:: nwsspc.sharp.calc.thermo.theta_level
   .. autofunction:: nwsspc.sharp.calc.thermo.lapse_rate
   .. autofunction:: nwsspc.sharp.calc.thermo.lapse_rate_max

   Dry Ascent
   ---------- 
   Dry thermodynamic routines that are directly related to vertical ascent.

   .. autofunction:: nwsspc.sharp.calc.thermo.drylift
   .. autofunction:: nwsspc.sharp.calc.thermo.lcl_temperature

   Moist Variables
   ---------------
   Moist thermodynamic variables.

   .. autofunction:: nwsspc.sharp.calc.thermo.mixratio
   .. autofunction:: nwsspc.sharp.calc.thermo.mixratio_ice
   .. autofunction:: nwsspc.sharp.calc.thermo.vapor_pressure
   .. autofunction:: nwsspc.sharp.calc.thermo.vapor_pressure_ice
   .. autofunction:: nwsspc.sharp.calc.thermo.specific_humidity
   .. autofunction:: nwsspc.sharp.calc.thermo.relative_humidity
   .. autofunction:: nwsspc.sharp.calc.thermo.temperature_at_mixratio
   .. autofunction:: nwsspc.sharp.calc.thermo.virtual_temperature
   .. autofunction:: nwsspc.sharp.calc.thermo.density
   .. autofunction:: nwsspc.sharp.calc.thermo.buoyancy
   .. autofunction:: nwsspc.sharp.calc.thermo.thetae
   .. autofunction:: nwsspc.sharp.calc.thermo.wetbulb
   .. autofunction:: nwsspc.sharp.calc.thermo.theta_wetbulb
   
   Moist Ascent 
   ------------
   Moist thermodunamic routines directly related to vertical ascent. 

   .. autoclass:: nwsspc.sharp.calc.thermo.adiabat 

      An enum class used to flag which type of parcel ascent is desired for nwsspc.sharp.calc.parcel.lifter_cm1

      .. autoattribute:: pseudo_liq

         Pseudoadiabatic liquid-only parcel ascent. 

      .. autoattribute:: pseudo_ice 

         Pseudoadiabatic liquid and ice based parcel ascent.

      .. autoattribute:: adiab_liq

         Adiabatic liquid-only parcel ascent.

      .. autoattribute:: adiab_ice

         Adiabatic liquid and ice based parcel ascent.

   .. autofunction:: nwsspc.sharp.calc.thermo.wobf
   .. autofunction:: nwsspc.sharp.calc.thermo.wetlift

   
