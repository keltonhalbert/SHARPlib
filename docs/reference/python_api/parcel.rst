parcel
======

.. currentmodule:: nwsspc.sharp.calc.parcel

Parcel Lifters
--------------
SHARPlib supports configurable parcel ascent via parcel lifter 'functors'. These can be passed as arguments to lifting routines in order to compute the desired moist adiabatic ascent. 

.. autoclass:: nwsspc.sharp.calc.parcel.lifter_wobus

   .. autofunction:: nwsspc.sharp.calc.parcel.lifter_wobus.__call__
   .. autofunction:: nwsspc.sharp.calc.parcel.lifter_wobus.parcel_virtual_temperature

.. autoclass:: nwsspc.sharp.calc.parcel.lifter_cm1

   .. autoproperty:: nwsspc.sharp.calc.parcel.lifter_cm1.converge
   .. autoproperty:: nwsspc.sharp.calc.parcel.lifter_cm1.ma_type
   .. autoproperty:: nwsspc.sharp.calc.parcel.lifter_cm1.pressure_incr
   .. autofunction:: nwsspc.sharp.calc.parcel.lifter_cm1.__call__
   .. autofunction:: nwsspc.sharp.calc.parcel.lifter_cm1.parcel_virtual_temperature

Parcels
-------
Parcel creation, parcel ascent, and integrated CAPE/CINH.

.. autoclass:: nwsspc.sharp.calc.parcel.LPL

   .. autoattribute:: SFC

      A tag for surface-based parcels. 

   .. autoattribute:: FCST

      A tag for forecast parcels.

   .. autoattribute:: MU

      A tag for most-unstable parcels. 

   .. autoattribute:: ML

      A tag for mixed-layer parcels.

   .. autoattribute:: USR

      A tag for user-defined parcels.


.. autoclass:: nwsspc.sharp.calc.parcel.Parcel 

   .. autoproperty:: nwsspc.sharp.calc.parcel.Parcel.pres
   .. autoproperty:: nwsspc.sharp.calc.parcel.Parcel.tmpk
   .. autoproperty:: nwsspc.sharp.calc.parcel.Parcel.dwpk
   .. autoproperty:: nwsspc.sharp.calc.parcel.Parcel.lcl_pressure
   .. autoproperty:: nwsspc.sharp.calc.parcel.Parcel.lfc_pressure
   .. autoproperty:: nwsspc.sharp.calc.parcel.Parcel.eql_pressure
   .. autoproperty:: nwsspc.sharp.calc.parcel.Parcel.cape
   .. autoproperty:: nwsspc.sharp.calc.parcel.Parcel.cinh

   Parcel Definitions
   ------------------
   .. autofunction:: nwsspc.sharp.calc.parcel.Parcel.surface_parcel
   .. autofunction:: nwsspc.sharp.calc.parcel.Parcel.mixed_layer_parcel
   .. autofunction:: nwsspc.sharp.calc.parcel.Parcel.most_unstable_parcel
   .. autofunction:: nwsspc.sharp.calc.parcel.Parcel.downdraft_parcel
   
   Parcel Operations 
   -----------------
   .. autofunction:: nwsspc.sharp.calc.parcel.Parcel.lift_parcel
   .. autofunction:: nwsspc.sharp.calc.parcel.Parcel.lower_parcel
   .. autofunction:: nwsspc.sharp.calc.parcel.Parcel.find_lfc_el
   .. autofunction:: nwsspc.sharp.calc.parcel.Parcel.cape_cinh

.. autoclass:: nwsspc.sharp.calc.parcel.DowndraftParcel

   .. autoproperty:: nwsspc.sharp.calc.parcel.DowndraftParcel.pres
   .. autoproperty:: nwsspc.sharp.calc.parcel.DowndraftParcel.tmpk
   .. autoproperty:: nwsspc.sharp.calc.parcel.DowndraftParcel.dwpk
   .. autoproperty:: nwsspc.sharp.calc.parcel.DowndraftParcel.cape
   .. autoproperty:: nwsspc.sharp.calc.parcel.DowndraftParcel.cinh

   Parcel Definitions
   ------------------
   .. autofunction:: nwsspc.sharp.calc.parcel.DowndraftParcel.min_thetae

   Parcel Operations
   -----------------
   .. autofunction:: nwsspc.sharp.calc.parcel.DowndraftParcel.lower_parcel
