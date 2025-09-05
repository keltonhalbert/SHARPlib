winds
=====

.. automodule:: nwsspc.sharp.calc.winds

   Vectors and Components 
   ----------------------
   Representations of atmospheric winds as vectors (speed, angle) or components (u, v), their magnitudes, and the various operatins to convert between representations.

   .. autoclass:: nwsspc.sharp.calc.winds.WindVector
      :members:

   .. autoclass:: nwsspc.sharp.calc.winds.WindComponents
      :members:

   .. autofunction:: nwsspc.sharp.calc.winds.components_to_vector
   .. autofunction:: nwsspc.sharp.calc.winds.vector_to_components
   .. autofunction:: nwsspc.sharp.calc.winds.vector_magnitude
   .. autofunction:: nwsspc.sharp.calc.winds.vector_angle
   .. autofunction:: nwsspc.sharp.calc.winds.u_component
   .. autofunction:: nwsspc.sharp.calc.winds.v_component

   Kinematic Variables
   -------------------
   Calculations derived from wind vectors or components related to mean winds, shear, helicity, vorticity, and more.

   .. autofunction:: nwsspc.sharp.calc.winds.helicity
   .. autofunction:: nwsspc.sharp.calc.winds.wind_shear
   .. autofunction:: nwsspc.sharp.calc.winds.mean_wind
