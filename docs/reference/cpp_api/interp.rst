interp
======

QC builds, the default, treat a level whose data is ``sharp::MISSING`` or NaN
as missing. ``interp_height`` and ``interp_pressure`` return a level's stored
value when the requested coordinate lands exactly on a level with valid data,
even if a neighbouring level is missing. Between levels, they interpolate
across missing levels to the nearest valid ones, and return ``MISSING`` when
one side has none. ``find_first_height`` and ``find_first_pressure`` skip
missing levels, so they still find a crossing between the valid levels around
them, and an exact match returns that level's coordinate even when it is the
only valid level. Builds configured with ``-DNO_QC=ON`` skip these checks.

.. doxygenfunction:: sharp::interp_pressure
.. doxygenfunction:: sharp::interp_height 
.. doxygenfunction:: sharp::find_first_pressure
.. doxygenfunction:: sharp::find_first_height 
