interp
======
.. automodule:: nwsspc.sharp.calc.interp 

    QC builds, the default, treat a level whose data is ``constants.MISSING``
    or NaN as missing. ``interp_height`` and ``interp_pressure`` return a
    level's stored value when the requested coordinate lands exactly on a
    level with valid data, even if a neighbouring level is missing. Between
    levels, they interpolate across missing levels to the nearest valid ones,
    and return ``MISSING`` when one side has none. ``find_first_height`` and
    ``find_first_pressure`` skip missing levels, so they still find a crossing
    between the valid levels around them, and an exact match returns that
    level's coordinate even when it is the only valid level.

    .. autofunction:: nwsspc.sharp.calc.interp.interp_height
    .. autofunction:: nwsspc.sharp.calc.interp.interp_pressure
    .. autofunction:: nwsspc.sharp.calc.interp.find_first_height
    .. autofunction:: nwsspc.sharp.calc.interp.find_first_pressure
