layer
=====

.. automodule:: nwsspc.sharp.calc.layer

    Layer Definitions
    -----------------
    Atmospheric layers can be defined in pressure or height coordinate space. These named structs are used to determine over which layer and in what coordinate space various calculations are performed in. There is an additional named struct that is used to store the index range the layer occupies. 

    .. autoclass:: nwsspc.sharp.calc.layer.PressureLayer
        :members: bottom, top, delta

    .. autoclass:: nwsspc.sharp.calc.layer.HeightLayer
        :members: bottom, top, delta

    .. autoclass:: nwsspc.sharp.calc.layer.LayerIndex
        :members:

    Layer Conversions
    -----------------
    Layers can be converted from one coordinate system to another using interpolation, or have their index values calculated for looping over arrays.

    .. autofunction:: nwsspc.sharp.calc.layer.height_layer_to_pressure
    .. autofunction:: nwsspc.sharp.calc.layer.pressure_layer_to_height
    .. autofunction:: nwsspc.sharp.calc.layer.get_layer_index

    Layer Calculations 
    ------------------
    Perform a search or calculation over a given layer type. 

    .. autofunction:: nwsspc.sharp.calc.layer.layer_min
    .. autofunction:: nwsspc.sharp.calc.layer.layer_max
    .. autofunction:: nwsspc.sharp.calc.layer.layer_mean
    .. autofunction:: nwsspc.sharp.calc.layer.integrate_layer_trapz

