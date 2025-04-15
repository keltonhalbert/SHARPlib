import pytest
import numpy as np

from nwsspc.sharp.calc import layer


def test_height_layer_construction():
    lyr1 = layer.HeightLayer(0, 3000)
    assert (lyr1.bottom == 0 and lyr1.top == 3000)

    with pytest.raises(ValueError):
        layer.HeightLayer(3000, 0)

    with pytest.raises(ValueError):
        layer.HeightLayer(0, np.nan)

    with pytest.raises(ValueError):
        layer.HeightLayer(np.nan, 3000)

    with pytest.raises(ValueError):
        layer.HeightLayer(np.nan, np.nan)

    with pytest.raises(ValueError):
        layer.HeightLayer(np.inf, 0)

    with pytest.raises(ValueError):
        layer.HeightLayer(0, np.inf)

    with pytest.raises(ValueError):
        layer.HeightLayer(np.inf, np.inf)

    with pytest.raises(ValueError):
        layer.HeightLayer(np.inf, np.nan)


def test_pressure_layer_construction():
    lyr1 = layer.PressureLayer(100000, 10000)
    assert (lyr1.bottom == 100000 and lyr1.top == 10000)

    with pytest.raises(ValueError):
        layer.PressureLayer(10000, 100000)

    with pytest.raises(ValueError):
        layer.PressureLayer(0, np.nan)

    with pytest.raises(ValueError):
        layer.PressureLayer(np.nan, 3000)

    with pytest.raises(ValueError):
        layer.PressureLayer(np.nan, np.nan)

    with pytest.raises(ValueError):
        layer.PressureLayer(np.inf, 0)

    with pytest.raises(ValueError):
        layer.PressureLayer(0, np.inf)

    with pytest.raises(ValueError):
        layer.PressureLayer(np.inf, np.inf)

    with pytest.raises(ValueError):
        layer.PressureLayer(np.inf, np.nan)


def test_layer_conversion():
    pres = np.arange(10000, 110000, 10000)[::-1]
    hght = np.array([0.0, 500.0, 1500.0, 2500.0, 4000.0,
                    5500.0, 7500.0, 8500.0, 10500.0, 12500.0])

    hght_lyr = layer.HeightLayer(0.0, 3000.0)
    pres_lyr = layer.PressureLayer(100000.0, 75000.0)

    new_pres_lyr = layer.height_layer_to_pressure(hght_lyr, pres, hght)
    new_hght_lyr = layer.pressure_layer_to_height(pres_lyr, pres, hght)

    assert (new_pres_lyr.bottom == 100000.0)
    assert (new_hght_lyr.bottom == 0.0)
    assert (new_pres_lyr.top == pytest.approx(66666.7))
    assert (new_hght_lyr.top == pytest.approx(1983.3248))

    # Check out of bounds behavior
    h_oob1 = layer.HeightLayer(-100, 250.0)
    h_oob2 = layer.HeightLayer(11500.0, 14000.0)
    p_oob1 = layer.PressureLayer(115000.0, 90000.0)
    p_oob2 = layer.PressureLayer(15000.0, 5000.0)

    oob1 = layer.height_layer_to_pressure(h_oob1, pres, hght)
    oob2 = layer.height_layer_to_pressure(h_oob2, pres, hght)
    oob3 = layer.pressure_layer_to_height(p_oob1, pres, hght)
    oob4 = layer.pressure_layer_to_height(p_oob2, pres, hght)

    assert (oob1.bottom == -9999.0)
    assert (oob2.bottom == -9999.0)
    assert (oob3.bottom == -9999.0)
    assert (oob4.bottom == -9999.0)
    assert (oob1.top == -9999.0)
    assert (oob2.top == -9999.0)
    assert (oob3.top == -9999.0)
    assert (oob4.top == -9999.0)

    # Make sure the enforcement of 1D arrays works
    with pytest.raises(TypeError):
        layer.height_layer_to_pressure(hght_lyr, pres.reshape((5, 2)), hght)

    with pytest.raises(TypeError):
        layer.pressure_layer_to_height(pres_lyr, pres.reshape((5, 2)), hght)


def test_layer_bounds_check():
    pres = np.arange(10000, 110000, 10000)[::-1]

    p_oob1 = layer.PressureLayer(100000.0, 5000.0)
    p_oob2 = layer.PressureLayer(110000.0, 10000.0)
    p_oob3 = layer.PressureLayer(110000.0, 5000.0)
    p_oob4 = layer.PressureLayer(5000.0, 2500.0)
    p_oob5 = layer.PressureLayer(115000.0, 110000.0)

    idx1 = layer.get_layer_index(p_oob1, pres)
    idx2 = layer.get_layer_index(p_oob2, pres)
    idx3 = layer.get_layer_index(p_oob3, pres)
    idx4 = layer.get_layer_index(p_oob4, pres)
    idx5 = layer.get_layer_index(p_oob5, pres)

    assert (idx1.kbot == 0)
    assert (idx1.ktop == 9)
    assert (idx2.kbot == 0)
    assert (idx2.ktop == 9)
    assert (idx3.kbot == 0)
    assert (idx3.ktop == 9)
    assert (idx4.kbot == 9)
    assert (idx4.ktop == 9)
    assert (idx5.kbot == 0)
    assert (idx5.ktop == 0)


def test_pressure_layer_min():
    pass


def test_pressure_layer_max():
    pass


def test_height_layer_min():
    pass


def test_height_layer_max():
    pass


def test_pressure_layer_mean():
    pass
