import pytest
import numpy as np

from nwsspc.sharp.calc import interp
from nwsspc.sharp.calc import constants

def test_interp_height():
    hght = np.arange(100.0, 1100.0, 100.0)
    data = np.arange(1.0, 11.0, 1.0)

    # Test missing/nan/inf behavior
    assert (interp.interp_height(0, hght, data) == constants.MISSING)
    assert (interp.interp_height(1100, hght, data) == constants.MISSING)
    assert (interp.interp_height(constants.MISSING, hght, data) == constants.MISSING)
    assert (interp.interp_height(np.inf, hght, data) == constants.MISSING)
    assert (interp.interp_height(np.nan, hght, data) == constants.MISSING)

    # Test exact values along the edges of the arrays
    assert (interp.interp_height(100, hght, data) == 1)
    assert (interp.interp_height(1000, hght, data) == 10)

    # Test exact value in the middle
    assert (interp.interp_height(500, hght, data) == 5)

    # Test between levels
    assert (interp.interp_height(550, hght, data) == 5.5)
    assert (interp.interp_height(110, hght, data) == pytest.approx(1.1))
    assert (interp.interp_height(391, hght, data) == pytest.approx(3.91))


def test_interp_pres():
    pres = np.arange(10000.0, 110000.0, 10000.0)[::-1]
    data = np.arange(1.0, 11.0, 1.0)

    # Test missing/nan/inf behavior
    assert (interp.interp_pressure(0, pres, data) == constants.MISSING)
    assert (interp.interp_pressure(110000.0, pres, data) == constants.MISSING)
    assert (interp.interp_pressure(constants.MISSING, pres, data) == constants.MISSING)
    assert (interp.interp_pressure(np.inf, pres, data) == constants.MISSING)
    assert (interp.interp_pressure(np.nan, pres, data) == constants.MISSING)

    # Test exact values along the edges of the array
    assert (interp.interp_pressure(100000.0, pres, data) == 1)
    assert (interp.interp_pressure(10000.0, pres, data) == 10)

    # Test an exact value in the middle
    assert (interp.interp_pressure(50000.0, pres, data) == 6)

    # Test between levels
    assert (interp.interp_pressure(97500.0, pres, data)
            == pytest.approx(1.2402969255))
    assert (interp.interp_pressure(95000.0, pres, data)
            == pytest.approx(1.4868382))
    assert (interp.interp_pressure(92500.0, pres, data)
            == pytest.approx(1.73995423))
