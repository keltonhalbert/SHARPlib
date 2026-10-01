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
    assert (interp.interp_height(constants.MISSING,
            hght, data) == constants.MISSING)
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
    assert (interp.interp_pressure(
        constants.MISSING, pres, data) == constants.MISSING)
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


def test_find_first_pres():
    pres = np.arange(10000.0, 110000.0, 10000.0)[::-1]
    data = np.arange(1.0, 11.0, 1.0)
    data[2] = constants.MISSING

    assert (interp.find_first_pressure(5.0, pres, data) == 60000.0)
    assert (interp.find_first_pressure(
        5.5, pres, data) == pytest.approx(54772.26))
    assert (interp.find_first_pressure(5.0, pres, data[::-1]) == 50000.0)
    assert (interp.find_first_pressure(
        5.5, pres, data[::-1]) == pytest.approx(54772.26))

    assert (interp.find_first_pressure(
        constants.MISSING, pres, data) == constants.MISSING)
    assert (interp.find_first_pressure(
        3, pres, data) == pytest.approx(79372.6))


def test_find_first_height():
    hght = np.arange(100.0, 1100.0, 100.0)
    data = np.arange(1.0, 11.0, 1.0)
    data[2] = constants.MISSING

    assert (interp.find_first_height(5.0, hght, data) == 500.0)
    assert (interp.find_first_height(5.5, hght, data) == 550.0)
    assert (interp.find_first_height(5.0, hght, data[::-1]) == 600.0)
    assert (interp.find_first_height(5.5, hght, data[::-1]) == 550.0)

    assert (interp.find_first_height(
        constants.MISSING, hght, data) == constants.MISSING)
    assert (interp.find_first_height(3, hght, data) == 300.0)


# QC builds handle NaN data like MISSING, and a query that lands exactly on a
# valid level returns that level's stored value. Each "was" comment is the
# output measured before that change.
def test_interp_height_nan_and_missing():
    hght = np.array([0.0, 100.0, 200.0], dtype="float32")

    # NaN neighbour is bridged; no valid level on one side gives MISSING
    data = np.array([1.0, np.nan, 9.0], dtype="float32")
    assert (interp.interp_height(150.0, hght, data) == 7.0)  # was 9
    assert (interp.interp_height(
        50.0, hght[:2], data[:2]) == constants.MISSING)  # was NaN

    # exact level next to a missing level
    for x in (np.nan, constants.MISSING):
        top = np.array([x, 280.0], dtype="float32")
        bot = np.array([280.0, x], dtype="float32")
        # was 280 (NaN) / MISSING (MISSING)
        assert (interp.interp_height(100.0, hght[:2], top) == 280.0)
        # was NaN (NaN) / MISSING (MISSING)
        assert (interp.interp_height(0.0, hght[:2], bot) == 280.0)

    hght_td = np.array([2950.0, 3000.0, 3050.0], dtype="float32")
    dwpk = np.array([270.0, 269.5, constants.MISSING], dtype="float32")
    assert (interp.interp_height(3000.0, hght_td, dwpk) == 269.5)  # was MISSING


def test_interp_pressure_nan_and_missing():
    pres = np.array([100000.0, 90000.0, 80000.0], dtype="float32")

    # NaN neighbour is bridged; no valid level on one side gives MISSING
    data = np.array([1.0, np.nan, 9.0], dtype="float32")
    assert (interp.interp_pressure(85000.0, pres, data)
            == pytest.approx(6.82651615))  # was 9
    assert (interp.interp_pressure(
        95000.0, pres[:2], data[:2]) == constants.MISSING)  # was NaN

    # exact level next to a missing level
    for x in (np.nan, constants.MISSING):
        top = np.array([x, 280.0], dtype="float32")
        bot = np.array([280.0, x], dtype="float32")
        # was 280 (NaN) / MISSING (MISSING)
        assert (interp.interp_pressure(90000.0, pres[:2], top) == 280.0)
        # was NaN (NaN) / MISSING (MISSING)
        assert (interp.interp_pressure(100000.0, pres[:2], bot) == 280.0)

    dwpk = np.array([270.0, 269.5, constants.MISSING], dtype="float32")
    assert (interp.interp_pressure(90000.0, pres, dwpk) == 269.5)  # was MISSING


def test_find_first_nan_and_missing():
    hght = np.array([0.0, 100.0, 200.0], dtype="float32")
    pres = np.array([100000.0, 90000.0, 80000.0], dtype="float32")

    # a crossing across a NaN level is found (was MISSING)
    data = np.array([1.0, np.nan, 9.0], dtype="float32")
    assert (interp.find_first_height(5.0, hght, data) == 100.0)
    assert (interp.find_first_pressure(5.0, pres, data)
            == pytest.approx(89442.67))

    # exact match on the only valid level
    # was the coordinate for NaN, MISSING for MISSING
    for x in (np.nan, constants.MISSING):
        for levels, idx in (([5.0, x], 0), ([x, 5.0], 1), ([x, 5.0, x], 1)):
            data = np.array(levels, dtype="float32")
            n = len(data)
            assert (interp.find_first_height(
                5.0, hght[:n], data) == hght[idx])
            assert (interp.find_first_pressure(
                5.0, pres[:n], data) == pres[idx])

    # single-level profile: an exact match returns its coordinate
    # (was MISSING)
    data = np.array([5.0], dtype="float32")
    assert (interp.find_first_height(5.0, hght[:1], data) == 0.0)
    assert (interp.find_first_pressure(5.0, pres[:1], data) == 100000.0)
    assert (interp.find_first_height(
        6.0, hght[:1], data) == constants.MISSING)
