
import pytest
import numpy as np

from nwsspc.sharp.calc import parcel
from nwsspc.sharp.calc import thermo
from nwsspc.sharp.calc import constants


def test_lifter_wobus():
    lifter = parcel.lifter_wobus()

    pres = 100000.0
    tmpk = 320.0
    new_pres = 50000.0

    assert (lifter(constants.MISSING, tmpk, new_pres) == constants.MISSING)
    assert (lifter(pres, constants.MISSING, new_pres) == constants.MISSING)
    assert (lifter(pres, tmpk, constants.MISSING) == constants.MISSING)

    assert (np.isnan(lifter(np.nan, tmpk, new_pres)))
    assert (np.isnan(lifter(pres, np.nan, new_pres)))
    assert (np.isnan(lifter(pres, tmpk, np.nan)))

    assert (lifter(pres, tmpk, new_pres) == pytest.approx(301.872528))


def test_lifter_cm1():
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    lifter.converge = 100.0

    pres = 100000.0
    tmpk = 320.0
    new_pres = 50000.0

    assert (lifter(constants.MISSING, tmpk, new_pres) == constants.MISSING)
    assert (lifter(pres, constants.MISSING, new_pres) == constants.MISSING)
    assert (lifter(pres, tmpk, constants.MISSING) == constants.MISSING)

    assert (np.isnan(lifter(np.nan, tmpk, new_pres)))
    assert (np.isnan(lifter(pres, np.nan, new_pres)))
