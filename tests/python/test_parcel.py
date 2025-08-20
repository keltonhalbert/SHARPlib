
import os
import pytest
import numpy as np
import pandas as pd

from nwsspc.sharp.calc import layer
from nwsspc.sharp.calc import parcel
from nwsspc.sharp.calc import thermo
from nwsspc.sharp.calc import constants


def load_parquet(filename):
    snd_df = pd.read_parquet(filename)
    snd_df = snd_df[snd_df["vwin"].notna()]
    snd_df = snd_df[snd_df["tmpc"].notna()]
    snd_df = snd_df[snd_df["relh"].notna()]
    snd_df = snd_df[snd_df["pres"] >= 50.0]

    pres = snd_df["pres"].to_numpy().astype('float32')*100.0
    hght = snd_df["hght"].to_numpy().astype('float32')
    tmpk = snd_df["tmpc"].to_numpy().astype('float32')+273.15
    dwpk = snd_df["dwpc"].to_numpy().astype('float32')+273.15
    wdir = snd_df["wdir"].to_numpy().astype('float32')
    wspd = snd_df["wspd"].to_numpy().astype('float32')
    uwin = snd_df["uwin"].to_numpy().astype('float32')
    vwin = snd_df["vwin"].to_numpy().astype('float32')

    # turn into height above ground level
    hght -= hght[0]

    # TO-DO - need a better interface to the API for doing this
    # uwin = np.empty(wspd.shape, dtype="float32")
    # vwin = np.empty(wspd.shape, dtype="float32")
    mixr = thermo.mixratio(pres, dwpk)
    vtmp = thermo.virtual_temperature(tmpk, mixr)

    return {
        "pres": pres, "hght": hght,
        "tmpk": tmpk, "mixr": mixr,
        "vtmp": vtmp, "dwpk": dwpk,
        "wdir": wdir, "wspd": wspd,
        "uwin": uwin, "vwin": vwin
    }


data_dir = os.path.join(
    os.path.dirname(__file__),
    "..",
    "..",
    "data",
    "test_snds"
)
filename = os.path.join(
    data_dir,
    "ddc.parquet"
)
snd_data = load_parquet(filename)


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


def test_surface_parcel():
    pres = snd_data["pres"][0]
    tmpk = snd_data["tmpk"][0]
    dwpk = snd_data["dwpk"][0]

    pcl = parcel.Parcel.surface_parcel(pres, tmpk, dwpk)

    # Wobus Lifter
    lifter = parcel.lifter_wobus()
    vtmpk = pcl.lift_parcel(lifter, snd_data["pres"])
    buoy = thermo.buoyancy(vtmpk, snd_data["vtmp"])
    cape, cinh = pcl.cape_cinh(snd_data["pres"], snd_data["hght"], buoy)

    assert (cape == pytest.approx(3353.6, abs=1e-1))
    assert (cinh == pytest.approx(-34.5697, abs=1e-1))
    assert (pcl.lfc_pressure == pytest.approx(71729, abs=1e-0))
    assert (pcl.eql_pressure == pytest.approx(17933, abs=1e-0))

    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    vtmpk = pcl.lift_parcel(lifter, snd_data["pres"])
    buoy = thermo.buoyancy(vtmpk, snd_data["vtmp"])
    cape, cinh = pcl.cape_cinh(snd_data["pres"], snd_data["hght"], buoy)

    assert (cape == pytest.approx(3107.6, abs=5e-1))
    assert (cinh == pytest.approx(-36.4, abs=5e-1))
    assert (pcl.lfc_pressure == pytest.approx(71482, abs=1e-0))
    assert (pcl.eql_pressure == pytest.approx(18969, abs=1e-0))


def test_mixed_layer_parcel():
    mix_lyr = layer.PressureLayer(
        snd_data["pres"][0], snd_data["pres"][0] - 10000.0)
    theta = thermo.theta(snd_data["pres"], snd_data["tmpk"])
    pcl = parcel.Parcel.mixed_layer_parcel(
        mix_lyr, snd_data["pres"], theta, snd_data["mixr"])

    assert (pcl.pres == pytest.approx(92043.0))
    assert (pcl.tmpk == pytest.approx(297.2543))
    assert (pcl.dwpk == pytest.approx(289.40469))

    lifter = parcel.lifter_wobus()
    vtmpk = pcl.lift_parcel(lifter, snd_data["pres"])
    buoy = thermo.buoyancy(vtmpk, snd_data["vtmp"])
    cape, cinh = pcl.cape_cinh(snd_data["pres"], snd_data["hght"], buoy)

    assert (cape == pytest.approx(2148.9, abs=1e-1))
    assert (cinh == pytest.approx(-128.49, abs=1e-1))
    assert (pcl.lfc_pressure == pytest.approx(67505, abs=1e-0))
    assert (pcl.eql_pressure == pytest.approx(20371, abs=1e-0))

    lifter = parcel.lifter_cm1()
    vtmpk = pcl.lift_parcel(lifter, snd_data["pres"])
    buoy = thermo.buoyancy(vtmpk, snd_data["vtmp"])
    cape, cinh = pcl.cape_cinh(snd_data["pres"], snd_data["hght"], buoy)

    assert (cape == pytest.approx(1929.36, abs=5e-1))
    assert (cinh == pytest.approx(-133.71, abs=5e-1))
    assert (pcl.lfc_pressure == pytest.approx(67100, abs=1e-0))
    assert (pcl.eql_pressure == pytest.approx(20968, abs=1e-0))


def test_most_unstable_parcel():
    search_layer = layer.PressureLayer(
        snd_data["pres"][0],
        snd_data["pres"][0] - 40000.0
    )
    lifter = parcel.lifter_wobus()

    pcl = parcel.Parcel.most_unstable_parcel(
        search_layer,
        lifter,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["vtmp"],
        snd_data["dwpk"]
    )

    assert (pcl.pres == pytest.approx(92043.0))
    assert (pcl.tmpk == pytest.approx(298.15))
    assert (pcl.dwpk == pytest.approx(291.532))
    assert (pcl.cape == pytest.approx(3353.6, abs=1e-1))
    assert (pcl.cinh == pytest.approx(-34.5697, abs=1e-1))
    assert (pcl.lfc_pressure == pytest.approx(71729, abs=1e-0))
    assert (pcl.eql_pressure == pytest.approx(17933, abs=1e-0))

    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq

    pcl = parcel.Parcel.most_unstable_parcel(
        search_layer,
        lifter,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["vtmp"],
        snd_data["dwpk"]
    )

    assert (pcl.pres == pytest.approx(92043.0))
    assert (pcl.tmpk == pytest.approx(298.15))
    assert (pcl.dwpk == pytest.approx(291.532))
    assert (pcl.cape == pytest.approx(3107.6, abs=5e-1))
    assert (pcl.cinh == pytest.approx(-36.4, abs=5e-1))
    assert (pcl.lfc_pressure == pytest.approx(71482, abs=1e0))
    assert (pcl.eql_pressure == pytest.approx(18969, abs=1e0))
