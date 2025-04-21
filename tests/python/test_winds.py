

import os
import pytest
import numpy as np
import pandas as pd

from nwsspc.sharp.calc import layer
from nwsspc.sharp.calc import winds
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


def test_helicity():
    hlyr = layer.HeightLayer(0, 3000.0)
    storm_mtn = winds.WindComponents()
    storm_mtn.u = 0
    storm_mtn.v = 0

    srh = winds.helicity(
        hlyr, storm_mtn, snd_data["hght"], snd_data["uwin"], snd_data["vwin"])

    assert (srh == pytest.approx(169.71))


def test_mean_wind():

    plyr = layer.PressureLayer(snd_data["pres"][0], 60000.0)
    mw = winds.mean_wind(
        plyr, snd_data["pres"], snd_data["uwin"], snd_data["vwin"])

    assert (mw.u == pytest.approx(4.39905))
    assert (mw.v == pytest.approx(11.97378))

    mw = winds.mean_wind(
        plyr, snd_data["pres"], snd_data["uwin"], snd_data["vwin"], True)

    assert (mw.u == pytest.approx(3.683177))
    assert (mw.v == pytest.approx(11.69769))


def test_wind_shear():
    hlyr = layer.HeightLayer(0, 6000.0)
    shr = winds.wind_shear(
        hlyr, snd_data["hght"], snd_data["uwin"], snd_data["vwin"])

    assert (shr.u == pytest.approx(14.6))
    assert (shr.v == pytest.approx(13.6))

    plyr = layer.PressureLayer(snd_data["pres"][0], 50000.0)
    shr = winds.wind_shear(
        plyr, snd_data["pres"], snd_data["uwin"], snd_data["vwin"])

    assert (shr.u == pytest.approx(19.810976))
    assert (shr.v == pytest.approx(8.3890247))
