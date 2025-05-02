

import os
import pytest
import numpy as np
import pandas as pd

from nwsspc.sharp.calc import layer
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


def test_lapse_rate():
    hlyr = layer.HeightLayer(0, 3000.0)
    plyr = layer.PressureLayer(70000.0, 50000.0)

    h_lr = thermo.lapse_rate(hlyr, snd_data["hght"], snd_data["tmpk"])
    p_lr = thermo.lapse_rate(
        plyr, snd_data["pres"], snd_data["hght"], snd_data["tmpk"])

    assert (h_lr == pytest.approx(6.63333))
    assert (p_lr == pytest.approx(8.45434))


def test_max_lapse_rate():
    search_hlyr = layer.HeightLayer(2000.0, 6000.0)
    search_plyr = layer.PressureLayer(80000.0, 50000.0)
    h_depth = 2000.0  # meters
    p_depth = 5000.0  # Pa

    h_max_lr, h_max_lyr = thermo.lapse_rate_max(
        search_hlyr,
        h_depth,
        snd_data["hght"],
        snd_data["tmpk"]
    )

    p_max_lr, p_max_lyr = thermo.lapse_rate_max(
        search_plyr,
        p_depth,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"]
    )

    assert (h_max_lr == pytest.approx(9.27))
    assert (h_max_lyr.bottom == 2000.0)
    assert (h_max_lyr.top == 4000.0)

    assert (p_max_lr == pytest.approx(9.4798536))
    assert (p_max_lyr.bottom == 68000.0)
    assert (p_max_lyr.top == 63000.0)
