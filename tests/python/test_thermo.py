

import os
import pytest
import numpy as np
import pandas as pd

from nwsspc.sharp.calc import layer
from nwsspc.sharp.calc import thermo
from nwsspc.sharp.calc import parcel
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


def test_theta_wetbulb():
    wobf = parcel.lifter_wobus()
    cm1 = parcel.lifter_cm1()
    cm1.ma_type = thermo.adiabat.pseudo_liq

    thetaw = thermo.theta_wetbulb(wobf, 99300.0, 300.0, 297.0)
    assert (thetaw == pytest.approx(298.10607))

    thetaw = thermo.theta_wetbulb(
        wobf,
        snd_data["pres"][:3000],
        snd_data["tmpk"][:3000],
        snd_data["dwpk"][:3000]
    )
    assert (thetaw.min() == pytest.approx(287.99442))
    assert (thetaw.max() == pytest.approx(300.6817))
    assert (thetaw.mean() == pytest.approx(291.79718))

    thetaw = thermo.theta_wetbulb(cm1, 99300.0, 300.0, 297.0)
    assert (thetaw == pytest.approx(298.0095825))

    # TO-DO: The solver fails to converge at very low
    # pressures and temperatures.
    # Need to look into how to remedy this.
    thetaw = thermo.theta_wetbulb(
        cm1,
        snd_data["pres"][:3000],
        snd_data["tmpk"][:3000],
        snd_data["dwpk"][:3000]
    )

    assert (thetaw.min() == pytest.approx(288.2039))
    assert (thetaw.max() == pytest.approx(300.84494))
    assert (thetaw.mean() == pytest.approx(292.04013))

    # thetaw = thermo.theta_wetbulb(
    #     cm1, snd_data["pres"], snd_data["tmpk"], snd_data["dwpk"])
    # print(thetaw.min(), thetaw.max(), thetaw.mean())


def test_thetae():
    thetae = thermo.thetae(99300.00, 300.0, 295.0)
    assert (thetae == pytest.approx(351.2240295))

    thetae = thermo.thetae(
        snd_data["pres"], snd_data["tmpk"], snd_data["dwpk"])
    assert (thetae.min() == pytest.approx(319.05536))
    assert (thetae.max() == pytest.approx(510.92166))
    assert (thetae.mean() == pytest.approx(362.22058))


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
