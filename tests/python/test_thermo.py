

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


def test_wobf():
    wobf = thermo.wobf(300.0)
    assert (wobf == pytest.approx(292.115967))

    wobf = thermo.wobf(snd_data["tmpk"])
    assert (wobf.min() == pytest.approx(273.2617))
    assert (wobf.max() == pytest.approx(291.03946))
    assert (wobf.mean() == pytest.approx(276.11484))


def test_lcl_temp():
    lcl_t = thermo.lcl_temperature(300.0, 297.0)
    assert (lcl_t == pytest.approx(296.27252))

    lcl_t = thermo.lcl_temperature(snd_data["tmpk"], snd_data["dwpk"])
    assert (lcl_t.min() == pytest.approx(175.60263))
    assert (lcl_t.max() == pytest.approx(289.98557))
    assert (lcl_t.mean() == pytest.approx(211.24747))


def test_vappres():
    vappres = thermo.vapor_pressure(100000.0, 300.0)
    assert (vappres == pytest.approx(3534.52))

    vappres = thermo.vapor_pressure(snd_data["pres"], snd_data["tmpk"])
    assert (vappres.min() == pytest.approx(0.9944834))
    assert (vappres.max() == pytest.approx(3167.4294))
    assert (vappres.mean() == pytest.approx(339.91318))


def test_vappres_ice():
    vappres = thermo.vapor_pressure_ice(40000.0, 250.0)
    assert (vappres == pytest.approx(75.6271286))

    vappres = thermo.vapor_pressure_ice(snd_data["pres"], snd_data["tmpk"])
    assert (vappres.min() == pytest.approx(0.5158234))
    assert (vappres.max() == pytest.approx(4015.757))
    assert (vappres.mean() == pytest.approx(380.8075))


def test_tmpk_at_mixr():
    dwpk = thermo.temperature_at_mixratio(snd_data["mixr"], snd_data["pres"])
    assert (dwpk == pytest.approx(snd_data["dwpk"], abs=1e-3))


def test_theta_lvl():
    plev = thermo.theta_level(309.1682, 300.0)
    assert (plev == pytest.approx(90000.0))


def test_theta():
    theta = thermo.theta(90000.0, 300.0)
    assert (theta == pytest.approx(309.1682))

    theta = thermo.theta(snd_data["pres"], snd_data["tmpk"])
    assert (theta.min() == pytest.approx(303.79684))
    assert (theta.max() == pytest.approx(510.80814))
    assert (theta.mean() == pytest.approx(358.42636))


def test_mixr():
    mixr = thermo.mixratio(100000.0, 300.0)
    assert (mixr == pytest.approx(0.022788666))

    mixr = thermo.mixratio(snd_data["pres"], snd_data["dwpk"])
    assert (mixr.max() == pytest.approx(0.014611241))
    assert (mixr.mean() == pytest.approx(0.0011161759))


def test_ice_mixr():
    mixr = thermo.mixratio_ice(50000.0, 250.0)
    assert (mixr == pytest.approx(0.0009421613067))

    mixr = thermo.mixratio_ice(snd_data["pres"], snd_data["tmpk"])
    assert (mixr.max() == pytest.approx(0.02837335))
    assert (mixr.mean() == pytest.approx(0.003338528))


def test_spechum():
    spfh = thermo.specific_humidity(0.02)
    assert (spfh == pytest.approx(0.01960784))

    spfh = thermo.specific_humidity(snd_data["mixr"])
    assert (spfh.max() == pytest.approx(0.014400828))
    assert (spfh.mean() == pytest.approx(0.0011059974))


def test_virtemp():
    vtmpk = thermo.virtual_temperature(300.0, 0.0)
    assert (vtmpk == 300.0)

    vtmpk = thermo.virtual_temperature(300.0, 0.015)
    assert (vtmpk == pytest.approx(302.69479))

    vtmpk = thermo.virtual_temperature(300.0, 0.015, 0.01)
    assert (vtmpk == pytest.approx(299.741699))

    vtmpk = thermo.virtual_temperature(300.0, 0.015, ri=0.01)
    assert (vtmpk == pytest.approx(299.741699))

    vtmpk = thermo.virtual_temperature(300.0, 0.015, rl=0.0015, ri=0.01)
    assert (vtmpk == pytest.approx(299.30368))

    vtmpk = thermo.virtual_temperature(snd_data["tmpk"], snd_data["mixr"])
    assert (vtmpk.min() == pytest.approx(208.25008))
    assert (vtmpk.max() == pytest.approx(300.75977))
    assert (vtmpk.mean() == pytest.approx(238.928))


def test_drylift():
    lcl_p, lcl_t = thermo.drylift(100000.0, 300.0, 297.0)
    assert (lcl_p == pytest.approx(95718.40625))
    assert (lcl_t == pytest.approx(296.27252))


def test_wetbulb():
    wobf = parcel.lifter_wobus()
    cm1 = parcel.lifter_cm1()
    cm1.ma_type = thermo.adiabat.pseudo_liq

    wetbulb = thermo.wetbulb(wobf, 99300.0, 300.0, 297.0)
    assert (wetbulb == pytest.approx(297.85162))

    wetbulb = thermo.wetbulb(
        wobf,
        snd_data["pres"][:3000],
        snd_data["tmpk"][:3000],
        snd_data["dwpk"][:3000]
    )
    assert (wetbulb.min() == pytest.approx(209.90778))
    assert (wetbulb.max() == pytest.approx(293.51944))
    assert (wetbulb.mean() == pytest.approx(245.64085))

    wetbulb = thermo.wetbulb(cm1, 99300.0, 300.0, 297.0)
    assert (wetbulb == pytest.approx(297.76919555))

    # TO-DO: The solver fails to converge at very low
    # pressures and temperatures.
    # Need to look into how to remedy this.
    wetbulb = thermo.wetbulb(
        cm1,
        snd_data["pres"][:3000],
        snd_data["tmpk"][:3000],
        snd_data["dwpk"][:3000]
    )

    assert (wetbulb.min() == pytest.approx(210.11125))
    assert (wetbulb.max() == pytest.approx(293.46805))
    assert (wetbulb.mean() == pytest.approx(245.90587))


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
