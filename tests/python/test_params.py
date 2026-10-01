
import functools
import os
import pytest
import numpy as np
import pandas as pd

from nwsspc.sharp.calc import interp
from nwsspc.sharp.calc import layer
from nwsspc.sharp.calc import parcel
from nwsspc.sharp.calc import params
from nwsspc.sharp.calc import thermo
from nwsspc.sharp.calc import winds
from nwsspc.sharp.calc import constants


def load_parquet(filename):
    snd_df = pd.read_parquet(filename)
    snd_df = snd_df[snd_df["vwin"].notna()]
    snd_df = snd_df[snd_df["tmpc"].notna()]
    snd_df = snd_df[snd_df["relh"].notna()]
    snd_df = snd_df[snd_df["pres"] >= 50.0]

    pres = snd_df["pres"].to_numpy().astype('float32')*np.float32(100.0)
    hght = snd_df["hght"].to_numpy().astype('float32')
    tmpk = snd_df["tmpc"].to_numpy().astype('float32')+np.float32(273.15)
    dwpk = snd_df["dwpc"].to_numpy().astype('float32')+np.float32(273.15)
    relh = snd_df["relh"].to_numpy().astype('float32')
    wdir = snd_df["wdir"].to_numpy().astype('float32')
    wspd = snd_df["wspd"].to_numpy().astype('float32')
    uwin = snd_df["uwin"].to_numpy().astype('float32')
    vwin = snd_df["vwin"].to_numpy().astype('float32')

    hght_msl = hght.copy()
    # turn into height above ground level
    hght -= hght[0]

    # TO-DO - need a better interface to the API for doing this
    # uwin = np.empty(wspd.shape, dtype="float32")
    # vwin = np.empty(wspd.shape, dtype="float32")
    mixr = thermo.mixratio(pres, dwpk)
    vtmp = thermo.virtual_temperature(tmpk, mixr)
    theta = thermo.theta(pres, tmpk)
    thetae = thermo.thetae(
        pres,
        tmpk,
        dwpk
    )

    return {
        "pres": pres, "hght": hght, "hght_msl": hght_msl,
        "tmpk": tmpk, "mixr": mixr,
        "relh": relh,
        "theta": theta,
        "thetae": thetae,
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


def test_effective_inflow_layer_wobus():

    lifter = parcel.lifter_wobus()
    mupcl = parcel.Parcel()
    eil = params.effective_inflow_layer(
        lifter,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["dwpk"],
        snd_data["vtmp"],
        mupcl=mupcl
    )

    assert (eil.bottom == pytest.approx(92043.0))
    assert (eil.top == pytest.approx(83384.0))
    assert (mupcl.cape == pytest.approx(3353.4, abs=5e-1))
    assert (mupcl.cinh == pytest.approx(-34.5707, abs=5e-4))


def test_effective_inflow_layer_cm1():
    mupcl = parcel.Parcel()
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    eil = params.effective_inflow_layer(
        lifter,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["dwpk"],
        snd_data["vtmp"],
        mupcl=mupcl
    )

    assert (eil.bottom == pytest.approx(92043.0))
    assert (eil.top == pytest.approx(83432.0))
    assert (mupcl.cape == pytest.approx(3107.6428, abs=5e-1))
    assert (mupcl.cinh == pytest.approx(-36.41, abs=5e-1))


def test_bunkers_motion_nonparcel():
    # Non parcel based Bunkers motion
    shr_lyr = layer.HeightLayer(0, 6000.0)
    mw_lyr = layer.HeightLayer(0, 6000)
    storm_mtn = params.storm_motion_bunkers(
        snd_data["pres"],
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"],
        mw_lyr, shr_lyr,
    )

    assert (storm_mtn.u == pytest.approx(10.197119))
    assert (storm_mtn.v == pytest.approx(5.7385821))


def test_bunkers_motion():
    # parcel based Bunkers motion
    mupcl = parcel.Parcel()
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    eil = params.effective_inflow_layer(
        lifter,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["dwpk"],
        snd_data["vtmp"],
        mupcl=mupcl
    )

    storm_mtn = params.storm_motion_bunkers(
        snd_data["pres"],
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"],
        eil, mupcl
    )

    # The parcel lifter's exp/log/pow differ in the last bits across
    # platforms; 1e-3 allows 0.5 m of drift in the MU EL.
    assert (storm_mtn.u == pytest.approx(9.701575, abs=1e-3))
    assert (storm_mtn.v == pytest.approx(5.622299, abs=1e-3))


def test_corfidi_vectors():
    upshear, downshear = params.mcs_motion_corfidi(
        snd_data["pres"],
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"]
    )

    assert (upshear.u == pytest.approx(12.5269, abs=1e-3))
    assert (upshear.v == pytest.approx(2.7875, abs=1e-3))
    assert (downshear.u == pytest.approx(23.1313, abs=1e-3))
    assert (downshear.v == pytest.approx(15.99528, abs=1e-3))


def interp_log10p(p, pres, arr):
    return np.interp(np.log10(p), np.log10(pres.astype("float64"))[::-1],
                     arr.astype("float64")[::-1])


def ebwd_from_definition(pres, hght, uwin, vwin, eil, eql_pres):
    z = hght.astype("float64")

    def z_agl(p):
        return interp_log10p(p, pres, z) - z[0]

    def wind(agl):
        return (np.interp(z[0] + agl, z, uwin), np.interp(z[0] + agl, z, vwin))

    bot = z_agl(eil.bottom)
    top = bot + 0.5 * (z_agl(eql_pres) - bot)
    (u_bot, v_bot), (u_top, v_top) = wind(bot), wind(top)
    return bot, top, (u_top - u_bot, v_top - v_bot)


@functools.cache
def ddc_at_station(station):
    if station == "file":
        hght = snd_data["hght_msl"]
    else:
        hght = snd_data["hght"] + np.float32(station)
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    mupcl = parcel.Parcel()
    eil = params.effective_inflow_layer(
        lifter, snd_data["pres"], hght, snd_data["tmpk"], snd_data["dwpk"],
        snd_data["vtmp"], mupcl=mupcl)
    return hght, eil, mupcl


@pytest.mark.parametrize("station", [0.0, "file"])
def test_effective_bulk_wind(station):
    hght, eil, mupcl = ddc_at_station(station)
    ebwd_cmp = params.effective_bulk_wind_difference(
        snd_data["pres"],
        hght,
        snd_data["uwin"],
        snd_data["vwin"],
        eil,
        mupcl.eql_pressure
    )

    bot, top, expected = ebwd_from_definition(
        snd_data["pres"], hght, snd_data["uwin"],
        snd_data["vwin"], eil, mupcl.eql_pressure)
    assert (bot == 0.0)
    # The parcel lifter's exp/log/pow differ in the last bits across
    # platforms; these pins allow 0.5 m of drift in the MU EL.
    assert (top == pytest.approx(5866.089, abs=0.5))
    assert (expected[0] == pytest.approx(14.6, abs=1e-2))
    assert (expected[1] == pytest.approx(13.32178, abs=1e-2))

    ebwd = winds.vector_magnitude(ebwd_cmp.u, ebwd_cmp.v)
    assert (ebwd_cmp.u == pytest.approx(expected[0], abs=1e-4))
    assert (ebwd_cmp.v == pytest.approx(expected[1], abs=1e-4))
    assert (ebwd == pytest.approx(np.hypot(*expected), abs=1e-4))

    shear = winds.wind_shear(layer.HeightLayer(bot, top), hght,
                             snd_data["uwin"], snd_data["vwin"])
    assert (ebwd_cmp.u == pytest.approx(shear.u, abs=1e-4))
    assert (ebwd_cmp.v == pytest.approx(shear.v, abs=1e-4))


@pytest.mark.parametrize("shift", [0.0, 1000.0, 1234.5, 762.3])
def test_effective_bulk_wind_station_height(shift):
    pres = np.array([100000, 90000, 80000, 70000, 60000, 50000],
                    dtype="float32")
    hght = np.array([0, 1000, 2000, 3000, 4000, 5000],
                    dtype="float32") + np.float32(shift)
    uwin = np.array([0, 10, 12, 13, 13, 13], dtype="float32")
    vwin = np.zeros(6, dtype="float32")
    ebwd = params.effective_bulk_wind_difference(
        pres, hght, uwin, vwin, layer.PressureLayer(100000, 90000), 60000)
    assert (ebwd.u == pytest.approx(12.0))
    assert (ebwd.v == 0.0)


def assert_missing_wind(wind):
    assert (wind.u == constants.MISSING and wind.v == constants.MISSING)


def test_wind_params_missing_layer():
    pres = np.array([100000, 95000, 90000, 85000, 80000], dtype="float32")
    hght = np.array([0, 500, 1000, 1500, 2000], dtype="float32")
    uwin = np.array([0, 5, 10, 15, 20], dtype="float32")
    vwin = np.array([0, 2, 4, 6, 8], dtype="float32")

    assert_missing_wind(params.effective_bulk_wind_difference(
        pres, hght, uwin, vwin, layer.PressureLayer(100000, 95000), 70000))

    assert_missing_wind(params.storm_motion_bunkers(
        pres, hght, uwin, vwin, layer.HeightLayer(0, 3000),
        layer.HeightLayer(0, 2000)))

    M = constants.MISSING
    assert_missing_wind(params.storm_motion_bunkers(
        pres, hght, uwin, vwin, layer.HeightLayer(0, 2000),
        layer.HeightLayer(M, M)))

    for vector in params.mcs_motion_corfidi(pres, hght / 2, uwin, vwin):
        assert_missing_wind(vector)

    pres6 = np.array([100000, 85000, 70000, 59000, 51000, 40000],
                     dtype="float32")
    hght6 = np.array([0, 1500, 3000, 4500, 5500, 7000], dtype="float32")
    uwin6 = np.array([0, 10, 20, 30, 35, 40], dtype="float32")
    vwin6 = np.zeros(6, dtype="float32")
    mu_pcl = parcel.Parcel()
    mu_pcl.cape = 3000.0
    mu_pcl.eql_pressure = 30000.0
    storm = winds.WindComponents()
    storm.u, storm.v = 5.0, 5.0
    hgz = layer.PressureLayer(59000.0, 51000.0)
    assert (params.large_hail_parameter(mu_pcl, 8.0, hgz, storm, pres6, hght6,
                                        uwin6, vwin6) == M)
    mu_pcl.eql_pressure = 51000.0
    assert (params.large_hail_parameter(mu_pcl, 8.0, hgz, storm, pres6, hght6,
                                        uwin6, vwin6) ==
            pytest.approx(126.17279, abs=1e-3))


def test_bunkers_motion_effective_fallback():
    pres = np.array([100000, 80000, 62000, 47000, 35000], dtype="float32")
    hght = np.array([0, 2000, 4000, 6000, 8000], dtype="float32")
    uwin = np.array([0, 10, 20, 30, 40], dtype="float32")
    vwin = np.zeros(5, dtype="float32")
    mupcl = parcel.Parcel()
    mupcl.eql_pressure = 80000.0
    motion = params.storm_motion_bunkers(
        pres, hght, uwin, vwin, layer.PressureLayer(105000, 95000), mupcl)
    fallback = params.storm_motion_bunkers(
        pres, hght, uwin, vwin, layer.HeightLayer(0, 6000),
        layer.HeightLayer(0, 6000))
    assert (motion.u == fallback.u and motion.v == fallback.v)
    assert (motion.u == pytest.approx(14.0566034))
    assert (motion.v == pytest.approx(-7.5))


def bunkers_from_definition(pres, hght, uwin, vwin, eil, eql_pres, left):
    p = pres.astype("float64")
    z_msl = hght.astype("float64")
    sfc = z_msl[0]

    def pres_at_agl(z_agl):
        return np.interp(sfc + z_agl, z_msl, p)

    def mean_wind(bot_agl, top_agl, weighted):
        pbot, ptop = pres_at_agl(bot_agl), pres_at_agl(top_agl)
        inside = (p < pbot) & (p > ptop)
        pp = np.concatenate(([pbot], p[inside], [ptop]))
        w = pp if weighted else np.ones_like(pp)
        mean = []
        for arr in (uwin, vwin):
            a = np.concatenate(([interp_log10p(pbot, pres, arr)], arr[inside],
                                [interp_log10p(ptop, pres, arr)]))
            mean.append(np.trapezoid(a * w, pp) / np.trapezoid(w, pp))
        return mean

    base = interp_log10p(eil.bottom, pres, z_msl) - sfc
    top = 0.65 * (interp_log10p(eql_pres, pres, z_msl) - sfc)
    if top - base < 3000.0:
        mw_layer, mean = None, mean_wind(0.0, 6000.0, False)
    else:
        mw_layer, mean = (base, top), mean_wind(base, top, True)
    lo = mean_wind(0.0, 500.0, False)
    hi = mean_wind(5500.0, 6000.0, False)
    shr_u, shr_v = hi[0] - lo[0], hi[1] - lo[1]
    k = (-7.5 if left else 7.5) / np.hypot(shr_u, shr_v)
    return mw_layer, (mean[0] + k * shr_v, mean[1] - k * shr_u)


def check_effective_bunkers(pres, hght, uwin, vwin, eil, mupcl, mw_layer):
    for left in (True, False):
        motion = params.storm_motion_bunkers(
            pres, hght, uwin, vwin, eil, mupcl, left)
        classic = params.storm_motion_bunkers(
            pres, hght, uwin, vwin,
            layer.HeightLayer(*(mw_layer or (0.0, 6000.0))),
            layer.HeightLayer(0.0, 6000.0), left, mw_layer is not None)
        assert (motion.u == pytest.approx(classic.u))
        assert (motion.v == pytest.approx(classic.v))

        oracle_layer, oracle = bunkers_from_definition(
            pres, hght, uwin, vwin, eil, mupcl.eql_pressure, left)
        if mw_layer is None:
            assert (oracle_layer is None)
        else:
            assert (oracle_layer == pytest.approx(mw_layer, abs=1e-2))
        assert (motion.u == pytest.approx(oracle[0], abs=1e-4))
        assert (motion.v == pytest.approx(oracle[1], abs=1e-4))
    return motion


@pytest.mark.parametrize("station", [0.0, "file", 1000.0, 762.3])
@pytest.mark.parametrize("eil_bottom, eil_top, expected", [
    (92043.0, 83432.0, (9.701575, 5.622300)),
    (85000.0, 75000.0, (12.387376, 6.366230)),
    (80000.0, 70000.0, (13.978964, 6.661992)),
])
def test_bunkers_motion_effective_layer(station, eil_bottom, eil_top,
                                        expected):
    pres = snd_data["pres"]
    hght, _, mupcl = ddc_at_station(station)
    assert (pres[0] == 92043.0)

    eil = layer.PressureLayer(eil_bottom, eil_top)
    base = layer.pressure_layer_to_height(eil, pres, hght, True).bottom
    el = interp.interp_pressure(mupcl.eql_pressure, pres, hght) - hght[0]
    # The parcel lifter's exp/log/pow differ in the last bits across
    # platforms; these pins allow 0.5 m of drift in the MU EL.
    assert (el == pytest.approx(11732.17, abs=0.5))

    motion = check_effective_bunkers(pres, hght, snd_data["uwin"],
                                     snd_data["vwin"], eil, mupcl,
                                     (base, 0.65 * el))
    assert (motion.u == pytest.approx(expected[0], abs=1e-3))
    assert (motion.v == pytest.approx(expected[1], abs=1e-3))


@pytest.mark.parametrize("station", [0.0, 1000.0, 762.3])
@pytest.mark.parametrize("base, el, mw_layer, expected", [
    (2000, 7000, None, (15.849664, -0.509418)),
    (6000, 14000, (6000.0, 9100.0), (27.916351, -0.846436)),
])
def test_bunkers_motion_effective_minimum_depth(station, base, el, mw_layer,
                                                expected):
    z = np.arange(33, dtype="float32") * np.float32(500.0)
    pres = np.float32(100000.0) * np.exp(-z / np.float32(8000.0))
    hght = z + np.float32(station)
    uwin = np.float32(30.0) * (np.float32(1.0) -
                               np.exp(-z / np.float32(4000.0)))
    vwin = np.float32(10.0) * np.sin(z / np.float32(3000.0))

    mupcl = parcel.Parcel()
    mupcl.eql_pressure = float(pres[el // 500])
    eil = layer.PressureLayer(float(pres[base // 500]),
                              float(pres[base // 500 + 2]))
    motion = check_effective_bunkers(pres, hght, uwin, vwin, eil, mupcl,
                                     mw_layer)
    # numpy's float32 exp and sin differ by an ULP across platforms.
    assert (motion.u == pytest.approx(expected[0], abs=1e-4))
    assert (motion.v == pytest.approx(expected[1], abs=1e-4))

def test_stp_scp_ship_dcp_lhp():
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    # get the mixed-layer parcel
    mix_lyr = layer.PressureLayer(
        snd_data["pres"][0], snd_data["pres"][0] - 10000.0)
    theta = thermo.theta(snd_data["pres"], snd_data["tmpk"])
    pcl = parcel.Parcel.mixed_layer_parcel(
        mix_lyr,
        snd_data["pres"],
        theta,
        snd_data["mixr"]
    )

    search_layer = layer.PressureLayer(
        snd_data["pres"][0],
        snd_data["pres"][0] - 40000.0
    )

    dcape_pcl = parcel.DowndraftParcel.min_thetae(
        search_layer,
        snd_data["pres"],
        snd_data["tmpk"],
        snd_data["dwpk"],
        snd_data["thetae"]
    )

    dpcl_t = dcape_pcl.lower_parcel(lifter, snd_data["pres"])
    dpcl_buoy = thermo.buoyancy(dpcl_t, snd_data["tmpk"])

    dcape, dcinh = dcape_pcl.cape_cinh(
        snd_data["pres"], snd_data["hght"], dpcl_buoy)

    # lift the parcel and get CAPE
    vtmpk = pcl.lift_parcel(lifter, snd_data["pres"])
    buoy = thermo.buoyancy(vtmpk, snd_data["vtmp"])
    cape, cinh = pcl.cape_cinh(snd_data["pres"], snd_data["hght"], buoy)

    # Get the effective inflow layer for effective SRH
    mupcl = parcel.Parcel()
    eil = params.effective_inflow_layer(
        lifter,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["dwpk"],
        snd_data["vtmp"],
        mupcl=mupcl
    )

    # Get the storm relative helicity for the effective inflow layer
    storm_mtn = params.storm_motion_bunkers(
        snd_data["pres"],
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"],
        eil, mupcl
    )
    esrh = winds.helicity(
        eil,
        storm_mtn,
        snd_data["pres"],
        snd_data["uwin"],
        snd_data["vwin"]
    )

    ebwd_cmp = params.effective_bulk_wind_difference(
        snd_data["pres"],
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"],
        eil,
        mupcl.eql_pressure
    )

    ebwd = winds.vector_magnitude(ebwd_cmp.u, ebwd_cmp.v)

    # Get the LCL height in meters AGL
    lcl_hght = interp.interp_pressure(
        pcl.lcl_pressure,
        snd_data["pres"],
        snd_data["hght"]
    ) - snd_data["hght"][0]

    stp = params.significant_tornado_parameter(
        pcl,
        lcl_hght,
        esrh,
        ebwd
    )
    assert (stp == pytest.approx(0.48378, abs=1e-4))

    scp = params.supercell_composite_parameter(mupcl.cape, esrh, ebwd)
    assert (scp == pytest.approx(7.9699, abs=1e-1))

    # get SHIP
    plyr = layer.PressureLayer(70000.0, 50000.0)
    hlyr = layer.HeightLayer(0.0, 6000.0)
    lr75 = thermo.lapse_rate(
        plyr, snd_data["pres"], snd_data["hght"], snd_data["tmpk"])
    t500 = interp.interp_pressure(50000.0, snd_data["pres"], snd_data["tmpk"])
    fzl = interp.find_first_height(
        constants.ZEROCNK, snd_data["hght"], snd_data["tmpk"])
    shr06 = winds.wind_shear(
        hlyr, snd_data["hght"], snd_data["uwin"], snd_data["vwin"])
    shr06 = winds.vector_magnitude(shr06.u, shr06.v)
    ship = params.significant_hail_parameter(mupcl, lr75, t500, fzl, shr06)
    assert (ship == pytest.approx(1.9521, abs=1e-3))

    hlyr_in_p = layer.height_layer_to_pressure(
        hlyr, snd_data["pres"], snd_data["hght"])
    mw06 = winds.mean_wind(
        hlyr_in_p, snd_data["pres"], snd_data["uwin"], snd_data["vwin"]
    )
    mw06 = winds.vector_magnitude(mw06.u, mw06.v)
    dcp = params.derecho_composite_parameter(dcape, mupcl.cape, shr06, mw06)
    assert (dcp == pytest.approx(7.63, abs=1e-1))

    hgz = params.hail_growth_layer(snd_data["pres"], snd_data["tmpk"])
    lhp = params.large_hail_parameter(
        mupcl,
        lr75,
        hgz,
        storm_mtn,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"]
    )
    assert (lhp == pytest.approx(14.1323, abs=1e-2))


def test_ehi():
    pres = snd_data["pres"][0]
    tmpk = snd_data["tmpk"][0]
    dwpk = snd_data["dwpk"][0]

    pcl = parcel.Parcel.surface_parcel(pres, tmpk, dwpk)

    # Wobus Lifter
    lifter = parcel.lifter_wobus()
    vtmpk = pcl.lift_parcel(lifter, snd_data["pres"])
    buoy = thermo.buoyancy(vtmpk, snd_data["vtmp"])
    cape, cinh = pcl.cape_cinh(snd_data["pres"], snd_data["hght"], buoy)

    srh_lyr = layer.HeightLayer(0.0, 3000.0)
    mupcl = parcel.Parcel()
    eil = params.effective_inflow_layer(
        lifter,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["dwpk"],
        snd_data["vtmp"],
        mupcl=mupcl
    )

    # Get the storm relative helicity for the effective inflow layer
    storm_mtn = params.storm_motion_bunkers(
        snd_data["pres"],
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"],
        eil, mupcl
    )
    srh = winds.helicity(
        srh_lyr,
        storm_mtn,
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"]
    )

    ehi = params.energy_helicity_index(pcl.cape, srh)
    # The parcel lifters' results differ in the last bits across
    # platforms; 1e-3 allows 0.5 m of drift in the EL.
    assert (ehi == pytest.approx(4.38889, abs=1e-3))


def test_convective_temperature():
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    cnvtv_tmpk = params.convective_temperature(
        lifter,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["vtmp"],
        snd_data["mixr"]
    )
    assert (cnvtv_tmpk == pytest.approx(304.65, abs=0.5))


def test_precipitable_water():
    plyr = layer.PressureLayer(snd_data["pres"][0], 40000.0)
    pwat = params.precipitable_water(plyr, snd_data["pres"], snd_data["mixr"])

    assert (pwat == pytest.approx(21.11469))


def test_hgz():
    hgz = params.hail_growth_layer(snd_data["pres"], snd_data["tmpk"])
    assert (hgz.bottom == pytest.approx(51430))
    assert (hgz.top == pytest.approx(35816))


def test_dgz():
    dgz = params.dendritic_layer(snd_data["pres"], snd_data["tmpk"])
    assert (dgz.bottom == pytest.approx(49598))
    assert (dgz.top == pytest.approx(45961))


def test_fwwi():
    fwwi = params.fosberg_fire_index(
        308,
        0.00001,
        13.5
    )
    assert (fwwi == 100.0)

    tmpk = np.array([308, 308, 308], dtype='float32')
    relh = np.array([0.00001, 0.00001, 0.00001], dtype='float32')
    wspd = np.array([13.5, 13.5, 13.5], dtype='float32')

    fwwi = params.fosberg_fire_index(tmpk, relh, wspd)
    assert (fwwi == np.array([100.0, 100.0, 100.0], dtype='float32')).all()


def test_pft():
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    mix_layer = layer.PressureLayer(
        snd_data["pres"][0], snd_data["pres"][0] - 10000.0)
    pft = params.pyrocumulonimbus_firepower_threshold(
        lifter,
        mix_layer,
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["mixr"],
        snd_data["vtmp"],
        snd_data["uwin"],
        snd_data["vwin"],
        snd_data["theta"]
    )
    assert (pft == pytest.approx(158187356160.0, abs=1e6))


def test_pft_missing():
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    M = constants.MISSING
    pres = snd_data["pres"]
    mix_layer = layer.PressureLayer(pres[0], pres[0] - 10000.0)
    theta = snd_data["theta"].copy()
    theta[pres < 75000.0] = M
    pcl = parcel.Parcel()
    pft = params.pyrocumulonimbus_firepower_threshold(
        lifter, mix_layer, pres, snd_data["hght"], snd_data["tmpk"],
        snd_data["mixr"], snd_data["vtmp"], snd_data["uwin"], snd_data["vwin"],
        theta, pcl=pcl)
    assert (pft == M)
    assert (pcl.pres == M)
