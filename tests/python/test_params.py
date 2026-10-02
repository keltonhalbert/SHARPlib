
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
    uwin = snd_data["uwin"].copy()
    uwin[pres > 80000.0] = M
    for u, th in ((snd_data["uwin"], theta), (uwin, snd_data["theta"])):
        pcl = parcel.Parcel()
        pft = params.pyrocumulonimbus_firepower_threshold(
            lifter, mix_layer, pres, snd_data["hght"], snd_data["tmpk"],
            snd_data["mixr"], snd_data["vtmp"], u, snd_data["vwin"], th,
            pcl=pcl)
        assert (pft == M)
        assert (pcl.pres == M)


# ===========================================================================
# Precipitation type: the modified Bourgouin method (Birk et al. 2021)
# ===========================================================================

def test_bourgouin_energy_struct():
    energy = params.BourgouinEnergy()
    assert energy.melting_energy_total == constants.MISSING
    assert energy.melting_energy_aloft == constants.MISSING
    assert energy.refreezing_energy == constants.MISSING

    energy = params.BourgouinEnergy(
        melting_energy_total=4.0,
        melting_energy_aloft=2.0,
        refreezing_energy=180.0
    )
    assert energy.melting_energy_total == 4.0
    assert energy.melting_energy_aloft == 2.0
    assert energy.refreezing_energy == 180.0

    energy.refreezing_energy = 100.0
    assert energy.refreezing_energy == 100.0


def test_precip_type_probabilities_struct():
    probs = params.PrecipTypeProbabilities()
    assert probs.rain == constants.MISSING
    assert probs.snow == constants.MISSING
    assert probs.freezing_rain == constants.MISSING
    assert probs.ice_pellets == constants.MISSING

    probs.rain = 0.5
    assert probs.rain == 0.5


# ---------------------------------------------------------------------------
# Wet-bulb melting and refreezing energies from a sounding
# ---------------------------------------------------------------------------

def _melting_area(hght, dtw):
    """
    Closed-form area of max(dtw, 0) under a piecewise-linear profile (K m):
    a trapezoid for each segment that doesn't cross 0, and the warm
    triangle of each segment that does.
    """
    dz = np.diff(hght)
    d0, d1 = dtw[:-1], dtw[1:]
    w0, w1 = np.maximum(d0, 0.0), np.maximum(d1, 0.0)
    crossing = d0 * d1 < 0.0
    span = np.where(crossing, np.abs(d0) + np.abs(d1), 1.0)
    triangle = 0.5 * dz * (w0**2 + w1**2) / span
    trapezoid = 0.5 * (w0 + w1) * dz
    return np.sum(np.where(crossing, triangle, trapezoid))


def test_bourgouin_pressure_min_constant():
    # The named default of pressure_min (Pa)
    assert params.BOURGOUIN_PRESSURE_MIN == 25000.0

    # The top level is above 250 hPa and the top segment is warm, so the
    # default cap drops all the melting, the same as passing the constant.
    pres = np.array([100000, 50000, 30000, 20000], dtype="float32")
    hght = np.array([0, 5500, 9000, 11500], dtype="float32")
    tw = (constants.ZEROCNK +
          np.array([-2.0, -1.0, -1.0, 5.0])).astype("float32")
    default = params.bourgouin_energy(pres, hght, tw)
    named = params.bourgouin_energy(
        pres, hght, tw, pressure_min=params.BOURGOUIN_PRESSURE_MIN)
    uncapped = params.bourgouin_energy(pres, hght, tw, pressure_min=0.0)
    assert default.melting_energy_total == 0.0
    assert named.melting_energy_total == default.melting_energy_total
    assert uncapped.melting_energy_total > 0.0


def test_bourgouin_energy_cap_profile():
    # Every level is at or above 250 hPa, and the top segment is warm, so
    # losing the top level would show up in the melting energy.
    pres = np.array([100000, 94500, 87000, 79500, 70000, 57500, 43000,
                     32000, 26500], dtype="float32")
    hght = np.array([0, 500, 1200, 2000, 3000, 4500, 6500, 8500, 10000],
                    dtype="float32")
    zerocnk = np.float32(constants.ZEROCNK)
    tmpk = zerocnk + np.array([-1.5, -0.5, 1.0, 2.0, -1.0, 0.5, -2.0, 0.25,
                               3.0], dtype="float32")
    # The departures the library sees, in double precision
    dtw = tmpk.astype("float64") - np.float64(zerocnk)
    z = hght.astype("float64")
    g_t0 = constants.GRAVITY / constants.ZEROCNK

    # The surface is cold, so the melting aloft is all of it. The
    # near-surface cold layer ends inside the second segment.
    melting = g_t0 * _melting_area(z, dtw)
    refreezing = g_t0 * _melting_area(z[:3], -dtw[:3])

    energy = params.bourgouin_energy(pres, hght, tmpk)
    assert energy.melting_energy_total == pytest.approx(melting, rel=1e-5)
    assert energy.melting_energy_aloft == pytest.approx(melting, rel=1e-5)
    assert energy.refreezing_energy == pytest.approx(refreezing, rel=1e-5)

    # The defaults are min_energy = 0 and pressure_min = 250 hPa, and with
    # every level at or above 250 hPa, removing the cap changes nothing.
    for other in (
        params.bourgouin_energy(pres, hght, tmpk, min_energy=0.0,
                                pressure_min=25000.0),
        params.bourgouin_energy(pres, hght, tmpk, pressure_min=0.0),
    ):
        assert other.melting_energy_total == energy.melting_energy_total
        assert other.melting_energy_aloft == energy.melting_energy_aloft
        assert other.refreezing_energy == energy.refreezing_energy

    # A cap at 400 hPa drops the top two levels. The wet-bulb above the cap
    # is never read, so it can be MISSING.
    melting = g_t0 * _melting_area(z[:7], dtw[:7])
    tmpk[7:] = constants.MISSING
    energy = params.bourgouin_energy(pres, hght, tmpk, pressure_min=40000.0)
    assert energy.melting_energy_total == pytest.approx(melting, rel=1e-5)
    assert energy.melting_energy_aloft == pytest.approx(melting, rel=1e-5)
    assert energy.refreezing_energy == pytest.approx(refreezing, rel=1e-5)


def test_bourgouin_energy_min_energy():
    # From the surface up: cold 100, warm 1.99, cold 80, and warm 20 J/kg,
    # each a triangle peaking 4 K from 0 C.
    energies = [-100.0, 1.99, -80.0, 20.0]
    peak = 4.0
    g_t0 = constants.GRAVITY / constants.ZEROCNK
    hght, dtw = [0.0], [0.0]
    for e in energies:
        depth = 2.0 * abs(e) / (g_t0 * peak)
        hght += [hght[-1] + depth / 2.0, hght[-1] + depth]
        dtw += [np.copysign(peak, e), 0.0]
    hght = np.array(hght, dtype="float32")
    pres = np.float32(100000.0) - np.float32(7.0) * hght
    tmpk = np.float32(constants.ZEROCNK) + np.array(dtw, dtype="float32")

    # With min_energy = 2 J/kg, the weak warm layer merges the two cold
    # layers into one near-surface cold layer.
    energy = params.bourgouin_energy(pres, hght, tmpk, min_energy=2.0)
    assert energy.melting_energy_total == pytest.approx(21.99, rel=1e-5)
    assert energy.melting_energy_aloft == pytest.approx(20.0, rel=1e-5)
    assert energy.refreezing_energy == pytest.approx(180.0, rel=1e-5)

    energy = params.bourgouin_energy(pres, hght, tmpk)
    assert energy.melting_energy_total == pytest.approx(21.99, rel=1e-5)
    assert energy.melting_energy_aloft == pytest.approx(21.99, rel=1e-5)
    assert energy.refreezing_energy == pytest.approx(100.0, rel=1e-5)


def test_bourgouin_energy_missing():
    empty = np.array([], dtype="float32")
    energy = params.bourgouin_energy(empty, empty, empty)
    assert energy.melting_energy_total == constants.MISSING
    assert energy.melting_energy_aloft == constants.MISSING
    assert energy.refreezing_energy == constants.MISSING

    # One level below the cap is not enough
    pres = np.array([30000, 20000], dtype="float32")
    hght = np.array([0, 1000], dtype="float32")
    tmpk = np.array([280, 280], dtype="float32")
    energy = params.bourgouin_energy(pres, hght, tmpk)
    assert energy.melting_energy_total == constants.MISSING
    assert energy.melting_energy_aloft == constants.MISSING
    assert energy.refreezing_energy == constants.MISSING

    with pytest.raises(BufferError):
        params.bourgouin_energy(pres, hght, tmpk[:1])


# ---------------------------------------------------------------------------
# Precipitation generation layer from a sounding
# ---------------------------------------------------------------------------

def relh_sounding(hght, relh, tmpk=283.15):
    """
    Pressure, height, temperature, and dewpoint arrays with the given
    relative humidities over liquid, at a temperature above 0 C.
    """
    hght = np.asarray(hght, dtype='float32')
    pres = np.float32(100000.0) - np.float32(10.0) * hght
    # Invert the saturation vapor pressure, 611.2 exp(17.67 Tc / (Tc + 243.5))
    tmpc = tmpk - 273.15
    log_vapr = np.log(relh) + 17.67 * tmpc / (tmpc + 243.5)
    dwpk = 273.15 + 243.5 * log_vapr / (17.67 - log_vapr)
    return (
        pres,
        hght,
        np.full(hght.shape, tmpk, dtype='float32'),
        np.asarray(dwpk, dtype='float32'),
    )


def test_precipitation_generation_layer_phase():
    pres = np.array([100000.0, 90000.0, 80000.0], dtype='float32')
    hght = np.array([0.0, 1000.0, 2000.0], dtype='float32')

    # Above 0 C, relative humidity is over liquid: moist, though dry over ice
    tmpk = np.full(3, 283.15, dtype='float32')
    dwpk = np.full(3, 280.0, dtype='float32')
    assert thermo.relative_humidity(
        pres[0], tmpk[0], dwpk[0]) == pytest.approx(0.808, abs=5e-4)
    assert thermo.relative_humidity_ice(
        pres[0], tmpk[0], dwpk[0]) == pytest.approx(0.733, abs=5e-4)
    lyr = params.precipitation_generation_layer(pres, hght, tmpk, dwpk)
    assert (lyr.bottom, lyr.top) == (0.0, 2000.0)

    # Below 0 C, relative humidity is over ice: moist, though dry over liquid
    tmpk = np.full(3, 263.15, dtype='float32')
    dwpk = np.full(3, 259.15, dtype='float32')
    assert thermo.relative_humidity(pres[0], tmpk[0], dwpk[0]) < 0.75
    assert thermo.relative_humidity_ice(pres[0], tmpk[0], dwpk[0]) > 0.75
    lyr = params.precipitation_generation_layer(pres, hght, tmpk, dwpk)
    assert (lyr.bottom, lyr.top) == (0.0, 2000.0)


def test_precipitation_generation_layer_elimination():
    # Dry below moist, crossing 75 % halfway between 0.5 and 1.0. A 1400 m
    # dry layer at the surface leaves the 1300 m moist layer above it.
    snd = relh_sounding([0.0, 1300.0, 1500.0, 2700.0], [0.5, 0.5, 1.0, 1.0])
    lyr = params.precipitation_generation_layer(*snd)
    assert lyr.bottom == pytest.approx(1400.0, abs=1e-2)
    assert lyr.top == 2700.0

    # A 1600 m dry layer eliminates it
    snd = relh_sounding([0.0, 1500.0, 1700.0, 2900.0], [0.5, 0.5, 1.0, 1.0])
    lyr = params.precipitation_generation_layer(*snd)
    assert lyr.bottom == constants.MISSING
    assert lyr.top == constants.MISSING


def test_precipitation_generation_layer_min_depth():
    # A 100 m dry sliver splits a moist layer into two 600 m layers
    snd = relh_sounding(
        [0.0, 550.0, 650.0, 750.0, 1300.0], [1.0, 1.0, 0.5, 1.0, 1.0])
    for lyr in (
        params.precipitation_generation_layer(*snd),
        params.precipitation_generation_layer(*snd, min_depth=0.0),
    ):
        assert lyr.bottom == constants.MISSING
        assert lyr.top == constants.MISSING

    # With min_depth = 150 m, the sliver is absorbed
    lyr = params.precipitation_generation_layer(*snd, min_depth=150.0)
    assert (lyr.bottom, lyr.top) == (0.0, 1300.0)


def test_precipitation_generation_layer_empty():
    empty = np.array([], dtype='float32')
    lyr = params.precipitation_generation_layer(empty, empty, empty, empty)
    assert lyr.bottom == constants.MISSING
    assert lyr.top == constants.MISSING


# ---------------------------------------------------------------------------
# Probability of ice, and precipitation-type probabilities from energies
# ---------------------------------------------------------------------------

def test_probability_of_ice():
    zerocnk = np.float32(constants.ZEROCNK)
    assert params.probability_of_ice(zerocnk - np.float32(15.0)) == 1.0
    assert params.probability_of_ice(zerocnk - np.float32(7.0)) == 0.0
    # float32 evaluation: the polynomial's terms cancel, so allow 1e-5
    assert params.probability_of_ice(
        zerocnk - np.float32(11.0)) == pytest.approx(0.583474, rel=1e-5)
    assert params.probability_of_ice(
        260.5) == pytest.approx(0.7286607, rel=1e-5)

    assert params.probability_of_ice(constants.MISSING) == constants.MISSING
    assert params.probability_of_ice(np.nan) == constants.MISSING


def test_modified_bourgouin_from_energies():
    cold_sfc = constants.ZEROCNK - 2.0
    warm_sfc = constants.ZEROCNK + 2.0

    # Eq. 7 is clamped before the weak-melting taper
    probs = params.modified_bourgouin(
        params.BourgouinEnergy(3.0, 0.0, 0.0), 1.0, warm_sfc)
    assert probs.rain == pytest.approx(0.60)
    assert probs.freezing_rain == 0.0
    assert probs.snow == 1.0
    assert probs.ice_pellets == 0.0

    # Freezing rain uses the total melting energy
    probs = params.modified_bourgouin(
        params.BourgouinEnergy(4.0, 2.01, 180.0), 1.0, cold_sfc)
    assert probs.freezing_rain == pytest.approx(0.6464)
    assert probs.rain == 0.0

    # Ice pellets jump to the Eq. 8 value as ME_aloft leaves 0
    energy = params.BourgouinEnergy(10.0, 0.0, 10.0)
    probs = params.modified_bourgouin(energy, 1.0, cold_sfc)
    assert probs.ice_pellets == 0.0
    energy.melting_energy_aloft = np.finfo(np.float32).tiny
    probs = params.modified_bourgouin(energy, 1.0, cold_sfc)
    assert probs.ice_pellets == pytest.approx(0.26)

    # MISSING or NaN in any input gives MISSING
    for probs in (
        params.modified_bourgouin(
            params.BourgouinEnergy(), 1.0, cold_sfc),
        params.modified_bourgouin(
            params.BourgouinEnergy(4.0, 2.01, 180.0), np.nan, cold_sfc),
        params.modified_bourgouin(
            params.BourgouinEnergy(4.0, 2.01, 180.0), 1.0, constants.MISSING),
    ):
        assert probs.rain == constants.MISSING
        assert probs.snow == constants.MISSING
        assert probs.freezing_rain == constants.MISSING
        assert probs.ice_pellets == constants.MISSING


# ---------------------------------------------------------------------------
# Precipitation type from a full sounding
# ---------------------------------------------------------------------------

def saturated_sounding(hght, tmpc):
    """
    Pressure, height, temperature, dewpoint, and wet-bulb arrays (Pa, m, K)
    of a saturated column, where all three temperatures are equal, from
    heights (m) and temperatures (C). Pressure falls 10 Pa per meter from
    1000 hPa.
    """
    hght = np.asarray(hght, dtype='float32')
    tmpk = np.float32(constants.ZEROCNK) + np.asarray(tmpc, dtype='float32')
    pres = np.float32(100000.0) - np.float32(10.0) * hght
    return pres, hght, tmpk, tmpk.copy(), tmpk.copy()


def assert_all_missing(probs):
    assert probs.rain == constants.MISSING
    assert probs.snow == constants.MISSING
    assert probs.freezing_rain == constants.MISSING
    assert probs.ice_pellets == constants.MISSING


def same_probs(a, b):
    return ((a.rain, a.snow, a.freezing_rain, a.ice_pellets) ==
            (b.rain, b.snow, b.freezing_rain, b.ice_pellets))


def test_modified_bourgouin_full_column():
    # Saturated, from the surface up: refreezing below 1400 m, a 2 C warm
    # nose at 1400-1800 m, then cold to -10 C (ProbIce 0.51) at the top.
    snd = saturated_sounding([0, 1400, 1600, 1800, 3000], [-8, 0, 2, 0, -10])
    pres, hght, tmpk, dwpk, wetbulb = snd
    probs = params.modified_bourgouin(pres, hght, tmpk, dwpk, wetbulb)
    assert isinstance(probs, params.PrecipTypeProbabilities)

    # Closed-form energies and ProbIce, through the overload that takes
    # energies
    g_t0 = constants.GRAVITY / constants.ZEROCNK
    melting = g_t0 * 0.5 * 2.0 * 400.0
    refreezing = g_t0 * 0.5 * 8.0 * 1400.0
    expected = params.modified_bourgouin(
        params.BourgouinEnergy(melting, melting, refreezing), 0.51,
        wetbulb[0])
    assert probs.rain == 0.0
    assert probs.snow == pytest.approx(expected.snow, rel=1e-4)
    assert probs.freezing_rain == pytest.approx(
        expected.freezing_rain, rel=1e-4)
    assert probs.ice_pellets == pytest.approx(0.51, rel=1e-4)

    # The steps called by hand give the same result, and so do the
    # defaults written out
    lyr = params.precipitation_generation_layer(pres, hght, tmpk, dwpk)
    tmin, _ = layer.layer_min(lyr, hght, tmpk)
    energy = params.bourgouin_energy(pres, hght, wetbulb)
    composed = params.modified_bourgouin(
        energy, params.probability_of_ice(tmin), wetbulb[0])
    assert same_probs(probs, composed)
    explicit = params.modified_bourgouin(
        pres, hght, tmpk, dwpk, wetbulb, min_depth=0.0, min_energy=0.0,
        pressure_min=25000.0)
    assert same_probs(probs, explicit)

    # The options reach the steps that use them: min_energy = 100 J/kg
    # merges the warm nose into the cold air around it
    probs = params.modified_bourgouin(
        pres, hght, tmpk, dwpk, wetbulb, min_energy=100.0)
    energy = params.bourgouin_energy(pres, hght, wetbulb, min_energy=100.0)
    composed = params.modified_bourgouin(
        energy, params.probability_of_ice(tmin), wetbulb[0])
    assert same_probs(probs, composed)
    assert probs.ice_pellets == 0.0


def test_modified_bourgouin_overloads():
    snd = saturated_sounding([0, 1500, 3000], [-2, -8, -16])
    energy = params.BourgouinEnergy(3.0, 0.0, 0.0)
    warm_sfc = constants.ZEROCNK + 2.0

    # Energies, ProbIce, and the surface wet-bulb: the scalar overload
    for probs in (
        params.modified_bourgouin(energy, 1.0, warm_sfc),
        params.modified_bourgouin(
            energy=energy, prob_ice=1.0, surface_wetbulb=warm_sfc),
    ):
        assert probs.rain == pytest.approx(0.60)

    # Five profiles: the full-column overload. The column is all snow.
    names = ("pressure", "height", "temperature", "dewpoint", "wetbulb")
    for probs in (
        params.modified_bourgouin(*snd),
        params.modified_bourgouin(**dict(zip(names, snd))),
        params.modified_bourgouin(*snd, 0.0, 0.0, 25000.0),
    ):
        assert (probs.rain, probs.snow, probs.freezing_rain,
                probs.ice_pellets) == (0.0, 1.0, 0.0, 0.0)

    with pytest.raises(TypeError):
        params.modified_bourgouin(energy, 1.0)
    with pytest.raises(TypeError):
        params.modified_bourgouin(*snd[:3])
    with pytest.raises(BufferError):
        params.modified_bourgouin(*snd[:4], snd[4][:2])


def test_precip_type_empty_arrays():
    # Every precipitation-type entry point that takes profiles returns
    # MISSING for empty arrays.
    empty = np.array([], dtype='float32')

    energy = params.bourgouin_energy(empty, empty, empty)
    assert energy.melting_energy_total == constants.MISSING
    assert energy.melting_energy_aloft == constants.MISSING
    assert energy.refreezing_energy == constants.MISSING

    lyr = params.precipitation_generation_layer(empty, empty, empty, empty)
    assert lyr.bottom == constants.MISSING
    assert lyr.top == constants.MISSING

    assert_all_missing(
        params.modified_bourgouin(empty, empty, empty, empty, empty))


def load_uppercase_parquet(filename):
    """
    Pressure (Pa), height (m), temperature (K), dewpoint (K), and wet-bulb
    temperature (K) from a 1 Hz sounding in the uppercase schema, which is
    in SI units. Rows without a temperature or dewpoint are dropped, which
    removes the mandatory 1000 hPa level when it is below ground.
    """
    snd_df = pd.read_parquet(filename)
    snd_df = snd_df[snd_df["TEMPERATURE"].notna()]
    snd_df = snd_df[snd_df["DEWPOINT"].notna()]

    pres = snd_df["PRESSURE"].to_numpy().astype('float32')
    hght = snd_df["GEOPOTENTIAL_HEIGHT"].to_numpy().astype('float32')
    tmpk = snd_df["TEMPERATURE"].to_numpy().astype('float32')
    dwpk = snd_df["DEWPOINT"].to_numpy().astype('float32')

    # The preconditions of the precipitation-type functions
    for arr in (pres, hght, tmpk, dwpk):
        assert np.all(np.isfinite(arr))
    assert np.all(np.diff(hght) > 0.0)
    assert np.all(np.diff(pres) < 0.0)

    # The wet-bulb temperature matters only up to the 250 hPa energy cap,
    # so compute it there and leave the rest MISSING. With the Wobus
    # lifter, that takes 0.27 ms for the 1812 levels of 2023-04-19 19Z up
    # to the cap, against 0.46 ms for all 2856 (cm1: 20 ms against 25 ms).
    wetbulb = np.full(tmpk.shape, constants.MISSING, dtype='float32')
    used = pres >= np.float32(25000.0)
    wetbulb[used] = thermo.wetbulb(
        parcel.lifter_wobus(), pres[used], tmpk[used], dwpk[used])

    return pres, hght, tmpk, dwpk, wetbulb


def _layers_numpy(pres, hght, tmpk, dwpk):
    """
    An independent version of the generation-layer rules, with no merging.
    Relative humidity is over ice below 0 C and over liquid otherwise, from
    the Bolton (1980) vapor pressure over liquid and its Magnus-form
    counterpart over ice. Layers are bounded by linearly interpolated 75 %
    crossings.

    Returns the layers from the surface up to the first dry layer deeper
    than 1500 m, as (bottom, top, moist) tuples, and the generation layer:
    the highest moist layer deeper than 1000 m among them, or None.
    """
    tmpc = tmpk.astype('float64') - 273.15
    dwpc = dwpk.astype('float64') - 273.15
    hght = hght.astype('float64')
    vapr = 611.2 * np.exp(17.67 * dwpc / (dwpc + 243.5))
    sat_liquid = 611.2 * np.exp(17.67 * tmpc / (tmpc + 243.5))
    sat_ice = 611.2 * np.exp(21.8745584 * tmpc / (tmpc + 265.49))
    relh = vapr / np.where(tmpc < 0.0, sat_ice, sat_liquid)
    # No level is exactly at 75 %, so no level continues a layer
    assert not np.any(relh == 0.75)

    moist = relh > 0.75
    k = np.nonzero(moist[1:] != moist[:-1])[0]
    crossings = hght[k] + (hght[k + 1] - hght[k]) * \
        (0.75 - relh[k]) / (relh[k + 1] - relh[k])
    edges = np.concatenate(([hght[0]], crossings, [hght[-1]]))
    sides = np.concatenate(([moist[0]], moist[k + 1]))

    layers = []
    generation = None
    for bottom, top, is_moist in zip(edges[:-1], edges[1:], sides):
        layers.append((bottom, top, is_moist))
        if is_moist and top - bottom > 1000.0:
            generation = (bottom, top)
        if not is_moist and top - bottom > 1500.0:
            break
    return layers, generation


def _check_1hz_missing(filename, dry_bottom, dry_top):
    """
    The numpy rules find no generation layer below a deep dry layer from
    about dry_bottom to dry_top (m), so the full column is MISSING, though
    the energies are finite. Returns the pressure profile.
    """
    pres, hght, tmpk, dwpk, wetbulb = load_uppercase_parquet(
        os.path.join(data_dir, filename))

    layers, generation = _layers_numpy(pres, hght, tmpk, dwpk)
    assert generation is None
    bottom, top, moist = layers[-1]
    assert not moist
    assert bottom == pytest.approx(dry_bottom, abs=1.0)
    assert top == pytest.approx(dry_top, abs=1.0)
    # Far from the 1000 m and 1500 m thresholds, so the float32 library
    # and this float64 version agree
    assert top - bottom > 1500.0 + 200.0
    assert max(t - b for b, t, m in layers if m) < 1000.0 - 200.0

    assert_all_missing(params.modified_bourgouin(
        pres, hght, tmpk, dwpk, wetbulb))

    energy = params.bourgouin_energy(pres, hght, wetbulb)
    for value in (energy.melting_energy_total, energy.melting_energy_aloft,
                  energy.refreezing_energy):
        assert np.isfinite(value)
        assert value >= 0.0
    return pres


def test_uppercase_parquet_loader():
    # Each file has one row below ground, the 1000 hPa mandatory level with
    # no temperature or dewpoint. The loader drops it.
    for filename in ("2023-04-19_19_72357.pq", "2023-04-20_00_72357.pq"):
        path = os.path.join(data_dir, filename)
        raw = pd.read_parquet(path)
        missing = raw["TEMPERATURE"].isna() | raw["DEWPOINT"].isna()
        assert missing.sum() == 1
        assert raw["PRESSURE"][missing].item() == 100000.0
        assert (raw["GEOPOTENTIAL_HEIGHT"][missing].item() <
                raw["STATION_HEIGHT"].iloc[0])

        pres, hght, tmpk, dwpk, wetbulb = load_uppercase_parquet(path)
        assert len(hght) == len(raw) - 1
        assert hght[0] == raw["STATION_HEIGHT"].iloc[0]


def test_modified_bourgouin_1hz_2023_04_19_19z():
    # A 698 m moist layer at 1149-1847 m under a 7.7 km dry layer
    _check_1hz_missing("2023-04-19_19_72357.pq", 1846.7, 9508.9)


def test_modified_bourgouin_1hz_2023_04_20_00z():
    # Smoke test: a 578 m moist layer at 937-1514 m under a 1751 m dry
    # layer, which eliminates the deep moist layer above it. Every level is
    # at or above 250 hPa.
    pres = _check_1hz_missing("2023-04-20_00_72357.pq", 1514.2, 3265.2)
    assert pres[-1] >= 25000.0


# ===========================================================================
# Precipitation type: the spectral bin classifier (Reeves et al. 2016)
# ===========================================================================

# ---------------------------------------------------------------------------
# Result types and the drop-size distribution
# ---------------------------------------------------------------------------

def test_precip_type_encoding():
    assert {member.name: int(member) for member in params.PrecipType} == {
        "missing": -9999,
        "rain": 1,
        "snow": 2,
        "rain_snow": 3,
        "freezing_rain": 4,
        "ice_pellets": 5,
        "freezing_rain_ice_pellets": 6,
        "rain_ice_pellets": 7,
    }
    assert float(params.PrecipType.missing) == constants.MISSING
    assert params.PrecipType(4) == params.PrecipType.freezing_rain


def test_spectral_bin_result():
    result = params.SpectralBinResult()
    assert result.precip_type == params.PrecipType.missing
    assert result.liquid_fraction == constants.MISSING
    assert result.supercooled_liquid_height == constants.MISSING

    for field in ("precip_type", "liquid_fraction",
                  "supercooled_liquid_height"):
        with pytest.raises(AttributeError):
            setattr(result, field, 1.0)


def test_spectral_bin_classifier_constants():
    assert params.SBC_MAX_BINS == 64
    assert params.SBC_ICE_NUCLEATION_TEMPERATURE == np.float32(267.15)


def test_spectral_bin_dsd_default():
    dsd = params.spectral_bin_dsd_default()
    assert dsd.nbins == 4
    assert dsd.rime_factor == 1.0
    assert dsd.diameter.dtype == np.float32
    assert dsd.concentration.dtype == np.float32
    np.testing.assert_array_equal(
        dsd.diameter, np.array([0.05, 0.75, 1.45, 2.15], dtype="float32"))
    np.testing.assert_array_equal(
        dsd.concentration,
        np.array([55.1843, 146.647, 11.6891, 3.60886], dtype="float32"))


def python_reference_dsd(deld):
    # run_sbc.py, cast to float32
    psd_orig = [55.1843, 66.0695, 130.272, 154.556, 203.649, 171.814,
                206.606, 146.647, 94.9404, 79.4013, 61.0083, 35.6567,
                25.4924, 16.2522, 11.6891, 7.49152, 3.60886]
    diameter_orig = [0.05 + 0.1 * i for i in range(len(psd_orig))]
    diameter = np.arange(0.05, 1.85 + deld, deld)
    concentration = np.interp(diameter, diameter_orig, psd_orig)
    return diameter.astype("float32"), concentration.astype("float32")


def test_spectral_bin_dsd():
    cxx_diameter = np.array([0.05, 0.65, 1.25, 1.85], dtype="float32")
    cxx_concentration = np.array([55.1843, 206.606, 25.4924, 3.60886],
                                 dtype="float32")

    for diameter, concentration in (
        python_reference_dsd(0.7),
        python_reference_dsd(0.1),
        (cxx_diameter, cxx_concentration),
    ):
        for rime_factor in (1.0, 5.0):
            dsd = params.spectral_bin_dsd(diameter, concentration,
                                          rime_factor=rime_factor)
            assert dsd.nbins == diameter.size
            assert dsd.rime_factor == rime_factor
            np.testing.assert_array_equal(dsd.diameter, diameter)
            np.testing.assert_array_equal(dsd.concentration, concentration)

    default = params.spectral_bin_dsd_default()
    diameter, concentration = python_reference_dsd(0.7)
    np.testing.assert_array_equal(default.diameter, diameter)
    np.testing.assert_array_equal(default.concentration, concentration)
    assert python_reference_dsd(0.1)[0].size == 19

    dsd = params.spectral_bin_dsd(cxx_diameter, cxx_concentration)
    assert dsd.rime_factor == 1.0
    for field in ("nbins", "rime_factor", "diameter", "concentration"):
        with pytest.raises(AttributeError):
            setattr(dsd, field, 1.0)

    diameter = (0.05 + 0.1 * np.arange(params.SBC_MAX_BINS + 1)).astype(
        "float32")
    concentration = np.ones(diameter.size, dtype="float32")
    assert params.spectral_bin_dsd(
        diameter[:-1], concentration[:-1]).nbins == params.SBC_MAX_BINS
    assert params.spectral_bin_dsd(diameter, concentration).nbins == 0


@pytest.mark.parametrize("diameter, concentration, rime_factor", [
    ([], [], 1.0),
    ([0.05, 0.65, 0.65, 1.85], [55.1843, 206.606, 25.4924, 3.60886], 1.0),
    ([0.0, 0.65, 1.25, 1.85], [55.1843, 206.606, 25.4924, 3.60886], 1.0),
    ([0.05, np.nan, 1.25, 1.85], [55.1843, 206.606, 25.4924, 3.60886], 1.0),
    ([0.05, 0.65, 1.25, 12.16], [55.1843, 206.606, 25.4924, 3.60886], 1.0),
    ([0.05, 0.65, 1.25, 1.85], [55.1843, -206.606, 25.4924, 3.60886], 1.0),
    ([0.05, 0.65, 1.25, 1.85], [0.0, 0.0, 0.0, 0.0], 1.0),
    ([0.05, 0.65, 1.25, 1.85], [55.1843, np.nan, 25.4924, 3.60886], 1.0),
    ([0.05, 0.65, 1.25, 1.85], [55.1843, 206.606, 25.4924, 3.60886], 0.99),
    ([0.05, 0.65, 1.25, 1.85], [55.1843, 206.606, 25.4924, 3.60886], 5.01),
    ([0.05, 0.65, 1.25, 1.85], [55.1843, 206.606, 25.4924, 3.60886], np.nan),
])
def test_spectral_bin_dsd_invalid(diameter, concentration, rime_factor):
    dsd = params.spectral_bin_dsd(np.array(diameter, dtype="float32"),
                                  np.array(concentration, dtype="float32"),
                                  rime_factor)
    assert dsd.nbins == 0
    assert dsd.rime_factor == constants.MISSING
    assert dsd.diameter.size == 0
    assert dsd.concentration.size == 0


def test_spectral_bin_dsd_sizes():
    diameter = np.array([0.05, 0.65, 1.25, 1.85], dtype="float32")
    with pytest.raises(BufferError):
        params.spectral_bin_dsd(diameter, diameter[:3])


# ---------------------------------------------------------------------------
# Cloud top from a sounding
# ---------------------------------------------------------------------------

sbc_reference_dir = os.path.join(
    os.path.dirname(__file__), "..", "..", "data", "sbc_reference")


def test_spectral_bin_cloud_top_golden():
    cases = pd.read_parquet(os.path.join(sbc_reference_dir, "cases.parquet"))
    levels = pd.read_parquet(
        os.path.join(sbc_reference_dir, "levels.parquet"))
    assert set(cases["group"]) == {"named", "sample", "corpus"}

    pres, hght, tmpk, dwpk, relh = (
        levels[name].to_numpy()
        for name in ("pressure", "height", "temperature", "dewpoint", "relh"))
    case_id = cases["case_id"].to_numpy()
    starts = np.searchsorted(levels["case_id"].to_numpy(), case_id)
    stops = np.searchsorted(levels["case_id"].to_numpy(), case_id,
                            side="right")
    assert np.all(stops > starts)

    tops = np.array([
        params.spectral_bin_cloud_top(pres[a:b], hght[a:b], tmpk[a:b],
                                      dwpk[a:b], relh[a:b])
        for a, b in zip(starts, stops)
    ], dtype="float32")
    np.testing.assert_array_equal(tops, cases["cloud_top_height"].to_numpy())

    level = cases["cloud_top_level"].to_numpy()
    cloud = level != constants.MISSING
    assert np.all(cloud == (tops != constants.MISSING))
    np.testing.assert_array_equal(tops[cloud],
                                  hght[starts[cloud] + level[cloud]])


def test_spectral_bin_cloud_top():
    hght = np.array([0.0, 1000.0, 2000.0], dtype="float32")
    pres = np.array([100000.0, 88250.0, 77880.0], dtype="float32")
    tmpk = np.full(3, 270.0, dtype="float32")
    dwpk = tmpk - np.array([1.0, 15.0, 2.0], dtype="float32")
    relh = np.array([0.9, 0.2, 0.9], dtype="float32")

    top = params.spectral_bin_cloud_top(pres, hght, tmpk, dwpk, relh)
    assert isinstance(top, float)
    assert top == 0.0

    missing = tmpk.copy()
    missing[0] = constants.MISSING
    assert params.spectral_bin_cloud_top(
        pres, hght, missing, dwpk, relh) == 2000.0
    missing = relh.copy()
    missing[2] = np.nan
    assert params.spectral_bin_cloud_top(
        pres, hght, tmpk, dwpk, missing) == 0.0

    dry = tmpk - np.float32(8.0)
    assert params.spectral_bin_cloud_top(
        pres, hght, tmpk, dry, relh) == 0.0
    assert params.spectral_bin_cloud_top(
        pres, hght, tmpk, dry, relh / np.float32(2.0)) == constants.MISSING

    empty = np.array([], dtype="float32")
    assert params.spectral_bin_cloud_top(
        empty, empty, empty, empty, empty) == constants.MISSING

    with pytest.raises(BufferError):
        params.spectral_bin_cloud_top(pres, hght, tmpk, dwpk, relh[:2])


# ---------------------------------------------------------------------------
# Precipitation type from a given cloud top: pre-classifier
# ---------------------------------------------------------------------------

sbc_reference_dir = os.path.join(
    os.path.dirname(__file__), "..", "..", "data", "sbc_reference")


@functools.cache
def sbc_reference():
    """
    The golden data of the spectral bin classifier: the cases table, and per
    case_id the level arrays (surface up) and the reference profile
    (levels x bins).
    """
    def read(name):
        return pd.read_parquet(os.path.join(sbc_reference_dir, name))

    cases = read("cases.parquet")
    levels = {
        case_id: tuple(group[field].to_numpy() for field in (
            "pressure", "height", "temperature", "dewpoint", "relh",
            "wetbulb"))
        for case_id, group in read("levels.parquet").groupby("case_id")
    }
    profiles = {
        case_id: group["liquid_fraction"].to_numpy().reshape(
            group["level"].nunique(), group["bin"].nunique())
        for case_id, group in read("profiles.parquet").groupby("case_id")
    }
    dsds = {
        name: (group["diameter"].to_numpy(), group["concentration"].to_numpy())
        for name, group in read("dsds.parquet").groupby("dsd_name")
    }
    return cases, levels, profiles, dsds


def sbc_reference_runs(stages):
    """
    Each golden case of the given stages, run with its reference cloud top,
    as (case, result, profile).
    """
    cases, levels, _, dsds = sbc_reference()
    for case in cases[cases["stage"].isin(stages)].itertuples():
        dsd = params.spectral_bin_dsd(*dsds[case.dsd_name],
                                      rime_factor=case.rime_factor)
        result, profile = params.spectral_bin_classifier(
            *levels[case.case_id], case.cloud_top_height, dsd,
            case.ice_nucleation_temperature, return_profile=True)
        without = params.spectral_bin_classifier(
            *levels[case.case_id], case.cloud_top_height, dsd,
            case.ice_nucleation_temperature)
        assert sbc_tuple(without) == sbc_tuple(result)
        yield case, result, profile


def sbc_tuple(result):
    return (result.precip_type, result.liquid_fraction,
            result.supercooled_liquid_height)


SBC_SN = (params.PrecipType.snow, 0.0, constants.MISSING)
SBC_FZRA = (params.PrecipType.freezing_rain, 1.0, 0.0)
SBC_RA = (params.PrecipType.rain, 1.0, constants.MISSING)
SBC_MISSING = (params.PrecipType.missing, constants.MISSING,
               constants.MISSING)


def test_spectral_bin_classifier_reference():
    _, _, profiles, _ = sbc_reference()
    counts = {}
    for case, result, profile in sbc_reference_runs(
            ("preclassifier", "no_cloud")):
        key = (case.group, case.stage)
        counts[key] = counts.get(key, 0) + 1
        assert sbc_tuple(result) == (
            params.PrecipType(case.precip_type), case.liquid_fraction,
            case.supercooled_liquid_height), case.case_id
        np.testing.assert_array_equal(profile, profiles[case.case_id])
        assert np.all(profile == constants.MISSING)
    assert counts == {
        ("named", "preclassifier"): 6,
        ("named", "no_cloud"): 1,
        ("corpus", "preclassifier"): 1207,
        ("corpus", "no_cloud"): 35,
    }


def saturated_sbc_profile(height, wetbulb):
    height = np.array(height, dtype="float32")
    wetbulb = np.array(wetbulb, dtype="float32")
    pressure = (100000.0 * np.exp(-height / 8000.0)).astype("float32")
    return (pressure, height, wetbulb, wetbulb.copy(),
            np.ones_like(wetbulb), wetbulb.copy())


def run_sbc(snd, cloud_top, dsd=None,
            tice=params.SBC_ICE_NUCLEATION_TEMPERATURE):
    """
    The result as a tuple, after checking that the profile is all MISSING
    and that asking for it does not change the result.
    """
    if dsd is None:
        dsd = params.spectral_bin_dsd_default()
    result, profile = params.spectral_bin_classifier(
        *snd, cloud_top, dsd, tice, return_profile=True)
    assert profile.dtype == np.float32
    assert profile.shape == (snd[1].size, dsd.nbins)
    assert np.all(profile == constants.MISSING)
    without = params.spectral_bin_classifier(*snd, cloud_top, dsd, tice)
    assert sbc_tuple(without) == sbc_tuple(result)
    return sbc_tuple(result)


SBC_TICE = float(params.SBC_ICE_NUCLEATION_TEMPERATURE)
SBC_HEIGHT_5 = [0.0, 1000.0, 2000.0, 3000.0, 4000.0]
SBC_HEIGHT_4 = SBC_HEIGHT_5[:4]
SBC_TIES_TW = [272.15, 271.15, 270.15, 268.15, 265.15, 262.15]


@pytest.mark.parametrize("height, wetbulb, cloud_top, expected", [
    (SBC_HEIGHT_5, [272.15, 273.15, 270.15, 265.15, 262.15], 4000.0, SBC_SN),
    (SBC_HEIGHT_4, [273.1, 275.15, 276.15, 274.15], 3000.0, SBC_FZRA),
    (SBC_HEIGHT_5, [270.15, 268.15, 262.15, 260.15, SBC_TICE], 4000.0,
     SBC_FZRA),
    (SBC_HEIGHT_5, [270.15, 268.15, 262.15, 260.15, 267.1], 4000.0, SBC_SN),
    (SBC_HEIGHT_5, [270.15, 268.15, SBC_TICE, 269.15, 262.15], 4000.0,
     SBC_FZRA),
    (SBC_HEIGHT_5, [270.15, 268.15, 267.1, 269.15, 262.15], 4000.0, SBC_SN),
    (SBC_HEIGHT_4, [270.15, 275.15, 276.15, 267.2], 3000.0, SBC_FZRA),
])
def test_spectral_bin_classifier_thresholds(height, wetbulb, cloud_top,
                                            expected):
    assert run_sbc(saturated_sbc_profile(height, wetbulb),
                   cloud_top) == expected


@pytest.mark.parametrize("height, wetbulb, cloud_top, expected", [
    ([0.0, 1000.0, 2000.0, 2900.0, 3100.0, 4000.0], SBC_TIES_TW, 4000.0,
     SBC_SN),
    ([0.0, 1000.0, 2000.0, 2901.0, 3100.0, 4000.0], SBC_TIES_TW, 4000.0,
     SBC_FZRA),
    ([0.0, 1000.0, 2000.0, 2900.0, 3099.0, 4000.0], SBC_TIES_TW, 4000.0,
     SBC_SN),
    ([1000.0, 2000.0, 3000.0, 3900.0, 4100.0, 5000.0], SBC_TIES_TW, 5000.0,
     SBC_SN),
    ([1000.0, 2000.0, 3000.0, 3901.0, 4100.0, 5000.0], SBC_TIES_TW, 5000.0,
     SBC_FZRA),
    ([0.0, 6000.0], [270.15, 262.15], 6000.0, SBC_SN),
    ([0.0, 6001.0], [270.15, 262.15], 6001.0, SBC_FZRA),
    ([0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0, 3500.0],
     [270.15, 269.15, 268.15, 268.15, 266.15, 275.15, 262.15, 260.15],
     2000.0, SBC_SN),
])
def test_spectral_bin_classifier_3km_window(height, wetbulb, cloud_top,
                                            expected):
    assert run_sbc(saturated_sbc_profile(height, wetbulb),
                   cloud_top) == expected


@pytest.mark.parametrize("wetbulb, tice, expected", [
    ([270.15, 268.15, 266.15, 264.15, 262.15], 267.15, SBC_SN),
    ([270.15, 268.15, 266.15, 264.15, 262.15], 263.15, SBC_FZRA),
    ([270.15, 268.15, 266.15, 262.15, 265.15], 267.15, SBC_SN),
    ([270.15, 268.15, 266.15, 262.15, 265.15], 263.15, SBC_FZRA),
])
def test_spectral_bin_classifier_tice(wetbulb, tice, expected):
    snd = saturated_sbc_profile(SBC_HEIGHT_5, wetbulb)
    assert run_sbc(snd, 4000.0, tice=tice) == expected


@pytest.mark.parametrize("surface", [0.0, 1500.0])
@pytest.mark.parametrize("cloud_top, expected", [
    (1000.0, SBC_FZRA),
    (1500.0, SBC_FZRA),
    (2000.0, SBC_SN),
    (2999.0, SBC_SN),
    (3000.0, SBC_FZRA),
    (3500.0, SBC_FZRA),
    (4000.0, SBC_SN),
    (20000.0, SBC_SN),
    (999.0, SBC_MISSING),
    (0.0, SBC_MISSING),
    (-1.0, SBC_MISSING),
])
def test_spectral_bin_classifier_cloud_top(surface, cloud_top, expected):
    snd = saturated_sbc_profile(
        np.array(SBC_HEIGHT_5) + surface,
        [272.15, 270.15, 266.15, 268.15, 262.15])
    assert run_sbc(snd, surface + cloud_top) == expected


def test_spectral_bin_classifier_invalid():
    snd = saturated_sbc_profile(SBC_HEIGHT_5,
                                [272.15, 270.15, 266.15, 268.15, 262.15])
    assert run_sbc(snd, 4000.0) == SBC_SN

    for cloud_top in (constants.MISSING, np.nan):
        assert run_sbc(snd, cloud_top) == SBC_MISSING

    empty = np.array([], dtype="float32")
    invalid_dsd = params.spectral_bin_dsd(empty, empty)
    assert run_sbc(snd, 4000.0, dsd=invalid_dsd) == SBC_MISSING

    for tice in (constants.MISSING, np.nan, 0.0, -1.0):
        assert run_sbc(snd, 4000.0, tice=tice) == SBC_MISSING

    for N in (0, 1):
        assert run_sbc(tuple(arr[:N] for arr in snd), 4000.0) == SBC_MISSING

    with pytest.raises(BufferError):
        params.spectral_bin_classifier(*snd[:5], snd[5][:4], 4000.0)


def test_spectral_bin_classifier_defaults():
    snd = saturated_sbc_profile(SBC_HEIGHT_5,
                                [270.15, 268.15, 266.15, 264.15, 262.15])
    result = params.spectral_bin_classifier(*snd, 4000.0)
    assert isinstance(result, params.SpectralBinResult)
    assert sbc_tuple(result) == SBC_SN
    assert sbc_tuple(params.spectral_bin_classifier(
        *snd, cloud_top=4000.0, ice_nucleation_temperature=263.15)) == SBC_FZRA

    deld_0_1 = params.spectral_bin_dsd(*python_reference_dsd(0.1))
    result, profile = params.spectral_bin_classifier(
        *snd, 4000.0, deld_0_1, return_profile=True)
    assert sbc_tuple(result) == SBC_SN
    assert profile.shape == (5, 19)

    signature = params.spectral_bin_classifier.__doc__.splitlines()[0]
    assert "dsd: nwsspc.sharp.calc.params.SpectralBinDSD = " \
        "spectral_bin_dsd_default()" in signature
    assert "ice_nucleation_temperature: float = " \
        "SBC_ICE_NUCLEATION_TEMPERATURE" in signature
    assert "return_profile: bool = False" in signature


@pytest.mark.parametrize("field", [2, 3, 4, 5])
@pytest.mark.parametrize("bad", [constants.MISSING, np.nan])
def test_spectral_bin_classifier_missing_levels(field, bad):
    def spoil(snd, level):
        snd = tuple(arr.copy() for arr in snd)
        snd[field][level] = bad
        return snd

    window = saturated_sbc_profile(
        [0.0, 1000.0, 2500.0, 3000.0, 4000.0, 5000.0],
        [272.15, 270.15, 268.15, 264.15, 262.15, 260.15])
    assert run_sbc(window, 5000.0) == SBC_SN
    assert run_sbc(spoil(window, 3), 5000.0) == SBC_FZRA

    surface = saturated_sbc_profile(
        [0.0, 500.0, 1000.0, 3200.0, 3600.0, 5000.0],
        [272.15, 271.15, 270.15, 268.15, 264.15, 260.15])
    assert run_sbc(surface, 5000.0) == SBC_FZRA
    assert run_sbc(spoil(surface, 0), 5000.0) == SBC_SN

    top = saturated_sbc_profile(SBC_HEIGHT_5,
                                [272.15, 270.15, 266.15, 268.15, 262.15])
    for cloud_top in (4000.0, 20000.0):
        assert run_sbc(top, cloud_top) == SBC_SN
        assert run_sbc(spoil(top, 4), cloud_top) == SBC_FZRA

    two = saturated_sbc_profile(SBC_HEIGHT_5[:3], [272.15, 270.15, 266.15])
    assert run_sbc(two, 2000.0) == SBC_SN
    assert run_sbc(spoil(two, 0), 2000.0) == SBC_SN
    assert run_sbc(spoil(two, 0), 1000.0) == SBC_MISSING
    assert run_sbc(spoil(spoil(two, 0), 1), 2000.0) == SBC_MISSING
    assert run_sbc(spoil(spoil(two, 1), 2), 2000.0) == SBC_MISSING


# ---------------------------------------------------------------------------
# Microphysics: frozen cloud tops and melting
# ---------------------------------------------------------------------------

SBC_THRESHOLDS = np.array([0.15, 0.60, 0.85])


def sbc_decision(liquid_fraction, crossings, warm):
    """
    The reference's decision tree on liquid fractions, with each case's 0 C
    crossings and surface class (Tw > 273.15 K), as PrecipType codes. It
    compares in float32, as the classifier does.
    """
    lf = np.asarray(liquid_fraction, dtype=np.float32)
    f32 = np.float32
    ice, liquid = lf == 0, lf == 1
    one_warm = np.where(liquid | (lf > f32(0.85)), 1,
                        np.where(ice | (lf < f32(0.60)), 2, 3))
    many_warm = np.where(ice | (lf < f32(0.15)), 5,
                         np.where(liquid | (f32(1) - lf < f32(0.15)), 1, 7))
    cold = np.where(ice | (lf < f32(0.15)), 5,
                    np.where(liquid | (lf > f32(0.85)), 4, 6))
    return np.where(warm, np.where(np.asarray(crossings) == 1, one_warm,
                                   many_warm), cold)


def sbc_same_level(agl, height, expected):
    """Whether two heights AGL are both MISSING, or both the same level."""
    if (height == constants.MISSING) or (expected == constants.MISSING):
        return height == expected
    level, expected_level = np.abs(agl - height), np.abs(agl - expected)
    return bool(level.argmin() == expected_level.argmin()
                and level.min() < 0.01 and expected_level.min() < 0.01)


def sbc_run_composed(case, snd, dsd):
    """The composed overload, from the case's reference cloud top."""
    return params.spectral_bin_classifier(
        *snd, case.cloud_top_height, dsd, case.ice_nucleation_temperature,
        return_profile=True)


def sbc_branches_within(letters):
    """Selects the core cases that take only the given branches."""
    return lambda cases: (cases["stage"] == "core") & cases["branches"].map(
        lambda branches: set(branches) <= set(letters))


def sbc_golden_check(select, run=sbc_run_composed):
    """
    Checks the golden cases that select(cases) picks under the golden-data
    comparison rules. run(case, snd, dsd) returns the result and the
    profile. A case without flags must pass every rule, and every case rule
    1. A near_discontinuity corpus case that breaks rules 2-4 is listed, and
    the listed cases stay at most 1 % of the corpus. Returns the selected
    cases with the returned values, the listed case_ids, and the largest
    liquid-fraction and profile errors outside the exemptions.
    """
    cases, levels, profiles, dsds = sbc_reference()
    chosen = cases[select(cases)].reset_index(drop=True)
    n = len(chosen)
    assert n > 0
    precip_type = np.empty(n, dtype=int)
    liquid_fraction = np.empty(n)
    slw_height = np.empty(n)
    same_slw_level = np.empty(n, dtype=bool)
    warm = np.empty(n, dtype=bool)
    ours, theirs = [], []
    for i, case in enumerate(chosen.itertuples()):
        snd = levels[case.case_id]
        dsd = params.spectral_bin_dsd(*dsds[case.dsd_name],
                                      rime_factor=case.rime_factor)
        result, profile = run(case, snd, dsd)
        precip_type[i] = int(result.precip_type)
        liquid_fraction[i] = result.liquid_fraction
        slw_height[i] = result.supercooled_liquid_height
        same_slw_level[i] = sbc_same_level(
            snd[1] - snd[1][0], result.supercooled_liquid_height,
            case.supercooled_liquid_height)
        warm[i] = snd[5][0] > np.float32(273.15)
        ours.append(profile.ravel())
        theirs.append(profiles[case.case_id].ravel())

    core = (chosen["stage"] == "core").to_numpy()
    corpus = (chosen["group"] == "corpus").to_numpy()
    ref_type = chosen["precip_type"].to_numpy()
    ref_lf = chosen["liquid_fraction"].to_numpy()
    ref_slw = chosen["supercooled_liquid_height"].to_numpy()

    exact = ((precip_type == ref_type) & (liquid_fraction == ref_lf)
             & (slw_height == ref_slw))
    consistent = precip_type == sbc_decision(
        liquid_fraction, chosen["crossings"].to_numpy(), warm)
    rule_1 = np.where(core, consistent, exact)

    nearest = SBC_THRESHOLDS[np.abs(ref_lf[:, None]
                                    - SBC_THRESHOLDS).argmin(axis=1)]
    other_side = sbc_decision(2 * nearest - ref_lf,
                              chosen["crossings"].to_numpy(), warm)
    rule_2 = (precip_type == ref_type) | (
        corpus & chosen["near_threshold"].to_numpy()
        & (precip_type == other_side))

    sizes = np.array([p.size for p in ours])
    case_of = np.repeat(np.arange(n), sizes)
    ours, theirs = np.concatenate(ours), np.concatenate(theirs)
    nbins = np.repeat(chosen["dsd_name"].map(
        lambda name: dsds[name][0].size).to_numpy(), sizes)
    position = np.arange(ours.size) - np.repeat(np.cumsum(sizes) - sizes,
                                                sizes)
    level, bin_ = position // nbins, position % nbins
    missing = theirs == constants.MISSING
    same_missing = (ours == constants.MISSING) == missing
    error = np.where(missing, 0.0, np.abs(ours - theirs))
    disc = chosen["near_discontinuity"].to_numpy()[case_of]
    exempt = disc & (level <= chosen["disc_level"].to_numpy()[case_of]) & (
        (chosen["disc_scope"].to_numpy()[case_of] == "column")
        | (bin_ == chosen["disc_bin"].to_numpy()[case_of]))

    def per_case(bad):
        return np.bincount(case_of, weights=bad, minlength=n) == 0

    rule_3 = ((np.abs(liquid_fraction - ref_lf) <= 0.005) & same_slw_level
              & per_case(~same_missing))
    rule_4 = per_case((error > 1e-3) & ~exempt)

    flagged = corpus & chosen["near_discontinuity"].to_numpy()
    broken = ~(rule_2 & rule_3 & rule_4)
    failing = chosen.loc[~rule_1 | (broken & ~flagged), "case_id"].tolist()
    assert not failing, f"cases that break the comparison rules: {failing}"
    listed = chosen.loc[broken & flagged, "case_id"].tolist()
    assert len(listed) <= 0.01 * (cases["group"] == "corpus").sum(), listed

    return {
        "cases": chosen.assign(result_precip_type=precip_type,
                               result_liquid_fraction=liquid_fraction,
                               result_supercooled_liquid_height=slw_height),
        "listed": listed,
        "liquid_fraction_error": float(np.abs(liquid_fraction
                                              - ref_lf)[core].max(initial=0)),
        "profile_error": float(error[~exempt].max(initial=0)),
    }


def test_spectral_bin_classifier_golden_frozen_tops():
    check = sbc_golden_check(sbc_branches_within("ABF"))
    cases = check["cases"]
    assert cases["group"].value_counts().to_dict() == {"corpus": 105,
                                                       "named": 8}
    assert set(cases["precip_type"]) == {1, 2, 3}
    assert set(cases["crossings"]) == {1}
    _, levels, _, _ = sbc_reference()
    assert all(levels[case_id][5][0] > np.float32(273.15)
               for case_id in cases["case_id"])
    assert set(cases["dsd_name"]) == {"python_default", "cpp_2.0.3",
                                      "python_deld_0.1"}
    assert set(cases["rime_factor"]) == {1.0, 5.0}
    assert set(cases["ice_nucleation_temperature"]) == {
        np.float32(263.15), np.float32(267.15)}
    assert check["listed"] == []


# ---------------------------------------------------------------------------
# Microphysics: refreezing
# ---------------------------------------------------------------------------

def test_spectral_bin_classifier_golden_refreezing():
    check = sbc_golden_check(sbc_branches_within("ABFG"))
    cases = check["cases"]
    refreezing = cases[cases["branches"].str.contains("G")]
    assert refreezing["group"].value_counts().to_dict() == {
        "corpus": 1001, "named": 7, "sample": 1}
    assert set(refreezing["precip_type"]) == {1, 4, 5, 6, 7}
    assert refreezing["tnuc_switched"].sum() == 642
    assert refreezing["refreeze_level_moved"].sum() == 249
    assert check["listed"] == []


def test_spectral_bin_classifier_golden_refreezing_switches():
    check = sbc_golden_check(lambda cases: sbc_branches_within("ABFG")(cases)
                             & (cases["group"] != "corpus")
                             & (cases["tnuc_switched"]
                                | cases["refreeze_level_moved"]))
    assert check["cases"]["case_id"].tolist() == [8, 10, 11, 14, 19, 1000]


def test_spectral_bin_classifier_golden_sample():
    check = sbc_golden_check(lambda cases: cases["group"] == "sample")
    sample = check["cases"].iloc[0]
    assert sample["result_precip_type"] == int(params.PrecipType.ice_pellets)
    assert sample["result_liquid_fraction"] == 0.0
    assert sample["result_supercooled_liquid_height"] == np.float32(1130.5469)


def test_spectral_bin_classifier_rule_2_at_tice():
    snd = saturated_sbc_profile(SBC_HEIGHT_4,
                                [270.15, 275.15, 276.15, SBC_TICE])
    result, profile = params.spectral_bin_classifier(*snd, 3000.0,
                                                     return_profile=True)
    assert sbc_tuple(result) == SBC_FZRA
    np.testing.assert_array_equal(profile, [[1.0] * 4] * 3 + [[0.0] * 4])


# ---------------------------------------------------------------------------
# Microphysics: liquid cloud tops
# ---------------------------------------------------------------------------

def test_spectral_bin_classifier_golden_liquid_tops():
    check = sbc_golden_check(sbc_branches_within("CDF"))
    cases = check["cases"]
    assert cases.groupby(["group", "branches"]).size().to_dict() == {
        ("corpus", "CDF"): 8, ("corpus", "CF"): 5, ("named", "CDF"): 1}
    _, levels, _, _ = sbc_reference()
    for case in cases.itertuples():
        wetbulb = levels[case.case_id][5]
        assert case.ice_nucleation_temperature < wetbulb[
            case.cloud_top_level] <= np.float32(273.15)
        assert wetbulb[0] > np.float32(273.15)
    assert set(cases["dsd_name"]) == {"python_default", "cpp_2.0.3",
                                      "python_deld_0.1"}
    assert set(cases["ice_nucleation_temperature"]) == {
        np.float32(263.15), np.float32(267.15)}
    assert check["listed"] == []


@pytest.mark.parametrize("wetbulb, expected", [
    ([273.15, 275.15, 276.15, 274.15], SBC_FZRA),
    ([278.15, 273.15, 275.15, 274.15], (params.PrecipType.rain, 1.0, 3000.0)),
])
def test_spectral_bin_classifier_warm_top_at_0c(wetbulb, expected):
    result, profile = params.spectral_bin_classifier(
        *saturated_sbc_profile(SBC_HEIGHT_4, wetbulb), 3000.0,
        params.spectral_bin_dsd_default(), return_profile=True)
    assert sbc_tuple(result) == expected
    np.testing.assert_array_equal(profile, np.ones((4, 4)))


# ---------------------------------------------------------------------------
# Precipitation type from a full sounding
# ---------------------------------------------------------------------------
