
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

    # turn into height above ground level, keeping the heights MSL as well
    hght_msl = hght.copy()
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

    assert (storm_mtn.u == pytest.approx(9.701575))
    assert (storm_mtn.v == pytest.approx(5.622299))


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


def ebwd_from_definition(pres, hght, uwin, vwin, eil, eql_pres):
    """The effective bulk wind difference from its definition, in numpy.

    The layer runs from the effective inflow base up half the distance from
    that base to the EL (Thompson et al. 2007), in meters AGL. Height is
    linear in log10(pressure), and wind is linear in height. Returns the
    layer bottom and top (m AGL) and the wind difference (u, v).
    """
    logp = np.log10(pres.astype("float64"))[::-1]
    z = hght.astype("float64")

    def z_agl(p):
        return np.interp(np.log10(p), logp, z[::-1]) - z[0]

    def wind(agl):
        return (np.interp(z[0] + agl, z, uwin), np.interp(z[0] + agl, z, vwin))

    bot = z_agl(eil.bottom)
    top = bot + 0.5 * (z_agl(eql_pres) - bot)
    (u_bot, v_bot), (u_top, v_top) = wind(bot), wind(top)
    return bot, top, (u_top - u_bot, v_top - v_bot)


# The sounding with heights AGL, as the other tests read it, and MSL, as the
# file has it, with the station at 790 m. Before SHARPlib-27b the AGL heights
# gave (14.6, 13.3217), as now. The MSL heights gave (16.7517, 5.89825), the
# shear over 790 to 6656.09 m AGL, a layer 790 m too high. The definition
# gives (14.6, 13.32178) over 0 to 5866.09 m AGL: the inflow layer starts at
# the ground and the EL is at 11732.18 m AGL.
@pytest.mark.parametrize("hght_key", ["hght", "hght_msl"])
def test_effective_bulk_wind(hght_key):
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    mupcl = parcel.Parcel()
    eil = params.effective_inflow_layer(
        lifter,
        snd_data["pres"],
        snd_data[hght_key],
        snd_data["tmpk"],
        snd_data["dwpk"],
        snd_data["vtmp"],
        mupcl=mupcl
    )

    ebwd_cmp = params.effective_bulk_wind_difference(
        snd_data["pres"],
        snd_data[hght_key],
        snd_data["uwin"],
        snd_data["vwin"],
        eil,
        mupcl.eql_pressure
    )

    bot, top, expected = ebwd_from_definition(
        snd_data["pres"], snd_data[hght_key], snd_data["uwin"],
        snd_data["vwin"], eil, mupcl.eql_pressure)
    assert (bot == 0.0)
    assert (top == pytest.approx(5866.089, abs=1e-3))
    assert (expected[0] == pytest.approx(14.6, abs=1e-5))
    assert (expected[1] == pytest.approx(13.32178, abs=1e-5))

    ebwd = winds.vector_magnitude(ebwd_cmp.u, ebwd_cmp.v)
    assert (ebwd_cmp.u == pytest.approx(expected[0], abs=1e-4))
    assert (ebwd_cmp.v == pytest.approx(expected[1], abs=1e-4))
    assert (ebwd == pytest.approx(np.hypot(*expected), abs=1e-4))

    # the same as wind_shear over the AGL layer from the definition
    shear = winds.wind_shear(layer.HeightLayer(bot, top), snd_data[hght_key],
                             snd_data["uwin"], snd_data["vwin"])
    assert (ebwd_cmp.u == pytest.approx(shear.u, abs=1e-4))
    assert (ebwd_cmp.v == pytest.approx(shear.v, abs=1e-4))


# The same sounding at another station height gives the same EBWD
# (SHARPlib-27b). The inflow layer base is the ground and the EL is 4000 m
# AGL, so the layer is 0 to 2000 m AGL, where u goes from 0 to 12. Before
# the fix, u was 12, 3, 2.531 and 5.1393 for these shifts.
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


# A layer that converts to MISSING (it extends past the profile, or its end
# has no valid data beyond it) gives MISSING winds instead of raising or
# returning garbage (SHARPlib-mut). Each "was" comment is the output measured
# before that change.
def test_wind_params_missing_layer():
    pres = np.array([100000, 95000, 90000, 85000, 80000], dtype="float32")
    hght = np.array([0, 500, 1000, 1500, 2000], dtype="float32")
    uwin = np.array([0, 5, 10, 15, 20], dtype="float32")
    vwin = np.array([0, 2, 4, 6, 8], dtype="float32")

    # the inflow layer top has no valid height above it: was ValueError
    hght_mis = hght.copy()
    hght_mis[-1] = constants.MISSING
    ebwd = params.effective_bulk_wind_difference(
        pres, hght_mis, uwin, vwin, layer.PressureLayer(100000, 80000), 95000)
    assert (ebwd.u == constants.MISSING and ebwd.v == constants.MISSING)

    # a mean wind layer past the top of the profile: was (-9996.2, -10006.0)
    storm_mtn = params.storm_motion_bunkers(
        pres, hght, uwin, vwin, layer.HeightLayer(0, 3000),
        layer.HeightLayer(0, 2000))
    assert (storm_mtn.u == constants.MISSING and
            storm_mtn.v == constants.MISSING)


def test_bunkers_motion_effective_fallback():
    # An inflow layer below the profile falls back to the non-parcel method
    # with 0-6 km layers: was ValueError
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
    assert (ehi == pytest.approx(4.38889, abs=1e-5))


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
