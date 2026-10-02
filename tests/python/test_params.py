
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
        "pres": pres, "hght": hght,
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


def test_effective_bulk_wind():
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
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

    ebwd_cmp = params.effective_bulk_wind_difference(
        snd_data["pres"],
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"],
        eil,
        mupcl.eql_pressure
    )

    ebwd = winds.vector_magnitude(ebwd_cmp.u, ebwd_cmp.v)
    assert (ebwd_cmp.u == pytest.approx(14.6, abs=1e-3))
    assert (ebwd_cmp.v == pytest.approx(13.321, abs=1e-3))
    assert (ebwd == pytest.approx(19.764, abs=1e-3))


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


@pytest.mark.parametrize("field, mask", [
    ("theta", snd_data["pres"] < 75000.0),
    ("uwin", snd_data["pres"] > 80000.0),
])
def test_pft_missing(field, mask):
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    pres = snd_data["pres"]
    mix_layer = layer.PressureLayer(pres[0], pres[0] - 10000.0)
    data = {"uwin": snd_data["uwin"], "theta": snd_data["theta"]}
    data[field] = np.where(mask, constants.MISSING,
                           data[field]).astype("float32")
    pcl = parcel.Parcel()
    pft = params.pyrocumulonimbus_firepower_threshold(
        lifter, mix_layer, pres, snd_data["hght"], snd_data["tmpk"],
        snd_data["mixr"], snd_data["vtmp"], data["uwin"], snd_data["vwin"],
        data["theta"], pcl=pcl)
    assert (pft == constants.MISSING)
    assert (pcl.pres == constants.MISSING)


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
