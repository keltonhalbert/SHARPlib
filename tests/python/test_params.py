
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

    pres = snd_df["pres"].to_numpy().astype('float32')*100.0
    hght = snd_df["hght"].to_numpy().astype('float32')
    tmpk = snd_df["tmpc"].to_numpy().astype('float32')+273.15
    dwpk = snd_df["dwpc"].to_numpy().astype('float32')+273.15
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

    return {
        "pres": pres, "hght": hght,
        "tmpk": tmpk, "mixr": mixr,
        "relh": relh,
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

    assert (storm_mtn.u == pytest.approx(10.232333))
    assert (storm_mtn.v == pytest.approx(5.7385511))


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

    assert (storm_mtn.u == pytest.approx(9.65811))
    assert (storm_mtn.v == pytest.approx(5.558156))


def test_corfidi_vectors():
    upshear, downshear = params.mcs_motion_corfidi(
        snd_data["pres"],
        snd_data["hght"],
        snd_data["uwin"],
        snd_data["vwin"]
    )

    assert (upshear.u == pytest.approx(12.7017, abs=1e-3))
    assert (upshear.v == pytest.approx(2.99329, abs=1e-3))
    assert (downshear.u == pytest.approx(23.2054, abs=1e-3))
    assert (downshear.v == pytest.approx(16.10519, abs=1e-3))


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


def test_stp_scp_ship():
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

    # lift the parcel and get CAPE
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
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
    assert (stp == pytest.approx(0.48329, abs=1e-4))

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
    assert (ehi == pytest.approx(4.41228, abs=1e-5))


def test_precipitable_water():
    plyr = layer.PressureLayer(snd_data["pres"][0], 40000.0)
    pwat = params.precipitable_water(plyr, snd_data["pres"], snd_data["mixr"])

    assert (pwat == pytest.approx(21.11469))


def test_hgz():
    hgz = params.hail_growth_layer(snd_data["pres"], snd_data["tmpk"])
    assert (hgz.bottom == 51787)
    assert (hgz.top == 35861)


def test_dgz():
    dgz = params.dendritic_layer(snd_data["pres"], snd_data["tmpk"])
    assert (dgz.bottom == 49598)
    assert (dgz.top == 46032)


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


def test_rainfall_velocity():
    v = params.rainfall_velocity(
        1.0,
        1.225,
        0.01
    )
    print(v)


def test_rainfall_efficiency():
    mix_lyr = layer.PressureLayer(
        snd_data["pres"][0], snd_data["pres"][0] - 10000.0)
    theta = thermo.theta(snd_data["pres"], snd_data["tmpk"])
    pcl = parcel.Parcel.mixed_layer_parcel(
        mix_lyr,
        snd_data["pres"],
        theta,
        snd_data["mixr"]
    )
    lifter = parcel.lifter_cm1()
    lifter.ma_type = thermo.adiabat.pseudo_liq
    vtmpk = pcl.lift_parcel(lifter, snd_data["pres"])
    buoy = thermo.buoyancy(vtmpk, snd_data["vtmp"])
    cape, cinh = pcl.cape_cinh(snd_data["pres"], snd_data["hght"], buoy)
    lcl_t = interp.interp_pressure(
        pcl.lcl_pressure, snd_data["pres"], snd_data["tmpk"])
    rl = thermo.mixratio(pcl.lcl_pressure, lcl_t)

    precip_eff = params.rainfall_efficiency(
        snd_data["pres"],
        snd_data["hght"],
        snd_data["tmpk"],
        snd_data["mixr"],
        pcl,
        rl
    )

    print(precip_eff)
