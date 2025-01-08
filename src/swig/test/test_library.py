import pandas as pd
import numpy as np

# fun little issue is that numpy and pandas
# needs to be imported first, or bad things
# happen
import nwsspc.sharp.calc.constants as constants
import nwsspc.sharp.calc.interp as interp
import nwsspc.sharp.calc.thermo as thermo
import nwsspc.sharp.calc.parcel as parcel
import nwsspc.sharp.calc.params as params
import nwsspc.sharp.calc.winds as winds
import nwsspc.sharp.calc.layer as layer


def load_snd1(filename):
    names = ["pres", "hght", "tmpk", "dwpk", "wdir", "wspd", "omeg"]
    # snd_df = pd.read_csv(filename, delimiter=",", comment="%", names=names, skiprows=7)
    snd_df = pd.read_parquet(filename)
    snd_df = snd_df[snd_df["vwin"].notna()]
    snd_df = snd_df[snd_df["tmpc"].notna()]
    snd_df = snd_df[snd_df["relh"].notna()]
    snd_df = snd_df[snd_df["pres"] >= 50.0]
    print(snd_df.columns)

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

    # for idx in range(len(uwin)):
    #     comp = winds.vector_to_components(float(wspd[idx]), float(wdir[idx]))
    #     if (comp.u != constants.MISSING):
    #         uwin[idx] = comp.u * 0.514444 ## convert to m/s
    #         vwin[idx] = comp.v * 0.514444 ## convert to m/s
    #     else:
    #         uwin[idx] = comp.u
    #         vwin[idx] = comp.v

    return {"pres": pres, "hght": hght, "tmpk": tmpk, "mixr": mixr, "vtmp": vtmp, "dwpk": dwpk, "uwin": uwin, "vwin": vwin}


def load_snd(filename):
    names = ["pres", "hght", "tmpk", "dwpk", "wdir", "wspd", "omeg"]
    snd_df = pd.read_csv(filename, delimiter=",",
                         comment="%", names=names, skiprows=7)

    pres = snd_df["pres"].to_numpy().astype('float32')*100.0
    hght = snd_df["hght"].to_numpy().astype('float32')
    tmpk = snd_df["tmpk"].to_numpy().astype('float32')+273.15
    dwpk = snd_df["dwpk"].to_numpy().astype('float32')+273.15
    wdir = snd_df["wdir"].to_numpy().astype('float32')
    wspd = snd_df["wspd"].to_numpy().astype('float32')

    # TO-DO - need a better interface to the API for doing this
    uwin = np.empty(wspd.shape, dtype="float32")
    vwin = np.empty(wspd.shape, dtype="float32")

    for idx in range(len(uwin)):
        comp = winds.vector_to_components(float(wspd[idx]), float(wdir[idx]))
        if (comp.u != constants.MISSING):
            uwin[idx] = comp.u * 0.514444  # convert to m/s
            vwin[idx] = comp.v * 0.514444  # convert to m/s
        else:
            uwin[idx] = comp.u
            vwin[idx] = comp.v

    return {"pres": pres, "hght": hght, "tmpk": tmpk, "dwpk": dwpk, "wdir": wdir, "wspd": wspd, "uwin": uwin, "vwin": vwin}


def test_interp(snd):
    print("====================")
    print("Testing the Interpolation bindings...")

    hght_1000hPa = interp.interp_pressure(100000.0, snd["pres"], snd["hght"])
    hght_500hPa = interp.interp_pressure(50000.0, snd["pres"], snd["hght"])
    print("1000 hPa Height: ", hght_1000hPa)
    print("500 hPa Height: ", hght_500hPa)

    tmpk_1000m = interp.interp_height(1000.0, snd["hght"], snd["tmpk"])
    tmpk_3000m = interp.interp_height(3000.0, snd["hght"], snd["tmpk"])
    print("1 km Temperature: ", tmpk_1000m)
    print("3 km Temperature: ", tmpk_3000m)

    hght_vals = np.arange(0, 3000.0, 100.0).astype('float32')
    pres_vals = np.arange(snd["pres"][0], snd["pres"]
                          [0] - 100.0*100.0, -10.0*100.0).astype('float32')

    tmpk_vals = interp.interp_height(hght_vals, snd["hght"], snd["tmpk"])
    dwpk_vals = interp.interp_pressure(pres_vals, snd["pres"], snd["dwpk"])
    print(tmpk_vals)
    print(dwpk_vals)
    print("====================")


def test_layer(snd):

    print("====================")
    print("Testing layer bindings...")
    pres_layer = layer.PressureLayer(100000.0, 50000.0)
    hght_layer = layer.HeightLayer(1000.0, 3000.0)

    idx1 = layer.get_layer_index(pres_layer, snd["pres"])
    idx2 = layer.get_layer_index(hght_layer, snd["hght"])

    print("Pressure Layer: 1000.0 hPa to 500.0 hPa")
    print("Pressure Layer kbot: ", idx1.kbot, "ktop: ", idx1.ktop)

    print("Height Layer: 1000.0 m to 3000.0 m")
    print("Height Layer kbot: ", idx2.kbot, "ktop: ", idx2.ktop)

    hlayer_from_player = layer.pressure_layer_to_height(
        pres_layer, snd["pres"], snd["hght"])
    player_from_hlayer = layer.height_layer_to_pressure(
        hght_layer, snd["pres"], snd["hght"])

    print(hlayer_from_player.bottom, hlayer_from_player.top)
    print(player_from_hlayer.bottom, player_from_hlayer.top)

    hlayer_min = layer.layer_min(hght_layer, snd["hght"], snd["dwpk"])
    player_min = layer.layer_min(pres_layer, snd["pres"], snd["tmpk"])
    print(hlayer_min)
    print(player_min)

    hlayer_max = layer.layer_max(hght_layer, snd["hght"], snd["dwpk"])
    player_max = layer.layer_max(pres_layer, snd["pres"], snd["tmpk"])
    print(hlayer_max)
    print(player_max)

    hlayer_mean = layer.layer_mean(
        hght_layer, snd["hght"], snd["pres"], snd["dwpk"])
    player_mean = layer.layer_mean(pres_layer, snd["pres"], snd["tmpk"])
    print(hlayer_mean)
    print(player_mean)

    print("====================")


def test_winds(snd):

    print("====================")
    print("Testing kinematic bindings...")
    pres_layer = layer.PressureLayer(100000.0, 50000.0)
    hght_layer = layer.HeightLayer(0.0, 3000.0)

    comp1 = winds.mean_wind(
        pres_layer, snd["pres"], snd["uwin"], snd["vwin"], True)
    comp2 = winds.mean_wind(
        pres_layer, snd["pres"], snd["uwin"], snd["vwin"], False)
    print("Pressure Weighted Mean Wind: u = ", comp1.u, "v = ", comp1.v)
    print("Non-Pressure Weighted Mean Wind: u = ", comp2.u, "v = ", comp2.v)

    shear = winds.wind_shear(pres_layer, snd["pres"], snd["uwin"], snd["vwin"])
    print("1000 hPa - 500 hPa wind shear: u = ", shear.u, "v = ", shear.v)

    strm_motnv = winds.WindComponents()
    strm_motnv.u = 0.0
    strm_motnv.v = 0.0
    helicity = winds.helicity(hght_layer, strm_motnv,
                              snd["hght"], snd["uwin"], snd["vwin"])
    print("0-3 km AGL Storm Relative Helicity: ", helicity)
    print("====================")


def test_thermo(snd):
    print("====================")
    print("Testing thermodynamic bindings...")

    hlayer = layer.HeightLayer(0, 3000.0)
    lapse_rate = thermo.lapse_rate(hlayer, snd["hght"], snd["tmpk"])
    print("0-3km AGL Lapse Rate: ", lapse_rate)

    print("Vapor Pressure (Pa)")
    vappres = thermo.vapor_pressure(snd["pres"], snd["tmpk"])
    print(thermo.vapor_pressure(100000.0, 300.0))
    print(vappres)

    print("Vapor Pressure (Pa) with respect to ice")
    vappres_ice = thermo.vapor_pressure_ice(snd["pres"], snd["tmpk"])
    print(vappres_ice)

    print("LCL Temperature (K)")
    print(thermo.lcl_temperature(300.0, 300.0))
    lcl_t = thermo.lcl_temperature(snd["tmpk"], snd["dwpk"])
    print(lcl_t)

    print("Potential Temperature (K)")
    print(thermo.theta(100000.0, 300.0))
    theta = thermo.theta(snd["pres"], snd["tmpk"])
    print(theta)

    print("Theta Level (Pa)")
    print(thermo.theta_level(300.0, 300.0))
    plev = thermo.theta_level(theta, snd["tmpk"])
    print(plev)

    print("Mixing Raio (kg/kg)")
    print(thermo.mixratio(100000.0, 300.0))
    mixr = thermo.mixratio(snd["pres"], snd["dwpk"])
    print(mixr)

    print("Mixing Raio (ice) (kg/kg)")
    print(thermo.mixratio_ice(100000.0, 300.0))
    mixr_i = thermo.mixratio_ice(snd["pres"], snd["dwpk"])
    print(mixr_i)

    print("Temperature At Mixratio (K)")
    print(thermo.temperature_at_mixratio(0.010, 90000.0))
    dwpk = thermo.temperature_at_mixratio(mixr, snd["pres"])
    print(dwpk)

    print("Specific Humidity")
    print(thermo.specific_humidity(0.010))
    spcf = thermo.specific_humidity(mixr)
    print(spcf)

    print("Virtual Temperature (simple) (K)")
    print(thermo.virtual_temperature(300.0, 0.010))
    vtmp = thermo.virtual_temperature(snd["tmpk"], mixr)
    print(vtmp)

    print("Virtual Temperature (full) (K)")
    print(thermo.virtual_temperature(300.0, 0.010, 0.0, 0.0))
    ql = np.zeros(vtmp.shape, dtype='float32')
    qi = np.zeros(vtmp.shape, dtype='float32')
    vtmp2 = thermo.virtual_temperature(snd["tmpk"], mixr, ql, qi)
    print(vtmp2)

    # print("Wetbulb (K)")
    # print(thermo.wetbulb(100000.0, 300.0, 298.0))
    # wetb = thermo.wetbulb(snd["pres"], snd["tmpk"], snd["dwpk"])
    # print(wetb)
    #
    # print("Theta Wetbulb (K)")
    # print(thermo.theta_wetbulb(100000.0, 300.0, 298.0))
    # twetb = thermo.theta_wetbulb(snd["pres"], snd["tmpk"], snd["dwpk"])
    # print(twetb)

    print("Theta E (K)")
    print(thermo.thetae(100000.0, 300.0, 298.0))
    thte = thermo.thetae(snd["pres"], snd["tmpk"], snd["dwpk"])
    print(thte)

    print("Buoyancy (m/s^2)")
    print(thermo.buoyancy(300.0, 300.0))
    buoy = thermo.buoyancy(snd["tmpk"]+2, snd["tmpk"])
    print(buoy)

    print("Moist Static Energy")
    print(thermo.moist_static_energy(0.0, 300.0, 0.010))
    mse = thermo.moist_static_energy(
        snd["hght"] - snd["hght"][0], snd["tmpk"], spcf)
    print(mse)

    print("Max Lapse Rate (Height AGL)")
    lyr_search = layer.HeightLayer(2000.0, 6000.0, 250.0)
    lyr_max = layer.HeightLayer(constants.MISSING, constants.MISSING)
    depth = 2000.0
    max_lr = thermo.lapse_rate_max(
        lyr_search, depth, snd["hght"], snd["tmpk"], lyr_max)
    print(max_lr)
    print(lyr_max.bottom, lyr_max.top)

    print("Max Lapse Rate (Pres)")
    lyr_search = layer.PressureLayer(80000.0, 50000.0, -500.0)
    lyr_max = layer.PressureLayer(constants.MISSING, constants.MISSING)
    depth = 20000.0
    max_lr = thermo.lapse_rate_max(
        lyr_search, depth, snd["pres"], snd["hght"], snd["tmpk"], lyr_max)
    print(max_lr)
    print(lyr_max.bottom, lyr_max.top)

    print("====================")


def test_parcel(snd):
    print("====================")
    print("Testing parcel bindings...")

    vtmp = np.zeros(snd["pres"].shape, dtype='float32')
    mixr = np.zeros(snd["pres"].shape, dtype='float32')
    mixr = thermo.mixratio(snd["pres"], snd["dwpk"])
    theta = thermo.theta(snd["pres"], snd["tmpk"])

    wobf = parcel.lifter_wobus()
    cm1 = parcel.lifter_cm1()
    cm1.pressure_incr = 1000.0
    cm1.converge = 0.01

    idx = 0
    for p, t, m in zip(snd["pres"], snd["tmpk"], mixr):
        vtmp[idx] = thermo.virtual_temperature(float(t), float(m))
        idx = idx + 1

    # define a mix layer
    mix_lyr_pr = layer.PressureLayer(
        float(snd["pres"][0]), float(snd["pres"][0] - 10000.0))
    mix_lyr_ht = layer.pressure_layer_to_height(
        mix_lyr_pr, snd["pres"], snd["hght"])

    # Create a surface-based parcel
    sfc_pcl = parcel.Parcel.surface_parcel(
        float(snd["pres"][0]), float(snd["tmpk"][0]), float(snd["dwpk"][0]))
    print("Surface-based parcel attributes")
    print(sfc_pcl.pres, sfc_pcl.tmpk, sfc_pcl.dwpk)

    # create a mixed-layer parcel
    ml_pcl1 = parcel.Parcel.mixed_layer_parcel(
        mix_lyr_pr, snd["pres"], snd["hght"], theta, mixr)
    print("Mixed-layer parcel attributes (PressureLayer)")
    print(ml_pcl1.pres, ml_pcl1.tmpk, ml_pcl1.dwpk)

    # test that the HeightLayer function works too
    ml_pcl2 = parcel.Parcel.mixed_layer_parcel(
        mix_lyr_ht, snd["pres"], snd["hght"], theta, mixr)
    print("Mixed-layer parcel attributes (HeightLayer)")
    print(ml_pcl2.pres, ml_pcl2.tmpk, ml_pcl2.dwpk)

    # Compute buoyancy from a surface-based parcel
    sfc_vtmp1 = sfc_pcl.lift_parcel(wobf, snd["pres"])
    sfc_buoy1 = thermo.buoyancy(sfc_vtmp1, vtmp)
    print("Surface-based parcel buoyancy (Wobus)")
    print(sfc_buoy1, sfc_buoy1.min(), sfc_buoy1.max())
    sfc_vtmp2 = sfc_pcl.lift_parcel(cm1, snd["pres"])
    sfc_buoy2 = thermo.buoyancy(sfc_vtmp2, vtmp)
    print("Surface-based parcel buoyancy (CM1)")
    print(sfc_buoy2, sfc_buoy2.min(), sfc_buoy2.max())

    # Compute buoyancy from a mixed-layer parcel
    ml_vtmp1 = ml_pcl1.lift_parcel(wobf, snd["pres"])
    ml_buoy1 = thermo.buoyancy(ml_vtmp1, vtmp)
    print("Mixed-layer parcel buoyancy (Wobus)")
    print(ml_buoy1, ml_buoy1.min(), ml_buoy1.max())
    ml_vtmp2 = ml_pcl1.lift_parcel(cm1, snd["pres"])
    ml_buoy2 = thermo.buoyancy(ml_vtmp2, vtmp)
    print("Mixed-layer parcel buoyancy (CM1)")
    print(ml_buoy2, ml_buoy2.min(), ml_buoy2.max())

    # Find the LFC and EL for a surface-based parcel
    sfc_lfc1, sfc_el1 = sfc_pcl.find_lfc_el(
        snd["pres"], snd["hght"], sfc_buoy1)
    print("Surface-based parcel LFC and EL pressure (Wobus)")
    print(sfc_lfc1, sfc_el1, sfc_pcl.lfc_pressure, sfc_pcl.eql_pressure)
    sfc_lfc2, sfc_el2 = sfc_pcl.find_lfc_el(
        snd["pres"], snd["hght"], sfc_buoy2)
    print("Surface-based parcel LFC and EL pressure (CM1)")
    print(sfc_lfc2, sfc_el2, sfc_pcl.lfc_pressure, sfc_pcl.eql_pressure)

    # Find the LFC and EL for a mixed-layer parcel
    ml_lfc1, ml_el1 = ml_pcl1.find_lfc_el(snd["pres"], snd["hght"], ml_buoy1)
    print("Mixed-layer parcel LFC and EL pressure (Wobus)")
    print(ml_lfc1, ml_el1, ml_pcl1.lfc_pressure, ml_pcl1.eql_pressure)
    ml_lfc2, ml_el2 = ml_pcl1.find_lfc_el(snd["pres"], snd["hght"], ml_buoy2)
    print("Mixed-layer parcel LFC and EL pressure (CM1)")
    print(ml_lfc2, ml_el2, ml_pcl1.lfc_pressure, ml_pcl1.eql_pressure)

    # Compute CAPE and CINH from a surface-based parcel
    sfc_cape1, sfc_cinh1 = sfc_pcl.cape_cinh(
        snd["pres"], snd["hght"], sfc_buoy1)
    print("Surface-based parcel CAPE and CINH (Wobus)")
    print(sfc_cape1, sfc_cinh1, sfc_pcl.cape, sfc_pcl.cinh)
    sfc_cape2, sfc_cinh2 = sfc_pcl.cape_cinh(
        snd["pres"], snd["hght"], sfc_buoy2)
    print("Surface-based parcel CAPE and CINH (CM1)")
    print(sfc_cape2, sfc_cinh2, sfc_pcl.cape, sfc_pcl.cinh)

    # Compute CAPE and CINH from a mixed-layer parcel
    ml_cape1, ml_cinh1 = ml_pcl1.cape_cinh(snd["pres"], snd["hght"], ml_buoy1)
    print("Mixed-layer parcel CAPE and CINH (Wobus)")
    print(ml_cape1, ml_cinh1, ml_pcl1.cape, ml_pcl1.cinh)
    ml_cape2, ml_cinh2 = ml_pcl1.cape_cinh(snd["pres"], snd["hght"], ml_buoy2)
    print("Mixed-layer parcel CAPE and CINH (CM1)")
    print(ml_cape2, ml_cinh2, ml_pcl1.cape, ml_pcl1.cinh)

    # find and compute the most-unstable parcel
    search_lyr = layer.PressureLayer(float(snd["pres"][0]), float(70000.0))
    mu_pcl1 = parcel.Parcel.most_unstable_parcel(
        wobf,
        search_lyr,
        snd["pres"],
        snd["hght"],
        snd["tmpk"],
        vtmp,
        snd["dwpk"]
    )
    print("Most-unstable parcel attributes (Wobus)")
    print("PCL PRES, TMPK, DWPK: ", mu_pcl1.pres, mu_pcl1.tmpk, mu_pcl1.dwpk)
    print("LFC, EL: ", mu_pcl1.lfc_pressure, mu_pcl1.eql_pressure)
    print("CAPE, CINH: ", mu_pcl1.cape, mu_pcl1.cinh)

    mu_pcl2 = parcel.Parcel.most_unstable_parcel(
        cm1,
        search_lyr,
        snd["pres"],
        snd["hght"],
        snd["tmpk"],
        vtmp,
        snd["dwpk"]
    )
    print("Most-unstable parcel attributes (CM1)")
    print("PCL PRES, TMPK, DWPK: ", mu_pcl2.pres, mu_pcl2.tmpk, mu_pcl2.dwpk)
    print("LFC, EL: ", mu_pcl2.lfc_pressure, mu_pcl2.eql_pressure)
    print("CAPE, CINH: ", mu_pcl2.cape, mu_pcl2.cinh)

    mu_pcl3 = parcel.Parcel()
    eil = params.effective_inflow_layer(wobf, snd["pres"], snd["hght"],
                                        snd["tmpk"], snd["dwpk"],
                                        vtmp, mu_pcl3)
    print("Effective Inflow Layer calculations and parcels (Wobus)")
    print("EIL bottom, top: ", eil.bottom, eil.top)
    print("PCL PRES, TMPK, DWPK: ", mu_pcl3.pres, mu_pcl3.tmpk, mu_pcl3.dwpk)
    print("LFC, EL: ", mu_pcl3.lfc_pressure, mu_pcl3.eql_pressure)
    print("CAPE, CINH: ", mu_pcl3.cape, mu_pcl3.cinh)

    mu_pcl4 = parcel.Parcel()
    eil = params.effective_inflow_layer(cm1, snd["pres"], snd["hght"],
                                        snd["tmpk"], snd["dwpk"],
                                        vtmp, mu_pcl4)
    print("Effective Inflow Layer calculations and parcels (CM1)")
    print("EIL bottom, top: ", eil.bottom, eil.top)
    print("PCL PRES, TMPK, DWPK: ", mu_pcl4.pres, mu_pcl4.tmpk, mu_pcl4.dwpk)
    print("LFC, EL: ", mu_pcl4.lfc_pressure, mu_pcl4.eql_pressure)
    print("CAPE, CINH: ", mu_pcl4.cape, mu_pcl4.cinh)

    mw_lyr = layer.HeightLayer(0, 6000.0)
    shr_lyr = layer.HeightLayer(0, 6000.0)
    strm_right_np = params.storm_motion_bunkers(snd["pres"], snd["hght"],
                                                snd["uwin"], snd["vwin"],
                                                mw_lyr, shr_lyr)
    strm_right_eff = params.storm_motion_bunkers(snd["pres"], snd["hght"],
                                                 snd["uwin"], snd["vwin"],
                                                 eil, mu_pcl1)
    print(
        "Bunkers Right (non-parcel): u = {0}\tv = {1}"
        .format(strm_right_np.u, strm_right_np.v)
    )
    print("Bunkers Right (parcel): u = {0}\tv = {1}".format(
        strm_right_eff.u, strm_right_eff.v)
    )

    print("====================")


def main():
    # snd = load_snd("./hires-SPC.txt")
    snd = load_snd1("./ddc.parquet")
    test_interp(snd)
    test_layer(snd)
    test_parcel(snd)
    test_thermo(snd)
    test_winds(snd)


if __name__ == "__main__":
    main()
