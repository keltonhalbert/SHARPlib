import pandas as pd
import numpy as np

## fun little issue is that numpy and pandas
## needs to be imported first, or bad things
## happen
import nwsspc.sharp.calc.constants as constants
import nwsspc.sharp.calc.profile as profile
import nwsspc.sharp.calc.interp as interp
import nwsspc.sharp.calc.thermo as thermo
import nwsspc.sharp.calc.parcel as parcel
import nwsspc.sharp.calc.params as params
import nwsspc.sharp.calc.winds as winds
import nwsspc.sharp.calc.layer as layer

def load_snd(filename):
    names = ["pres", "hght", "tmpk", "dwpk", "wdir", "wspd", "omeg"]
    snd_df = pd.read_csv(filename, delimiter=",", comment="%", names=names, skiprows=7)

    pres = snd_df["pres"].to_numpy().astype('float32')*100.0
    hght = snd_df["hght"].to_numpy().astype('float32')
    tmpk = snd_df["tmpk"].to_numpy().astype('float32')+273.15
    dwpk = snd_df["dwpk"].to_numpy().astype('float32')+273.15
    wdir = snd_df["wdir"].to_numpy().astype('float32')
    wspd = snd_df["wspd"].to_numpy().astype('float32')

    ## TO-DO - need a better interface to the API for doing this
    uwin = np.empty(wspd.shape, dtype="float32")
    vwin = np.empty(wspd.shape, dtype="float32")

    for idx in range(len(uwin)):
        comp = winds.vector_to_components(float(wspd[idx]), float(wdir[idx]))
        if (comp.u != constants.MISSING):
            uwin[idx] = comp.u * 0.514444 ## convert to m/s
            vwin[idx] = comp.v * 0.514444 ## convert to m/s
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
    pres_vals = np.arange(snd["pres"][0], snd["pres"][0] - 100.0*100.0, -10.0*100.0).astype('float32')

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

    hlayer_from_player = layer.pressure_layer_to_height(pres_layer, snd["pres"], snd["hght"])
    player_from_hlayer = layer.height_layer_to_pressure(hght_layer, snd["pres"], snd["hght"])

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

    hlayer_mean = layer.layer_mean(hght_layer, snd["hght"], snd["pres"], snd["dwpk"])
    player_mean = layer.layer_mean(pres_layer, snd["pres"], snd["tmpk"])
    print(hlayer_mean)
    print(player_mean)

    print("====================")

def test_winds(snd):

    print("====================")
    print("Testing kinematic bindings...")
    prof = profile.create_profile(snd["pres"], snd["hght"], 
                                  snd["tmpk"], snd["dwpk"], 
                                  snd["wspd"], snd["wdir"], 
                                  profile.Source_Observed, False) 

    pres_layer = layer.PressureLayer(100000.0, 50000.0)
    hght_layer = layer.HeightLayer(0.0, 3000.0)

    comp1 = winds.mean_wind(pres_layer, snd["pres"], snd["uwin"], snd["vwin"])
    comp2 = winds.mean_wind_npw(pres_layer, snd["pres"], snd["uwin"], snd["vwin"])
    print("Pressure Weighted Mean Wind: u = ", comp1.u, "v = ", comp1.v)
    print("Non-Pressure Weighted Mean Wind: u = ", comp2.u, "v = ", comp2.v)

    shear = winds.wind_shear(pres_layer, snd["pres"], snd["uwin"], snd["vwin"])
    print("1000 hPa - 500 hPa wind shear: u = ", shear.u, "v = ", shear.v)

    helicity = winds.helicity(hght_layer, strm_motnv, snd["hght"], snd["uwin"], snd["vwin"])
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

    print("Wetbulb (K)")
    print(thermo.wetbulb(100000.0, 300.0, 298.0))
    wetb = thermo.wetbulb(snd["pres"], snd["tmpk"], snd["dwpk"])
    print(wetb)

    print("Theta Wetbulb (K)")
    print(thermo.theta_wetbulb(100000.0, 300.0, 298.0))
    twetb = thermo.theta_wetbulb(snd["pres"], snd["tmpk"], snd["dwpk"])
    print(twetb)

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
    mse = thermo.moist_static_energy(snd["hght"] - snd["hght"][0], snd["tmpk"], spcf)
    print(mse)


    print("====================")

def test_parcel(snd):
    print("====================")
    print("Testing parcel bindings...")

    vtmp = np.zeros(snd["pres"].shape, dtype='float32') 
    mixr = np.zeros(snd["pres"].shape, dtype='float32')
    qi = np.zeros(snd["pres"].shape, dtype='float32')
    ql = np.zeros(snd["pres"].shape, dtype='float32')

    mixr = thermo.mixratio(snd["pres"], snd["dwpk"])

    idx = 0
    for p, t, m in zip(snd["pres"], snd["tmpk"], mixr):
        vtmp[idx] = thermo.virtual_temperature(float(t), float(m))
        idx = idx + 1 

    mupcl = parcel.Parcel()
    eil = params.effective_inflow_layer(snd["pres"], snd["hght"], 
                                        snd["tmpk"], snd["dwpk"], 
                                        vtmp, mupcl)
    print(eil.bottom, eil.top)
    print(mupcl.pres, mupcl.lcl_pressure, mupcl.lfc_pressure, mupcl.eql_pressure, mupcl.cape, mupcl.cinh)

    mw_lyr = layer.HeightLayer(0, 6000.0)
    shr_lyr = layer.HeightLayer(0, 6000.0)
    strm_right_np = params.storm_motion_bunkers(snd["pres"], snd["hght"], 
                                                snd["uwin"], snd["vwin"],
                                                mw_lyr, shr_lyr)
    strm_right_eff = params.storm_motion_bunkers(snd["pres"], snd["hght"],
                                                 snd["uwin"], snd["vwin"],
                                                 eil, mupcl)

    print("Bunkers Right (non-parcel): u = {0}\tv = {1}".format(strm_right_np.u, strm_right_np.v))
    print("Bunkers Right (parcel): u = {0}\tv = {1}".format(strm_right_eff.u, strm_right_eff.v))

    print("====================")


def main():
    snd = load_snd("./hires-SPC.txt")
    test_interp(snd)
    test_layer(snd)
    test_parcel(snd)
    test_thermo(snd)
    test_winds(snd)

if __name__ == "__main__":
    main()
