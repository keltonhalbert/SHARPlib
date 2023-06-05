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

def load_sounding(filename):
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


def test_interp(sounding):
    print("====================")
    print("Testing the Interpolation bindings...")

    hght_1000hPa = interp.interp_pressure(100000.0, sounding["pres"], sounding["hght"])
    hght_500hPa = interp.interp_pressure(50000.0, sounding["pres"], sounding["hght"])
    print("1000 hPa Height: ", hght_1000hPa)
    print("500 hPa Height: ", hght_500hPa)

    tmpk_1000m = interp.interp_height(1000.0, sounding["hght"], sounding["tmpk"])
    tmpk_3000m = interp.interp_height(3000.0, sounding["hght"], sounding["tmpk"])
    print("1 km Temperature: ", tmpk_1000m)
    print("3 km Temperature: ", tmpk_3000m)

    hght_vals = np.arange(0, 3000.0, 100.0).astype('float32')
    pres_vals = np.arange(sounding["pres"][0], sounding["pres"][0] - 100.0*100.0, -10.0*100.0).astype('float32')

    tmpk_vals = interp.interp_height(hght_vals, sounding["hght"], sounding["tmpk"])
    dwpk_vals = interp.interp_pressure(pres_vals, sounding["pres"], sounding["dwpk"])
    print(tmpk_vals)
    print(dwpk_vals)
    print("====================")

def test_layer(sounding):

    print("====================")
    print("Testing layer bindings...")
    pres_layer = layer.PressureLayer(100000.0, 50000.0)
    hght_layer = layer.HeightLayer(1000.0, 3000.0)

    idx1 = layer.get_layer_index(pres_layer, sounding["pres"])
    idx2 = layer.get_layer_index(hght_layer, sounding["hght"])

    print("Pressure Layer: 1000.0 hPa to 500.0 hPa")
    print("Pressure Layer kbot: ", idx1.kbot, "ktop: ", idx1.ktop)

    print("Height Layer: 1000.0 m to 3000.0 m")
    print("Height Layer kbot: ", idx2.kbot, "ktop: ", idx2.ktop)

    hlayer_from_player = layer.pressure_layer_to_height(pres_layer, sounding["pres"], sounding["hght"])
    player_from_hlayer = layer.height_layer_to_pressure(hght_layer, sounding["pres"], sounding["hght"])

    print(hlayer_from_player.bottom, hlayer_from_player.top)
    print(player_from_hlayer.bottom, player_from_hlayer.top)

    hlayer_min = layer.layer_min(hght_layer, sounding["hght"], sounding["dwpk"])
    player_min = layer.layer_min(pres_layer, sounding["pres"], sounding["tmpk"])
    print(hlayer_min)
    print(player_min)

    hlayer_max = layer.layer_max(hght_layer, sounding["hght"], sounding["dwpk"])
    player_max = layer.layer_max(pres_layer, sounding["pres"], sounding["tmpk"])
    print(hlayer_max)
    print(player_max)

    hlayer_mean = layer.layer_mean(hght_layer, sounding["hght"], sounding["pres"], sounding["dwpk"])
    player_mean = layer.layer_mean(pres_layer, sounding["pres"], sounding["tmpk"])
    print(hlayer_mean)
    print(player_mean)

    print("====================")

def test_winds(sounding):

    print("====================")
    print("Testing kinematic bindings...")
    prof = profile.create_profile(sounding["pres"], sounding["hght"], 
                                  sounding["tmpk"], sounding["dwpk"], 
                                  sounding["wspd"], sounding["wdir"], 
                                  profile.Source_Observed, False) 

    pres_layer = layer.PressureLayer(100000.0, 50000.0)
    hght_layer = layer.HeightLayer(0.0, 3000.0)
    strm_motnv = params.storm_motion_bunkers(prof, False)

    comp1 = winds.mean_wind(pres_layer, sounding["pres"], sounding["uwin"], sounding["vwin"])
    comp2 = winds.mean_wind_npw(pres_layer, sounding["pres"], sounding["uwin"], sounding["vwin"])
    print("Pressure Weighted Mean Wind: u = ", comp1.u, "v = ", comp1.v)
    print("Non-Pressure Weighted Mean Wind: u = ", comp2.u, "v = ", comp2.v)

    shear = winds.wind_shear(pres_layer, sounding["pres"], sounding["uwin"], sounding["vwin"])
    print("1000 hPa - 500 hPa wind shear: u = ", shear.u, "v = ", shear.v)

    helicity = winds.helicity(hght_layer, strm_motnv, sounding["hght"], sounding["uwin"], sounding["vwin"])
    print("0-3 km AGL Storm Relative Helicity: ", helicity)
    print("====================")

def test_thermo(sounding):
    print("====================")
    print("Testing thermodynamic bindings...")

    hlayer = layer.HeightLayer(0, 3000.0)
    lapse_rate = thermo.lapse_rate(hlayer, sounding["hght"], sounding["tmpk"])
    print("0-3km AGL Lapse Rate: ", lapse_rate)
    print("====================")

def test_parcel(sounding):
    print("====================")
    print("Testing parcel bindings...")

    vtmp = np.zeros(sounding["pres"].shape, dtype='float32') 
    mixr = np.zeros(sounding["pres"].shape, dtype='float32')
    qi = np.zeros(sounding["pres"].shape, dtype='float32')
    ql = np.zeros(sounding["pres"].shape, dtype='float32')

    thermo.mixratio(sounding["pres"], sounding["dwpk"], mixr)

    idx = 0
    for p, t, m in zip(sounding["pres"], sounding["tmpk"], mixr):
        vtmp[idx] = thermo.virtual_temperature(float(t), float(m))
        idx = idx + 1 

    mupcl = parcel.Parcel()
    eil = params.effective_inflow_layer(sounding["pres"], sounding["hght"], 
                                        sounding["tmpk"], sounding["dwpk"], 
                                        vtmp, mupcl)
    print(eil.bottom, eil.top)
    print(mupcl.pres, mupcl.lcl_pressure, mupcl.lfc_pressure, mupcl.eql_pressure, mupcl.cape, mupcl.cinh)

    print("====================")


def main():
    sounding = load_sounding("./hires-SPC.txt")
    test_interp(sounding)
    test_layer(sounding)
    test_parcel(sounding)
    test_winds(sounding)
    test_thermo(sounding)

if __name__ == "__main__":
    main()
