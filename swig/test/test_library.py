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
import nwsspc.sharp.calc.utils as utils

def load_sounding(filename):
    names = ["pres", "hght", "tmpc", "dwpc", "wdir", "wspd"]
    snd_df = pd.read_csv(filename, delimiter=",", comment="%", names=names, skiprows=5)

    pres = snd_df["pres"].to_numpy().astype('float32')
    hght = snd_df["hght"].to_numpy().astype('float32')
    tmpc = snd_df["tmpc"].to_numpy().astype('float32')
    dwpc = snd_df["dwpc"].to_numpy().astype('float32')
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
    
    return {"pres": pres, "hght": hght, "tmpc": tmpc, "dwpc": dwpc, "wdir": wdir, "wspd": wspd, "uwin": uwin, "vwin": vwin}


def test_interp(sounding):
    print("====================")
    print("Testing the Interpolation bindings...")

    hght_1000hPa = interp.interp_pressure(1000.0, sounding["pres"], sounding["hght"])
    hght_500hPa = interp.interp_pressure(500.0, sounding["pres"], sounding["hght"])
    print("1000 hPa Height: ", hght_1000hPa)
    print("500 hPa Height: ", hght_500hPa)

    tmpc_1000m = interp.interp_height(1000.0, sounding["hght"], sounding["tmpc"])
    tmpc_3000m = interp.interp_height(3000.0, sounding["hght"], sounding["tmpc"])
    print("1 km Temperature: ", tmpc_1000m)
    print("3 km Temperature: ", tmpc_3000m)
    print("====================")

def test_utils(sounding):

    print("====================")
    print("Testing utility bindings...")
    pres_layer = utils.PressureLayer(1000.0, 500.0)
    hght_layer = utils.HeightLayer(1000.0, 3000.0)

    idx1 = utils.get_layer_index(pres_layer, sounding["pres"])
    idx2 = utils.get_layer_index(hght_layer, sounding["hght"])

    print("Pressure Layer: 1000.0 hPa to 500.0 hPa")
    print("Pressure Layer kbot: ", idx1.kbot, "ktop: ", idx1.ktop)

    print("Height Layer: 1000.0 m to 3000.0 m")
    print("Height Layer kbot: ", idx2.kbot, "ktop: ", idx2.ktop)
    print("====================")

def test_winds(sounding):

    print("====================")
    print("Testing kinematic bindings...")
    prof = profile.create_profile(sounding["pres"], sounding["hght"], 
                                  sounding["tmpc"], sounding["dwpc"], 
                                  sounding["wspd"], sounding["wdir"], 
                                  profile.Source_Observed, False) 

    pres_layer = utils.PressureLayer(1000.0, 500.0)
    hght_layer = utils.HeightLayer(0.0, 3000.0)
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

    layer = utils.HeightLayer(0, 3000.0)
    lapse_rate = thermo.lapse_rate(layer, sounding["hght"], sounding["tmpc"])
    print("0-3km AGL Lapse Rate: ", lapse_rate)
    print("====================")

def test_parcel(sounding):
    print("====================")
    print("Testing parcel bindings...")

    ## example of processing a single parcel
    prof = profile.create_profile(sounding["pres"], sounding["hght"], 
                                  sounding["tmpc"], sounding["dwpc"], 
                                  sounding["wspd"], sounding["wdir"], 
                                  profile.Source_Observed, False) 
    pcl = parcel.Parcel()
    parcel.define_parcel(prof, pcl, parcel.LPL_MU)
    parcel.parcel_wobf(prof, pcl)
    print("Most Unstable Parcel: cape = ", pcl.cape, "cinh = ", pcl.cinh)


    ## example of processing an array of parcels
    prof_arr = np.array([prof, prof])
    pcl_arr = np.array([parcel.Parcel(), parcel.Parcel()])

    ## define the LPL of each parcel
    for pc, pf in zip(pcl_arr, prof_arr):
        parcel.define_parcel(pf, pc, parcel.LPL_MU)

    ## turn the parcel lifting routine into 
    ## something vectorizable by numpy
    cape_func = np.vectorize(parcel.parcel_wobf)

    ## call the parcel routines on each profile/parcel combo
    cape_func(prof_arr, pcl_arr)
    print("Vectorized Most Unstable CAPE: cape = ", pcl_arr[0].cape, ",", pcl_arr[1].cape, 
                                          "cinh = ", pcl_arr[0].cinh, ",", pcl_arr[1].cinh)
    print("====================")


def main():
    sounding = load_sounding("./hires-SPC.txt")
    test_interp(sounding)
    test_utils(sounding)
    test_winds(sounding)
    test_thermo(sounding)
    test_parcel(sounding)

if __name__ == "__main__":
    main()
