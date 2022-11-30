import numpy as np
import nwsspc.sharp.calc.winds as winds
import nwsspc.sharp.calc.interp as interp
import nwsspc.sharp.calc.utils as utils


def test_interp():

    pres = np.arange(50.0, 1010.0, 50.0, dtype='float32')[::-1]
    hght = np.linspace(0, 15000, pres.shape[0]).astype('float32')
    data = np.linspace(1, 100, pres.shape[0]).astype('float32')
    

    test_1000hPa = interp.interp_pressure(1000.0, pres, data)
    test_500hPa = interp.interp_pressure(500.0, pres, data)
    print(test_1000hPa)
    print(test_500hPa)
    print(data[pres==500.0])

    test_100m = interp.interp_height(100.0, hght, data)
    test_1000m = interp.interp_height(1000.0, hght, data)
    print(test_100m)
    print(test_1000m)

def test_utils():
    pres_layer = utils.PressureLayer(1000.0, 500.0)
    hght_layer = utils.HeightLayer(100.0, 5000.0)

    pres = np.arange(50.0, 1010.0, 50.0, dtype='float32')[::-1]
    hght = np.linspace(0, 15000, pres.shape[0]).astype('float32')
    data = np.linspace(1, 100, pres.shape[0]).astype('float32')

    idx1 = utils.get_layer_index(pres_layer, pres)
    idx2 = utils.get_layer_index(hght_layer, hght)

    print(idx1.kbot, idx1.ktop)
    print(idx2.kbot, idx2.ktop)

def test_winds():
    pres = np.arange(50.0, 1010.0, 50.0, dtype='float32')[::-1]
    hght = np.linspace(0, 15000, pres.shape[0]).astype('float32')
    uwin = np.linspace(1, 100, pres.shape[0]).astype('float32')
    vwin = np.linspace(1, 100, pres.shape[0]).astype('float32')

    pres_layer = utils.PressureLayer(1000.0, 500.0)
    hght_layer = utils.HeightLayer(0.0, 3000.0)
    strm_motnv = winds.WindComponents()
    strm_motnv.u = 10.0
    strm_motnv.v = 0.0

    comp1 = winds.mean_wind(pres_layer, pres, uwin, vwin)
    comp2 = winds.mean_wind_npw(pres_layer, pres, uwin, vwin)
    print(comp1.u, comp1.v)
    print(comp2.u, comp2.v)

    shear = winds.wind_shear(pres_layer, pres, uwin, vwin)
    print(shear.u, shear.v)

    helicity = winds.helicity(hght_layer, strm_motnv, hght, uwin, vwin)
    print(helicity)


def main():

    test_interp()
    test_utils()
    test_winds()

if __name__ == "__main__":
    main()
