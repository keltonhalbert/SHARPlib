import numpy as np
import nwsspc.sharp.calc.interp as interp


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

def main():

    test_interp()

if __name__ == "__main__":
    main()
