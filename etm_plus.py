"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io.ascii import read
from astropy import table
from mpl_toolkits.axes_grid1 import make_axes_locatable
from bandaveraging import split_spectrum, bandaverage_multi

bands = [1,2,3,4,8]
responses = [np.loadtxt(f"spectral_response/ETM+/spectral_b{x}.dat", skiprows=3, unpack=True) for x in bands]
colours = ["b", "g", "r", "xkcd:dark red", "k"]

for response, band, colour in zip(responses, bands, colours):
    plt.plot(response[0], response[1], label=f"ETM+ band {band}", c=colour)
plt.xlim(400, 920)
plt.ylim(0, 1.01)
plt.xlabel("Wavelength [nm]")
plt.ylabel("Relative response")
plt.legend(loc="best")
plt.show()

data_norcohab = read("data/norcohab_processed.tab")
data_archemhab = read("data/archemhab_processed.tab")

data_all = table.vstack([data_norcohab, data_archemhab])

wavelengths, Ed = split_spectrum(data_all, "Ed")
wavelengths, Lw = split_spectrum(data_all, "Lw")
wavelengths, R_rs = split_spectrum(data_all, "R_rs")

reflectance_space = np.array([bandaverage_multi(response[0], response[1], wavelengths, R_rs) for response in responses])
radiance_space = np.array([bandaverage_multi(response[0], response[1], wavelengths, Lw) for response in responses]) / np.array([bandaverage_multi(response[0], response[1], wavelengths, Ed) for response in responses])
difference = 100*(reflectance_space - radiance_space) / radiance_space

for diff, band, colour in zip(difference, bands, colours):
    plt.hist(diff, bins=np.linspace(-5, 5, 100), color=colour)
    plt.xlim(-5, 5)
    plt.xlabel("Difference (%)")
    plt.ylabel("Frequency")
    plt.title(f"ETM+ band {band}")
    plt.show()
