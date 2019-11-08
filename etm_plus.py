"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from bandaveraging import load_data, bandaverage_multi, boxplot_absolute, boxplot_relative

wavelengths_data, Ed, Lw, R_rs = load_data()

bands = [1,2,3,4,8]
band_labels = [f"ETM+\nband {band}" for band in bands]
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

reflectance_space = np.array([bandaverage_multi(response[0], response[1], wavelengths_data, R_rs) for response in responses])
radiance_space = np.array([bandaverage_multi(response[0], response[1], wavelengths_data, Lw) for response in responses]) / np.array([bandaverage_multi(response[0], response[1], wavelengths_data, Ed) for response in responses])
difference_absolute = reflectance_space - radiance_space
difference_relative = 100*difference_absolute / radiance_space

boxplot_relative(difference_relative, colours=colours, band_labels=band_labels, sensor_label="ETM+")
boxplot_absolute(difference_absolute, colours=colours, band_labels=band_labels, sensor_label="ETM+")
