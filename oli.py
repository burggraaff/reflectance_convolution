"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from bandaveraging import load_data, bandaverage_multi, boxplot_absolute, boxplot_relative
from matplotlib import pyplot as plt

wavelengths_data, Ed, Lw, R_rs = load_data()

band_labels = ["CA", "Blue", "Green", "Red", "NIR", "Pan"]
responses = [np.loadtxt(f"spectral_response/OLI/{band}.txt", skiprows=1, unpack=True, usecols=[0,1]) for band in band_labels]
colours = ["xkcd:dark blue", "xkcd:blue", "xkcd:green", "xkcd:red", "xkcd:dark brown", "k"]

for response, band, colour in zip(responses, band_labels, colours):
    plt.plot(response[0], response[1], label=f"OLI band {band}", c=colour)
plt.xlim(400, 920)
plt.ylim(0, 1.01)
plt.xlabel("Wavelength [nm]")
plt.ylabel("Relative response")
plt.legend(loc="best")
plt.grid(ls="--", color="0.5")
plt.savefig(f"results/OLI_bands.pdf")
plt.show()

reflectance_space = np.array([bandaverage_multi(response[0], response[1], wavelengths_data, R_rs) for response in responses])
radiance_space = np.array([bandaverage_multi(response[0], response[1], wavelengths_data, Lw) for response in responses]) / np.array([bandaverage_multi(response[0], response[1], wavelengths_data, Ed) for response in responses])
difference_absolute = reflectance_space - radiance_space
difference_relative = 100*difference_absolute / radiance_space

boxplot_relative(difference_relative, colours=colours, band_labels=band_labels, sensor_label="OLI")
boxplot_absolute(difference_absolute, colours=colours, band_labels=band_labels, sensor_label="OLI")
