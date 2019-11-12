"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from bandaveraging import plot_bands, load_data, calculate_differences, boxplot_absolute, boxplot_relative

wavelengths_data, Ed, Lw, R_rs = load_data()

for sat, filename in zip(["A", "B"], ["spectral_response/MSI/MSI_S2A.txt", "spectral_response/MSI/MSI_S2B.txt"]):
    satellite_label = f"Sentinel-2{sat}"
    band_labels = [f"B{nr}" for nr in range(1,8)]
    wavelengths_msi, *responses = np.loadtxt(filename, skiprows=1, unpack=True)
    responses = np.array([responses])[0]
    colours = ["xkcd:blue", "xkcd:cyan", "xkcd:green", "xkcd:red", "xkcd:dark red", "xkcd:dark brown", "k"]

    plot_bands(wavelengths_msi, responses, band_labels=band_labels, colours=colours, sensor_label=satellite_label)

    difference_absolute, difference_relative = calculate_differences(wavelengths_msi, responses, wavelengths_data, Ed, Lw, R_rs)

    boxplot_relative(difference_relative, colours=colours, band_labels=band_labels, sensor_label=satellite_label)
    boxplot_absolute(difference_absolute, colours=colours, band_labels=band_labels, sensor_label=satellite_label)
