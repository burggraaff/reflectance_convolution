"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from bandaveraging import plot_bands, load_data, calculate_differences, boxplot_absolute, boxplot_relative

wavelengths_data, Ed, Lw, R_rs = load_data()

for sat, filename in zip(["Aqua", "Terra"], ["spectral_response/HMODISA_RSRs.txt", "spectral_response/HMODIST_RSRs.txt"]):
    satellite_label = f"MODIS {sat}"
    band_labels = [f"{wvl} nm" for wvl in [412, 443, 469, 488, 531, 551, 555, 645, 667, 678, 748]]
    wavelengths_modis, *responses = np.loadtxt(filename, skiprows=8, unpack=True, usecols=np.arange(12))
    responses = np.array([responses])[0]
    ind = np.where(wavelengths_modis <= 900)[0]
    wavelengths_modis = wavelengths_modis[ind]
    responses = responses[:, ind]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:blue", "xkcd:cyan", "xkcd:bright green", "xkcd:green", "xkcd:forest green", "xkcd:red", "xkcd:dark red", "xkcd:dark brown", "k"]

    plot_bands(wavelengths_modis, responses, band_labels=band_labels, colours=colours, sensor_label=satellite_label)

    difference_absolute, difference_relative = calculate_differences(wavelengths_modis, responses, wavelengths_data, Ed, Lw, R_rs)

    boxplot_relative(difference_relative, colours=colours, band_labels=band_labels, sensor_label=satellite_label)
    boxplot_absolute(difference_absolute, colours=colours, band_labels=band_labels, sensor_label=satellite_label)