"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from bandaveraging import plot_bands, load_data, calculate_differences, boxplot_absolute, boxplot_relative

wavelengths_data, Ed, Lw, R_rs = load_data()

band_labels = [f"M{j}" for j in np.arange(1,7)]
wavelengths_viirs, *responses_raw = np.loadtxt("spectral_response/VIIRSN_IDPSv3_RSRs.txt", skiprows=5, unpack=True, usecols=np.arange(7))
responses_raw = np.array(responses_raw)
responses = responses_raw / responses_raw.max(axis=1)[:, np.newaxis]
colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:cyan", "xkcd:forest green", "xkcd:dark red", "k"]

plot_bands(wavelengths_viirs, responses, band_labels=band_labels, colours=colours, sensor_label="VIIRS")

difference_absolute, difference_relative = calculate_differences(wavelengths_viirs, responses, wavelengths_data, Ed, Lw, R_rs)

boxplot_relative(difference_relative, colours=colours, band_labels=band_labels, sensor_label="VIIRS")
boxplot_absolute(difference_absolute, colours=colours, band_labels=band_labels, sensor_label="VIIRS")
