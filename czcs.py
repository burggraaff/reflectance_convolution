"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from bandaveraging import plot_bands, load_data, calculate_differences, boxplot_absolute, boxplot_relative

wavelengths_data, Ed, Lw, R_rs = load_data()

band_labels = [f"{wvl} nm" for wvl in [443, 520, 550, 670]]
wavelengths_czcs, *responses = np.loadtxt("spectral_response/CZCS_RSRs.txt", skiprows=56, unpack=True)
responses = np.array(responses)
responses[responses < 0] = 0
colours = ["xkcd:dark blue", "xkcd:lime green", "xkcd:forest green", "xkcd:dark red"]

plot_bands(wavelengths_czcs, responses, band_labels=band_labels, colours=colours, sensor_label="CZCS")

difference_absolute, difference_relative = calculate_differences(wavelengths_czcs, responses, wavelengths_data, Ed, Lw, R_rs)

boxplot_relative(difference_relative, colours=colours, band_labels=band_labels, sensor_label="CZCS")
boxplot_absolute(difference_absolute, colours=colours, band_labels=band_labels, sensor_label="CZCS")
