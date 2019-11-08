"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from bandaveraging import plot_bands, load_data, calculate_differences, boxplot_absolute, boxplot_relative

wavelengths_data, Ed, Lw, R_rs = load_data()

band_labels = [f"{wvl} nm" for wvl in [412, 443, 490, 510, 555, 670, 765]]
wavelengths_seawifs, *responses_raw = np.loadtxt("spectral_response/SeaWiFS_RSRs.txt", skiprows=9, unpack=True, usecols=np.arange(8))
responses_raw = np.array(responses_raw)
responses = responses_raw / responses_raw.max(axis=1)[:, np.newaxis]
colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:cyan", "xkcd:lime green", "xkcd:forest green", "xkcd:dark red", "k"]

plot_bands(wavelengths_seawifs, responses, band_labels=band_labels, colours=colours, sensor_label="SeaWiFS")

difference_absolute, difference_relative = calculate_differences(wavelengths_seawifs, responses, wavelengths_data, Ed, Lw, R_rs)

boxplot_relative(difference_relative, colours=colours, band_labels=band_labels, sensor_label="SeaWiFS")
boxplot_absolute(difference_absolute, colours=colours, band_labels=band_labels, sensor_label="SeaWiFS")
