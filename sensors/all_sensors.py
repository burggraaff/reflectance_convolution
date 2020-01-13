"""
Generate boxcar and gaussian spectral response functions
"""

from sba.bandaveraging import calculate_differences
from sba.io import load_data
from sba.response_curves import load_all_sensors
from sba import plotting as p
import numpy as np

label, wavelengths_data, Ed, Lw, R_rs = load_data()

sensors = load_all_sensors()

reflectance_space = np.vstack([sensor.band_average(wavelengths_data, R_rs) for sensor in sensors])
radiance_space = np.vstack([sensor.band_average(wavelengths_data, Lw) / sensor.band_average(wavelengths_data, Ed)  for sensor in sensors])

difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)

labels = [sensor.sensor_band_labels for sensor in sensors]
labels = [label for sublist in labels for label in sublist]

colours = [sensor.colours for sensor in sensors]
colours = [colour for sublist in colours for colour in sublist]

p.boxplot_absolute(difference_absolute, band_labels=labels, colours=colours, sensor_label="All", data_label=label)
p.boxplot_relative(difference_relative, band_labels=labels, colours=colours, sensor_label="All", data_label=label)
