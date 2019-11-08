"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from bandaveraging import plot_bands, load_data, calculate_differences, boxplot_absolute, boxplot_relative
from pathlib import Path

wavelengths_data, Ed, Lw, R_rs = load_data()

files = list(Path("spectral_response/SPECTACLE/").glob("*.npy"))
devices = [file.stem for file in files]
bands = [[f"{device} {c}" for c in "RGB"] for device in devices]
band_labels = [band for sublist in bands for band in sublist]
colours = [*"rgb"] * len(devices)

wavelengths_spectacle = np.load(files[0])[0]
responses = np.vstack([np.load(file)[1:4] for file in files])

plot_bands(wavelengths_spectacle, responses, band_labels=band_labels, colours=colours, sensor_label="SPECTACLE")

difference_absolute, difference_relative = calculate_differences(wavelengths_spectacle, responses, wavelengths_data, Ed, Lw, R_rs)

boxplot_relative(difference_relative, colours=colours, band_labels=band_labels, sensor_label="SPECTACLE")
boxplot_absolute(difference_absolute, colours=colours, band_labels=band_labels, sensor_label="SPECTACLE")
