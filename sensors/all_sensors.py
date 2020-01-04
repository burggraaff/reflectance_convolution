"""
Generate boxcar and gaussian spectral response functions
"""

from sba.bandaveraging import load_data, calculate_differences
import sba.response_curves as rc
from sba import plotting as p
from matplotlib import pyplot as plt
import numpy as np

label, wavelengths_data, Ed, Lw, R_rs = load_data()

sensors = [func() for func in rc.functions]

# plot all response curves
fig, axs = plt.subplots(nrows=len(sensors), sharex=True, sharey=True, tight_layout=True, figsize=(3, 15), gridspec_kw={"hspace": 0, "wspace": 0})
for sensor, ax in zip(sensors, axs):
    sensor.plot(ax)
    ax.set_ylabel(sensor.name)
    ax.tick_params(axis="y", left=False, labelleft=False)

for ax in axs[:-1]:
    ax.tick_params(bottom=False, labelbottom=False)
axs[-1].set_xlabel("Wavelength [nm]")
axs[0].set_title("Relative spectral responses")
plt.savefig("results/all_bands.pdf")
plt.show()
plt.close()

reflectance_space = np.vstack([sensor.band_average(wavelengths_data, R_rs) for sensor in sensors])
radiance_space = np.vstack([sensor.band_average(wavelengths_data, Lw) / sensor.band_average(wavelengths_data, Ed)  for sensor in sensors])

difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)

labels = [sensor.sensor_band_labels for sensor in sensors]
labels = [label for sublist in labels for label in sublist]

colours = [sensor.colours for sensor in sensors]
colours = [colour for sublist in colours for colour in sublist]

p.boxplot_absolute(difference_absolute, band_labels=labels, colours=colours, sensor_label="", data_label=label)
p.boxplot_relative(difference_relative, band_labels=labels, colours=colours, sensor_label="", data_label=label)
