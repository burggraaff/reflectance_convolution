"""
Generate boxcar and gaussian spectral response functions
"""

from bandaveraging import load_data
from matplotlib import pyplot as plt
import response_curves as rc

wavelengths_data, Ed, Lw, R_rs = load_data()

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

for sensor in sensors:
    sensor.plot()

    reflectance_space = sensor.band_average(wavelengths_data, R_rs)
    radiance_space = sensor.band_average(wavelengths_data, Lw) / sensor.band_average(wavelengths_data, Ed)
    difference_absolute = reflectance_space - radiance_space
    difference_relative = 100*difference_absolute / radiance_space

    sensor.boxplot_relative(difference_relative)
    sensor.boxplot_absolute(difference_absolute)
