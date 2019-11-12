"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
import xarray as xr
from bandaveraging import plot_bands, load_data, calculate_differences, boxplot_absolute, boxplot_relative, bandaverage_multi
from matplotlib import pyplot as plt

wavelengths_data, Ed, Lw, R_rs = load_data()

for sat, filename in zip(["A", "B"], ["spectral_response/OLCI/S3A_OL_SRF_20160713_mean_rsr.nc4", "spectral_response/OLCI/S3B_OL_SRF_0_20180109_mean_rsr.nc4"]):
    satellite_label = f"OLCI S3-{sat}"
    band_labels = [f"Oa{nr:02d}" for nr in range(1,17)]
    sr_olci = xr.open_dataset(filename)
    wavelengths_olci = sr_olci.mean_spectral_response_function_wavelength.data[:16]
    responses = sr_olci.mean_spectral_response_function.data[:16]
    colours = ["xkcd:dark violet", "xkcd:violet", "xkcd:purple", "xkcd:light blue", "xkcd:lime green", "xkcd:bright green", "xkcd:bright red", "xkcd:dark red", "xkcd:reddish brown", "xkcd:brown", "xkcd:dark brown", "k", "k", "k", "k", "k"]

    for wavelengths, response, band, colour in zip(wavelengths_olci, responses, band_labels, colours):
        plt.plot(wavelengths, response, label=f"{band}", c=colour)
    plt.xlim(380, 800)
    plt.ylim(0, 1.01)
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Relative response")
    plt.legend(loc="best")
    plt.grid(ls="--", color="0.5")
    plt.savefig(f"results/{satellite_label}_bands.pdf")
    plt.show()

    reflectance_space = np.array([bandaverage_multi(wavelengths, response, wavelengths_data, R_rs) for wavelengths, response in zip(wavelengths_olci, responses)])
    radiance_space = np.array([bandaverage_multi(wavelengths, response, wavelengths_data, Lw) for wavelengths, response in zip(wavelengths_olci, responses)]) / np.array([bandaverage_multi(wavelengths, response, wavelengths_data, Ed) for wavelengths, response in zip(wavelengths_olci, responses)])
    difference_absolute = reflectance_space - radiance_space
    difference_relative = 100*difference_absolute / radiance_space

    boxplot_relative(difference_relative, colours=colours, band_labels=band_labels, sensor_label=satellite_label)
    boxplot_absolute(difference_absolute, colours=colours, band_labels=band_labels, sensor_label=satellite_label)
