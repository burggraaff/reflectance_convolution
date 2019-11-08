"""
Module with functions for band-averaging
"""

import numpy as np
from matplotlib import pyplot as plt

def split_spectrum(data_table, label):
    keys_relevant = [key for key in data_table.keys() if label in key]
    wavelengths = np.array([float(key.split("_")[-1]) for key in keys_relevant])
    try:
        spectra = np.array([data_table[key]._data for key in keys_relevant]).T
    except AttributeError:
        spectra = np.array([data_table[key].data for key in keys_relevant]).T
    return wavelengths, spectra

def bandaverage(band_wavelengths, band_response, data_wavelengths, data_response):
    response_interpolated = np.interp(band_wavelengths, data_wavelengths, data_response, left=0, right=0)
    response_multiplied = response_interpolated * band_response
    response_sum = np.trapz(response_multiplied, x=band_wavelengths)
    weight_sum = np.trapz(band_response, x=band_wavelengths)
    response_average = response_sum / weight_sum
    return response_average

def bandaverage_multi(band_wavelengths, band_response, data_wavelengths, data_response_multi):
    response_interpolated = np.array([np.interp(band_wavelengths, data_wavelengths, data_response, left=0, right=0) for data_response in data_response_multi])
    response_multiplied = response_interpolated * band_response
    response_sum = np.trapz(response_multiplied, x=band_wavelengths, axis=1)
    weight_sum = np.trapz(band_response, x=band_wavelengths)
    response_average = response_sum / weight_sum
    return response_average

def plot_bands(wavelengths, responses, band_labels=None, colours=None, sensor_label=None):
    if colours is None:
        colours = ["k"] * len(responses)
    if band_labels is None:
        band_labels = [None] * len(responses)

    for response, band_label, colour in zip(responses, band_labels, colours):
        plt.plot(wavelengths, response, label=band_label, c=colour)
    plt.xlim(380, 800)
    plt.ylim(0, 1.01)
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Relative response")
    plt.title(sensor_label)
    plt.legend(loc="best")
    plt.grid(ls="--", color="0.5")
    plt.savefig(f"results/{sensor_label}_bands.pdf")
    plt.show()
    plt.close()
