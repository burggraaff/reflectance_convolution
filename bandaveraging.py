"""
Module with functions for band-averaging
"""

import numpy as np

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
    wavelength_step = band_wavelengths[1] - band_wavelengths[0]
    response_sum = wavelength_step * response_multiplied.sum()
    weight_sum = wavelength_step * band_response.sum()
    response_average = response_sum / weight_sum
    return response_average

def bandaverage_multi(band_wavelengths, band_response, data_wavelengths, data_response_multi):
    response_interpolated = np.array([np.interp(band_wavelengths, data_wavelengths, data_response, left=0, right=0) for data_response in data_response_multi])
    response_multiplied = response_interpolated * band_response
    wavelength_step = band_wavelengths[1] - band_wavelengths[0]
    response_sum = wavelength_step * response_multiplied.sum(axis=1)
    weight_sum = wavelength_step * band_response.sum()
    response_average = response_sum / weight_sum
    return response_average
