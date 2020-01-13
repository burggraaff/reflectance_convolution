"""
Module with functions for band-averaging
"""

import numpy as np
from scipy.integrate import simps

integrate = np.trapz

def interpolate_spectral_data(band_wavelengths, data_wavelengths, data_response, extrapolation_value=np.nan):
    data_interpolated = np.interp(band_wavelengths, data_wavelengths, data_response, left=extrapolation_value, right=extrapolation_value)
    return data_interpolated


def bandaverage(band_wavelengths, band_response, data_wavelengths, data_response):
    response_interpolated = interpolate_spectral_data(band_wavelengths, data_wavelengths, data_response)
    response_multiplied = response_interpolated * band_response
    response_sum = integrate(response_multiplied, x=band_wavelengths)
    weight_sum = integrate(band_response, x=band_wavelengths)
    response_average = response_sum / weight_sum
    return response_average


def bandaverage_multi(band_wavelengths, band_response, data_wavelengths, data_response_multi):
    response_interpolated = np.array([interpolate_spectral_data(band_wavelengths, data_wavelengths, data_response) for data_response in data_response_multi])
    response_multiplied = response_interpolated * band_response
    response_sum = integrate(response_multiplied, x=band_wavelengths, axis=1)
    weight_sum = integrate(band_response, x=band_wavelengths)
    response_average = response_sum / weight_sum
    return response_average


def calculate_differences(reflectance_space, radiance_space):
    difference_absolute = reflectance_space - radiance_space
    difference_relative = 100 * difference_absolute / radiance_space
    difference_relative[difference_absolute == 0] = 0

    return difference_absolute, difference_relative


def calculate_median_and_errors(differences):
    lower_percentile = np.nanpercentile(differences, 15.9, axis=1)
    medians = np.nanmedian(differences, axis=1)
    upper_percentile = np.nanpercentile(differences, 84.1, axis=1)
    lower_error = medians - lower_percentile
    upper_error = upper_percentile - medians
    return medians, lower_error, upper_error
