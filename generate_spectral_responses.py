"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io.ascii import read
from astropy import table
from mpl_toolkits.axes_grid1 import make_axes_locatable

wavelengths_band = np.arange(320, 800, 0.1)

wavelengths_central = np.arange(350, 780, 3)
FWHMs = np.arange(1, 75, 1)

boxcar_result = np.tile(np.nan, [len(FWHMs), len(wavelengths_central)])
gauss_result = boxcar_result.copy()

wavelengths_interp = np.arange(380, 800.5, 0.5)

def interpolate_row(row, spectrum):
    spectrum_data = np.array([row[f"{spectrum}_{wvl}"] for wvl in wavelengths])
    spectrum_interpolated = np.interp(wavelengths_interp, wavelengths, spectrum_data, left=0, right=0)
    return spectrum_interpolated

def interpolate_table(data_table, spectrum):
    interpolated_data = np.array([interpolate_row(row, spectrum) for row in data_table])
    table_names = [f"{spectrum}_{wvl}" for wvl in wavelengths_interp]
    interpolated_table = table.Table(data=interpolated_data, names=table_names)
    return interpolated_table

def generate_boxcar(x, center, fwhm):
    response = np.zeros_like(x)
    response[np.abs(x - center) <= fwhm/2] = 1.
    return response

def generate_gaussian(x, center, fwhm):
    response = np.exp(-(x-center)**2 / (2 * fwhm**2))
    return response

data_norcohab = read("data/norcohab_processed.tab")
data_archemhab = read("data/archemhab_processed.tab")

#data_all = table.vstack([data_norcohab, data_archemhab])
data_all = data_norcohab

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

wavelengths, Ed = split_spectrum(data_all, "Ed")
wavelengths, Lw = split_spectrum(data_all, "Lw")
wavelengths, R_rs = split_spectrum(data_all, "R_rs")

for i,center in enumerate(wavelengths_central):
    for j,fwhm in enumerate(FWHMs):

        # Boxcar
        if center-fwhm/2 < wavelengths_band[0] or center+fwhm/2 > wavelengths_band[-1]:
            # Skip combination if the boxcar response does not fall entirely
            # within the data wavelength range
            continue
        boxcar_response = generate_boxcar(wavelengths_band, center, fwhm)
        boxcar_reflectance_space = bandaverage_multi(wavelengths_band, boxcar_response, wavelengths, R_rs)
        boxcar_radiance_space = bandaverage_multi(wavelengths_band, boxcar_response, wavelengths, Lw) / bandaverage_multi(wavelengths_band, boxcar_response, wavelengths, Ed)
        boxcar_result[j,i] = np.median((boxcar_radiance_space - boxcar_reflectance_space) / boxcar_radiance_space)

        if center-1.5*fwhm < wavelengths_band[0] or center+1.5*fwhm > wavelengths_band[-1]:
            # Skip combination if the gaussian response does not fall
            # within the data wavelength range up to 3 stds out
            continue
        gaussian_response = generate_gaussian(wavelengths_band, center, fwhm)
        gauss_reflectance_space = bandaverage_multi(wavelengths_band, gaussian_response, wavelengths, R_rs)
        gauss_radiance_space = bandaverage_multi(wavelengths_band, gaussian_response, wavelengths, Lw) / bandaverage_multi(wavelengths_band, gaussian_response, wavelengths, Ed)
        gauss_result[j,i] = np.median((gauss_radiance_space - gauss_reflectance_space) / gauss_radiance_space)

for result in [boxcar_result, gauss_result]:
    im = plt.imshow(100*result, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], aspect="auto")
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    divider = make_axes_locatable(im.axes)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel("Difference (Rad. space - Refl. space, %)")
    plt.show()

# Loop over central wavelengths
    # Loop over bandwidths
    # Skip if (significantly) out-of-range compared to data
    # Generate boxcar, gaussian
    # Store in big array

# Separate script: import responses
# Load spectra
# Convolve, calculate difference