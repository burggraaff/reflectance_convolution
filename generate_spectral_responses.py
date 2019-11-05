"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io.ascii import read
from astropy import table

data_wavelengths = np.arange(300, 800, 0.1)

central_wavelengths = np.arange(350, 780, 3)
FWHMs = np.arange(1, 100, 1)

boxcar_result = np.tile(np.nan, [len(FWHMs), len(central_wavelengths)])
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

data_all = table.vstack([data_norcohab, data_archemhab])

for i,center in enumerate(central_wavelengths):
    for j,fwhm in enumerate(FWHMs):

        # Boxcar
        if center-fwhm/2 < data_wavelengths[0] or center+fwhm/2 > data_wavelengths[-1]:
            # Skip combination if the boxcar response does not fall entirely
            # within the data wavelength range
            continue
        boxcar_response = generate_boxcar(data_wavelengths, center, fwhm)
        boxcar_result[j,i] = 1

        if center-1.5*fwhm < data_wavelengths[0] or center+1.5*fwhm > data_wavelengths[-1]:
            # Skip combination if the gaussian response does not fall
            # within the data wavelength range up to 3 stds out
            continue
        gaussian_response = generate_gaussian(data_wavelengths, center, fwhm)
        gauss_result[j,i] = 1


plt.imshow(boxcar_result, origin="lower", extent=[central_wavelengths[0], central_wavelengths[-1], FWHMs[0], FWHMs[-1]])
plt.show()

plt.imshow(gauss_result, origin="lower", extent=[central_wavelengths[0], central_wavelengths[-1], FWHMs[0], FWHMs[-1]])
plt.show()

# Loop over central wavelengths
    # Loop over bandwidths
    # Skip if (significantly) out-of-range compared to data
    # Generate boxcar, gaussian
    # Store in big array

# Separate script: import responses
# Load spectra
# Convolve, calculate difference