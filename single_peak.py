"""
Bias as function of distance between peak center and band center, at several FWHMs
"""

import numpy as np
from matplotlib import pyplot as plt
from sba.response_curves import generate_boxcar
from sba.bandaveraging import calculate_differences

wavelengths = np.arange(0, 100, 0.1)

fwhm = 5
sigma = fwhm/2.355

R_rs = 1.2-wavelengths/100 + np.exp(-(wavelengths-50)**2 / (2 * sigma**2))

Ed = 1-0.5*(np.exp(-(wavelengths-25)**2 / (2 * sigma**2)) + np.exp(-(wavelengths-75)**2 / (2 * sigma**2)))

Lw = Ed * R_rs

central_wavelengths = np.arange(0, 101)
FWHMs = np.arange(1, 15)

results = np.tile(np.nan, (len(central_wavelengths), len(FWHMs)))
for i,center_sensor in enumerate(central_wavelengths):
    for j,fwhm_sensor in enumerate(FWHMs):
        sensor = generate_boxcar(center_sensor, fwhm_sensor)

        reflectance_space = sensor.band_average(wavelengths, [R_rs])
        radiance_space = sensor.band_average(wavelengths, [Lw]) / sensor.band_average(wavelengths, [Ed])

        difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)

        results[i,j] = difference_relative[0,0]

plt.plot(central_wavelengths, results[:,0])
plt.plot(central_wavelengths, results[:,-1])
