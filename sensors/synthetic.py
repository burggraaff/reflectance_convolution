"""
Apply synthetic (boxcar and gaussian) spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from sba.bandaveraging import calculate_differences
from sba.io import load_data
from sba.response_curves import read_synthetic_sensor_type, load_synthetic_sensor
from sba.plotting import synthetic_sensor_contourf, synthetic_sensor_contourf_combined

label, wavelengths_data, Ed, Lw, R_rs = load_data()
sensor_type = read_synthetic_sensor_type()

wavelengths_central = np.arange(330, 800, 5)
FWHMs = np.concatenate([np.arange(1, 30, 1), np.arange(30, 82, 2)])

# Arrays for storing results - per FWHM, per central wavelength, absolute/relative
medians = np.tile(np.nan, [2, len(FWHMs), len(wavelengths_central)])  # Median
perc5s = np.copy(medians)  # 5th percentile
perc95s = np.copy(medians)  # 95th percentile

for i,center in enumerate(wavelengths_central):
    print(f"Central wavelength: {center} nm")
    for j,fwhm in enumerate(FWHMs):
        boxcar = load_synthetic_sensor(sensor_type, center, fwhm)
        reflectance_space = boxcar.band_average(wavelengths_data, R_rs)
        radiance_space = boxcar.band_average(wavelengths_data, Lw) / boxcar.band_average(wavelengths_data, Ed)

        difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)
        difference_combined = np.array([difference_absolute, difference_relative])[:,0]

        medians[:,j,i] = np.median(difference_combined, axis=1)
        perc5s[:,j,i] = np.percentile(difference_combined, 5, axis=1)
        perc95s[:,j,i] = np.percentile(difference_combined, 95, axis=1)

results_stacked = np.stack([perc5s, medians, perc95s])
results_absrel = np.moveaxis(results_stacked, 1, 0)
results_absrel[0] *= 1e6  # Convert to 10^-6 sr^-1
quantities = ["P5", "Median", "P95"]

for results_combined, absrel in zip(results_absrel, ["abs", "rel"]):
    for result, quantity in zip(results_combined, quantities):
        synthetic_sensor_contourf(wavelengths_central, FWHMs, result, sensor_type=sensor_type, absrel=absrel, label=label, quantity=quantity)

    synthetic_sensor_contourf_combined(wavelengths_central, FWHMs, results_combined, sensor_type=sensor_type, absrel=absrel, label=label, quantities=quantities)
