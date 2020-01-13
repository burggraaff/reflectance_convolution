"""
Apply synthetic (boxcar and gaussian) spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from sba.bandaveraging import calculate_differences
from sba.io import load_data
from sba.response_curves import read_synthetic_sensor_type, load_synthetic_sensor
from sba.plotting import synthetic_sensor_contourf

label, wavelengths_data, Ed, Lw, R_rs = load_data()
sensor_type = read_synthetic_sensor_type()

wavelengths_central = np.arange(330, 800, 5)
FWHMs = np.concatenate([np.arange(1, 30, 1), np.arange(30, 82, 2)])

result_absolute = np.tile(np.nan, [len(FWHMs), len(wavelengths_central)])
result_relative = result_absolute.copy()

for i,center in enumerate(wavelengths_central):
    print(f"Central wavelength: {center} nm")
    for j,fwhm in enumerate(FWHMs):
        boxcar = load_synthetic_sensor(sensor_type, center, fwhm)
        reflectance_space = boxcar.band_average(wavelengths_data, R_rs)
        radiance_space = boxcar.band_average(wavelengths_data, Lw) / boxcar.band_average(wavelengths_data, Ed)

        difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)

        result_absolute[j,i] = np.median(difference_absolute)
        result_relative[j,i] = np.median(difference_relative)

for result, absrel in zip([result_absolute, result_relative], ["abs", "rel"]):
    synthetic_sensor_contourf(wavelengths_central, FWHMs, result, sensor_type=sensor_type, absrel=absrel, label=label)