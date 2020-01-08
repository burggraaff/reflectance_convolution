"""
Generate boxcar and gaussian spectral response functions
"""

from sba.io import load_data
from sba.response_curves import load_from_name

label, wavelengths_data, Ed, Lw, R_rs = load_data()
sensors = load_from_name()

for sensor in sensors:
    print(sensor)
    sensor.plot()

    reflectance_space = sensor.band_average(wavelengths_data, R_rs)
    radiance_space = sensor.band_average(wavelengths_data, Lw) / sensor.band_average(wavelengths_data, Ed)
    difference_absolute = reflectance_space - radiance_space
    difference_relative = 100*difference_absolute / radiance_space

    sensor.boxplot_relative(difference_relative, data_label=label)
    sensor.boxplot_absolute(difference_absolute, data_label=label)
