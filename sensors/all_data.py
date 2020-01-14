"""
Generate boxcar and gaussian spectral response functions
"""

from sba.bandaveraging import calculate_differences
from sba.io import load_data_file
from sba.response_curves import load_selected_sensors
import sys
from pathlib import Path

sensors = load_selected_sensors(*sys.argv[1:])

data_files = Path("data").glob("*processed.tab")

for file in data_files:
    for sensor in sensors:
        label, wavelengths_data, Ed, Lw, R_rs = load_data_file(file)

        reflectance_space = sensor.band_average(wavelengths_data, R_rs)
        radiance_space = sensor.band_average(wavelengths_data, Lw) / sensor.band_average(wavelengths_data, Ed)

        difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)

        sensor.boxplot_absolute(difference_absolute, data_label=label)
        sensor.boxplot_relative(difference_relative, data_label=label)
