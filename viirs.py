"""
Generate boxcar and gaussian spectral response functions
"""

from bandaveraging import load_data, bandaverage_multi, boxplot_absolute, boxplot_relative
import response_curves as rc

wavelengths_data, Ed, Lw, R_rs = load_data()

VIIRS = rc.load_VIIRS()

VIIRS.plot()

reflectance_space = VIIRS.band_average(wavelengths_data, R_rs)
radiance_space = VIIRS.band_average(wavelengths_data, Lw) / VIIRS.band_average(wavelengths_data, Ed)
difference_absolute = reflectance_space - radiance_space
difference_relative = 100*difference_absolute / radiance_space

VIIRS.boxplot_relative(difference_relative)
VIIRS.boxplot_absolute(difference_absolute)
