"""
Generate boxcar and gaussian spectral response functions
"""

from bandaveraging import load_data, bandaverage_multi, boxplot_absolute, boxplot_relative
import response_curves as rc

wavelengths_data, Ed, Lw, R_rs = load_data()

ETM_plus = rc.load_ETM_plus()

ETM_plus.plot()

reflectance_space = ETM_plus.band_average(wavelengths_data, R_rs)
radiance_space = ETM_plus.band_average(wavelengths_data, Lw) / ETM_plus.band_average(wavelengths_data, Ed)
difference_absolute = reflectance_space - radiance_space
difference_relative = 100*difference_absolute / radiance_space

ETM_plus.boxplot_relative(difference_relative)
ETM_plus.boxplot_absolute(difference_absolute)
