"""
Generate boxcar and gaussian spectral response functions
"""

from sba.bandaveraging import load_data
import sba.response_curves as rc

label, wavelengths_data, Ed, Lw, R_rs = load_data()

CZCS = rc.load_CZCS()

CZCS.plot()

reflectance_space = CZCS.band_average(wavelengths_data, R_rs)
radiance_space = CZCS.band_average(wavelengths_data, Lw) / CZCS.band_average(wavelengths_data, Ed)
difference_absolute = reflectance_space - radiance_space
difference_relative = 100*difference_absolute / radiance_space

CZCS.boxplot_relative(difference_relative, data_label=label)
CZCS.boxplot_absolute(difference_absolute, data_label=label)
