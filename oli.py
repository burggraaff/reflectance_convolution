"""
Generate boxcar and gaussian spectral response functions
"""

from sba.bandaveraging import load_data, bandaverage_multi, boxplot_absolute, boxplot_relative
import sba.response_curves as rc

wavelengths_data, Ed, Lw, R_rs = load_data()

OLI = rc.load_OLI()

OLI.plot()

reflectance_space = OLI.band_average(wavelengths_data, R_rs)
radiance_space = OLI.band_average(wavelengths_data, Lw) / OLI.band_average(wavelengths_data, Ed)
difference_absolute = reflectance_space - radiance_space
difference_relative = 100*difference_absolute / radiance_space

OLI.boxplot_relative(difference_relative)
OLI.boxplot_absolute(difference_absolute)
