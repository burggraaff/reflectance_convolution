"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io.ascii import read
from astropy import table
from mpl_toolkits.axes_grid1 import make_axes_locatable
from bandaveraging import split_spectrum, bandaverage_multi

bands = [1,2,3,4,8]
band_labels = [f"ETM+\nband {band}" for band in bands]
responses = [np.loadtxt(f"spectral_response/ETM+/spectral_b{x}.dat", skiprows=3, unpack=True) for x in bands]
colours = ["b", "g", "r", "xkcd:dark red", "k"]

for response, band, colour in zip(responses, bands, colours):
    plt.plot(response[0], response[1], label=f"ETM+ band {band}", c=colour)
plt.xlim(400, 920)
plt.ylim(0, 1.01)
plt.xlabel("Wavelength [nm]")
plt.ylabel("Relative response")
plt.legend(loc="best")
plt.show()

data_norcohab = read("data/norcohab_processed.tab")
data_archemhab = read("data/archemhab_processed.tab")

data_all = table.vstack([data_norcohab, data_archemhab])

wavelengths, Ed = split_spectrum(data_all, "Ed")
wavelengths, Lw = split_spectrum(data_all, "Lw")
wavelengths, R_rs = split_spectrum(data_all, "R_rs")

reflectance_space = np.array([bandaverage_multi(response[0], response[1], wavelengths, R_rs) for response in responses])
radiance_space = np.array([bandaverage_multi(response[0], response[1], wavelengths, Lw) for response in responses]) / np.array([bandaverage_multi(response[0], response[1], wavelengths, Ed) for response in responses])
difference_absolute = reflectance_space - radiance_space
difference_relative = 100*difference_absolute / radiance_space

for difference_set, unit in zip([difference_absolute, difference_relative], ["sr$^{-1}$", "%"]):
    bplot = plt.boxplot(difference_set.T, vert=False, showfliers=False, whis=[5,95], patch_artist=True, labels=band_labels)
    for patch, colour in zip(bplot["boxes"], colours):
        patch.set_facecolor(colour)
    plt.xlabel(f"Difference [{unit}]")
    plt.title("Landsat 7 ETM+")
    plt.grid(ls="--", color="0.5")
    plt.show()

#    for diff, band, colour in zip(difference_set, bands, colours):
#        vmin, vmax = np.nanpercentile(diff, 0.5), np.nanpercentile(diff, 99.5)
#        ME = np.median(diff)
#        MAE = np.median(np.abs(diff))
#        plt.hist(diff, bins=np.linspace(vmin, vmax, 100), color=colour)
#        plt.xlim(vmin, vmax)
#        plt.xlabel(f"Difference [{unit}]")
#        plt.ylabel("Frequency")
#        plt.title(f"ETM+ band {band} \n Med. Err. {ME:+.6f} {unit}; Med. Abs. Err. {MAE:.6f} {unit}")
#        plt.show()