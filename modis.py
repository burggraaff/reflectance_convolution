"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io.ascii import read
from astropy import table
from bandaveraging import split_spectrum, bandaverage_multi

data_norcohab = read("data/norcohab_processed.tab")
data_archemhab = read("data/archemhab_processed.tab")

data_all = table.vstack([data_norcohab, data_archemhab])

wavelengths, Ed = split_spectrum(data_all, "Ed")
wavelengths, Lw = split_spectrum(data_all, "Lw")
wavelengths, R_rs = split_spectrum(data_all, "R_rs")

for sat, filename in zip(["Aqua", "Terra"], ["spectral_response/HMODISA_RSRs.txt", "spectral_response/HMODIST_RSRs.txt"]):
    satellite_label = f"MODIS {sat}"
    band_labels = [f"{wvl} nm" for wvl in [412, 443, 469, 488, 531, 551, 555, 645, 667, 678, 748]]
    wavelengths_modis, *responses = np.loadtxt(filename, skiprows=8, unpack=True, usecols=np.arange(11))
    responses = np.array([responses])[0]
    ind = np.where(wavelengths_modis <= 900)[0]
    wavelengths_modis = wavelengths_modis[ind]
    responses = responses[:, ind]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:blue", "xkcd:cyan", "xkcd:bright green", "xkcd:green", "xkcd:forest green", "xkcd:red", "xkcd:dark red", "xkcd:dark brown", "k"]

    for response, band_label, colour in zip(responses, band_labels, colours):
        plt.plot(wavelengths_modis, response, label=band_label, c=colour)
    plt.xlim(380, 720)
    plt.ylim(0, 1.01)
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Relative response")
    plt.title(satellite_label)
    plt.legend(loc="best")
    plt.show()

    reflectance_space = np.array([bandaverage_multi(wavelengths_modis, response, wavelengths, R_rs) for response in responses])
    radiance_space = np.array([bandaverage_multi(wavelengths_modis, response, wavelengths, Lw) for response in responses]) / np.array([bandaverage_multi(wavelengths_modis, response, wavelengths, Ed) for response in responses])
    difference_absolute = reflectance_space - radiance_space
    difference_relative = 100*difference_absolute / radiance_space

    for difference_set, unit in zip([difference_absolute, difference_relative], ["sr$^{-1}$", "%"]):
        bplot = plt.boxplot(difference_set.T, vert=False, showfliers=False, whis=[5,95], patch_artist=True)
        for patch, colour in zip(bplot["boxes"], colours):
            patch.set_facecolor(colour)
        plt.xlabel(f"Difference [{unit}]")
        plt.yticks(np.arange(1,11), band_labels)
        plt.title(satellite_label)
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
