"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from bandaveraging import bandaverage_multi, plot_bands, load_data

wavelengths_data, Ed, Lw, R_rs = load_data()

for sat, filename in zip(["Aqua", "Terra"], ["spectral_response/HMODISA_RSRs.txt", "spectral_response/HMODIST_RSRs.txt"]):
    satellite_label = f"MODIS {sat}"
    band_labels = [f"{wvl} nm" for wvl in [412, 443, 469, 488, 531, 551, 555, 645, 667, 678, 748]]
    wavelengths_modis, *responses = np.loadtxt(filename, skiprows=8, unpack=True, usecols=np.arange(12))
    responses = np.array([responses])[0]
    ind = np.where(wavelengths_modis <= 900)[0]
    wavelengths_modis = wavelengths_modis[ind]
    responses = responses[:, ind]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:blue", "xkcd:cyan", "xkcd:bright green", "xkcd:green", "xkcd:forest green", "xkcd:red", "xkcd:dark red", "xkcd:dark brown", "k"]

    plot_bands(wavelengths_modis, responses, band_labels=band_labels, colours=colours, sensor_label=satellite_label)

    reflectance_space = np.array([bandaverage_multi(wavelengths_modis, response, wavelengths_data, R_rs) for response in responses])
    radiance_space = np.array([bandaverage_multi(wavelengths_modis, response, wavelengths_data, Lw) for response in responses]) / np.array([bandaverage_multi(wavelengths_modis, response, wavelengths_data, Ed) for response in responses])
    difference_absolute = reflectance_space - radiance_space
    difference_relative = 100*difference_absolute / radiance_space

    for difference_set, unit, label in zip([difference_absolute, difference_relative], ["sr$^{-1}$", "%"], ["abs", "rel"]):
        bplot = plt.boxplot(difference_set.T, vert=False, showfliers=False, whis=[5,95], patch_artist=True, labels=band_labels)
        for patch, colour in zip(bplot["boxes"], colours):
            patch.set_facecolor(colour)
        plt.xlabel(f"Difference [{unit}]")
        plt.title(satellite_label)
        plt.grid(ls="--", color="0.5")
        plt.savefig(f"results/{satellite_label}_{label}.pdf")
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
