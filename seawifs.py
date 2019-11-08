"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io.ascii import read
from astropy import table
from bandaveraging import split_spectrum, bandaverage_multi, plot_bands

band_labels = [f"{wvl} nm" for wvl in [412, 443, 490, 510, 555, 670, 765]]
wavelengths_seawifs, *responses_raw = np.loadtxt("spectral_response/SeaWiFS_RSRs.txt", skiprows=9, unpack=True, usecols=np.arange(8))
responses_raw = np.array(responses_raw)
responses = responses_raw / responses_raw.max(axis=1)[:, np.newaxis]
colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:cyan", "xkcd:lime green", "xkcd:forest green", "xkcd:dark red", "k"]

plot_bands(wavelengths_seawifs, responses, band_labels=band_labels, colours=colours, sensor_label="SeaWiFS")

data_norcohab = read("data/norcohab_processed.tab")
data_archemhab = read("data/archemhab_processed.tab")

data_all = table.vstack([data_norcohab, data_archemhab])
data_all = data_norcohab

wavelengths, Ed = split_spectrum(data_all, "Ed")
wavelengths, Lw = split_spectrum(data_all, "Lw")
wavelengths, R_rs = split_spectrum(data_all, "R_rs")

reflectance_space = np.array([bandaverage_multi(wavelengths_seawifs, response, wavelengths, R_rs) for response in responses])
radiance_space = np.array([bandaverage_multi(wavelengths_seawifs, response, wavelengths, Lw) for response in responses]) / np.array([bandaverage_multi(wavelengths_seawifs, response, wavelengths, Ed) for response in responses])
difference_absolute = reflectance_space - radiance_space
difference_relative = 100*difference_absolute / radiance_space

for difference_set, unit, label in zip([difference_absolute, difference_relative], ["sr$^{-1}$", "%"], ["abs", "rel"]):
    bplot = plt.boxplot(difference_set.T, vert=False, showfliers=False, whis=[5,95], patch_artist=True, labels=band_labels)
    for patch, colour in zip(bplot["boxes"], colours):
        patch.set_facecolor(colour)
    plt.xlabel(f"Difference [{unit}]")
    plt.title("SeaWiFS")
    plt.grid(ls="--", color="0.5")
    plt.savefig(f"results/SeaWiFS_{label}.pdf")
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