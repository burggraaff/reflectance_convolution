"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io.ascii import read
from astropy import table
from bandaveraging import split_spectrum, bandaverage_multi

band_labels = [f"{wvl} nm" for wvl in [443, 520, 550, 670]]
wavelengths_czcs, *responses = np.loadtxt("spectral_response/CZCS_RSRs.txt", skiprows=56, unpack=True)
responses = np.array(responses)
responses[responses < 0] = 0
colours = ["xkcd:dark blue", "xkcd:lime green", "xkcd:forest green", "xkcd:dark red"]

for response, band_label, colour in zip(responses, band_labels, colours):
    plt.plot(wavelengths_czcs, response, label=band_label, c=colour)
plt.xlim(380, 800)
plt.ylim(0, 1.01)
plt.xlabel("Wavelength [nm]")
plt.ylabel("Relative response")
plt.title("CZCS")
plt.legend(loc="best")
plt.show()

data_norcohab = read("data/norcohab_processed.tab")
data_archemhab = read("data/archemhab_processed.tab")

data_all = table.vstack([data_norcohab, data_archemhab])
data_all = data_norcohab

wavelengths, Ed = split_spectrum(data_all, "Ed")
wavelengths, Lw = split_spectrum(data_all, "Lw")
wavelengths, R_rs = split_spectrum(data_all, "R_rs")

reflectance_space = np.array([bandaverage_multi(wavelengths_czcs, response, wavelengths, R_rs) for response in responses])
radiance_space = np.array([bandaverage_multi(wavelengths_czcs, response, wavelengths, Lw) for response in responses]) / np.array([bandaverage_multi(wavelengths_czcs, response, wavelengths, Ed) for response in responses])
difference_absolute = reflectance_space - radiance_space
difference_relative = 100*difference_absolute / radiance_space

for difference_set, unit, label in zip([difference_absolute, difference_relative], ["sr$^{-1}$", "%"], ["abs", "rel"]):
    lower_percentile = np.nanpercentile(difference_set, 15.9, axis=1)
    medians = np.nanmedian(difference_set, axis=1)
    upper_percentile = np.nanpercentile(difference_set, 84.1, axis=1)
    lower_error = medians - lower_percentile
    upper_error = upper_percentile - medians
    for band, med, low, up in zip(band_labels, medians, lower_error, upper_error):
        print(f"CZCS {band} band: {med} (+{up}, -{low}) {unit}")

    bplot = plt.boxplot(difference_set.T, vert=False, showfliers=False, whis=[5,95], patch_artist=True, labels=band_labels)
    for patch, colour in zip(bplot["boxes"], colours):
        patch.set_facecolor(colour)
    plt.xlabel(f"Difference [{unit}]")
    plt.title("CZCS")
    plt.grid(ls="--", color="0.5")
    plt.savefig(f"results/CZCS_{label}.pdf")
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