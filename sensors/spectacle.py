from sba.bandaveraging import calculate_differences
from sba.io import load_data_file
from sba.response_curves import load_SPECTACLE
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import warnings

sensor = load_SPECTACLE()

data_files = Path("data").glob("*processed.tab")

labels, wavelengths, Eds, Lws, R_rss = zip(*[load_data_file(file) for file in data_files])

def get_differences(band):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        RrsR = [band.convolve(wavelengths_data, R_rs) for wavelengths_data, R_rs in zip(wavelengths, R_rss)]
        RrsL = [band.convolve(wavelengths_data, Lw)/band.convolve(wavelengths_data, Ed) for wavelengths_data, Ed, Lw in zip(wavelengths, Eds, Lws)]

        difference_absolute, difference_relative = zip(*[calculate_differences(reflectance_space, radiance_space) for reflectance_space, radiance_space in zip(RrsR, RrsL)])

    # Get all differences into one array
    difference_absolute = np.array([x for y in difference_absolute for x in y]) * 1e6
    difference_relative = np.array([x for y in difference_relative for x in y])
    difference_combined = np.stack([difference_absolute, difference_relative])

    print(band)
    return difference_combined

differences = np.array([get_differences(band) for band in sensor.bands])
differences = np.moveaxis(differences, 1, 0)

labels =["R", "G\nDJI Phantom Pro 4", "B", "R", "G\nApple iPhone SE", "B", "R", "G\nSamsung Galaxy S8", "B"]

fig, axs = plt.subplots(nrows=2, figsize=(7,3), tight_layout=True, sharex=True, gridspec_kw={"hspace": 0, "wspace": 0})
for ax, diff in zip(axs, differences):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bplot = ax.boxplot(diff.T, showfliers=False, whis=[5,95], patch_artist=True, labels=labels)
    for patch, colour in zip(bplot["boxes"], sensor.get_band_colours()):
        patch.set_facecolor(colour)

    ax.grid(ls="--")
    ax.axhline(0, c="k", ls="--")

    axs[0].set_ylabel("[$10^{-6}$ sr$^{-1}$]")
    axs[1].set_ylabel(r"$\Delta \bar R_{rs}$ [%]")

    axs[0].set_title("SPECTACLE low-cost sensors, all data")

fig.align_labels()

plt.savefig("results/all/SPECTACLE_nice.pdf")
plt.show()
plt.close()
