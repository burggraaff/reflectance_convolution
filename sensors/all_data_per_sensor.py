from sba.bandaveraging import calculate_differences
from sba.io import load_data_file
from sba.response_curves import load_all_sensors
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import warnings

sensors = load_all_sensors()

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

differences = [[get_differences(band) for band in sensor.bands] for sensor in sensors]

for sensor, diffs in zip(sensors, differences):
    fig, axs = plt.subplots(nrows=2, figsize=(7,3), tight_layout=True, sharex=True, gridspec_kw={"hspace": 0, "wspace": 0})
    diffs_array = np.array(diffs)
    diffs_array = np.moveaxis(diffs_array, 1, 0)
    for ax, diff in zip(axs, diffs_array):
        without_nan = [d[~np.isnan(d)] for d in diff]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bplot = ax.boxplot(without_nan, showfliers=False, whis=[5,95], patch_artist=True, labels=sensor.get_band_labels())
        for patch, colour in zip(bplot["boxes"], sensor.get_band_colours()):
            patch.set_facecolor(colour)

        ax.grid(ls="--")
        ax.axhline(0, c="k", ls="--")

        if len("".join(sensor.get_band_labels())) >= 40:
            ax.tick_params(axis="x", rotation=90)

    axs[0].set_ylabel("[$10^{-6}$ sr$^{-1}$]")
    axs[1].set_ylabel("$\Delta R_{rs}$ [%]")

    axs[0].set_title(f"All data ; {sensor.name}")

    fig.align_labels()

    plt.savefig(f"results/all/{sensor.name}.pdf")
    plt.show()
    plt.close()
