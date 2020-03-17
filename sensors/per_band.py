from sba.bandaveraging import calculate_differences
from sba.io import load_data_file
from sba.response_curves import load_all_sensors
from pathlib import Path
from matplotlib import pyplot as plt
import warnings
import numpy as np

sensors = load_all_sensors()

data_files = Path("data").glob("*processed.tab")

labels, wavelengths, Eds, Lws, R_rss = zip(*[load_data_file(file) for file in data_files])

def boxplot(band, diff_abs, diff_rel, labels, saveto="boxplot.pdf", sensor_name=""):
    fig, axs = plt.subplots(nrows=2, figsize=(7,2), sharex=True, gridspec_kw={"hspace": 0.05, "wspace": 0})
    for diff, ax in zip([diff_abs, diff_rel], axs):
        diffs_including_all = [[x for y in diff for x in y], *diff]
        as_arrays = [np.array(d) for d in diffs_including_all]
        without_nan = [d[~np.isnan(d)] for d in as_arrays]
        labels_all = ["all", *labels]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            bplot = ax.boxplot(without_nan, showfliers=False, whis=[5,95], patch_artist=True, labels=labels_all)
        for patch in bplot["boxes"]:
            patch.set_facecolor(band.colour)

        ax.grid(ls="--")
        ax.axhline(0, c="k", ls="--")
        ax.tick_params(axis="x", rotation=90)

    axs[0].set_ylabel("[$10^{-6}$ sr$^{-1}$]")
    axs[0].tick_params(bottom=False, labelbottom=False)
    axs[1].set_ylabel(r"$\Delta \bar R_{rs}$ [%]")

    axs[0].set_title(f"Convolution error in the {sensor_name} {band} band")

    fig.align_labels()

    plt.savefig(saveto, bbox_inches="tight")
    plt.show()
    plt.close()


for sensor in sensors:
    print(sensor)
    for band in sensor.bands:
        print(f"     {band}")

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            RrsR = [band.convolve(wavelengths_data, R_rs) for wavelengths_data, R_rs in zip(wavelengths, R_rss)]
            RrsL = [band.convolve(wavelengths_data, Lw)/band.convolve(wavelengths_data, Ed) for wavelengths_data, Ed, Lw in zip(wavelengths, Eds, Lws)]

            difference_absolute, difference_relative = zip(*[calculate_differences(reflectance_space, radiance_space) for reflectance_space, radiance_space in zip(RrsR, RrsL)])

        # Convert to 10^-6 sr
        difference_absolute = [1e6 * diff for diff in difference_absolute]

        short_name = band.label.replace('\n', '_').replace(" ", "_")

        boxplot(band, difference_absolute, difference_relative, labels, saveto=f"results/per_band/{sensor.name}_{short_name}.pdf", sensor_name=sensor.name)
