from sba.bandaveraging import calculate_differences
from sba.io import load_data_file
from sba.chla import HydroColor
from matplotlib import pyplot as plt
from pathlib import Path

data_files = Path("data").glob("*processed.tab")
labels, wavelengths, Eds, Lws, R_rss = zip(*[load_data_file(file) for file in data_files])

difference_absolute, difference_relative = zip(*[calculate_differences(*HydroColor(wavelengths_data, Ed, Lw, R_rs)) for wavelengths_data, Ed, Lw, R_rs in zip(wavelengths, Eds, Lws, R_rss)])

for difference, absrel in zip([difference_absolute, difference_relative], ["abs", "rel"]):
    unit = "NTU" if absrel == "abs" else "%"

    plt.figure(figsize=(7,2.5), tight_layout=True)
    bp = plt.boxplot(difference, whis=[5,95], showfliers=False, labels=labels, patch_artist=True, medianprops={"c": "k"})
    for patch in bp["boxes"]:
        patch.set_facecolor("xkcd:tan")
    plt.ylabel(f"$\Delta$ Turbidity [{unit}]")
    plt.grid(ls="--")
    plt.axhline(0, ls="--", c="k")
    plt.tick_params(axis="x", rotation=90)
    plt.title(f"Convolution error in HydroColor turbidity estimates")
    plt.savefig(f"results/chla/HydroColor_{absrel}.pdf")
    plt.show()

# Zoom
for difference, absrel in zip([difference_absolute, difference_relative], ["abs", "rel"]):
    unit = "NTU" if absrel == "abs" else "%"

    plt.figure(figsize=(7,2.5), tight_layout=True)
    bp = plt.boxplot(difference, whis=[5,95], showfliers=False, labels=labels, patch_artist=True, medianprops={"c": "k"})
    for patch in bp["boxes"]:
        patch.set_facecolor("xkcd:tan")
    plt.ylabel(f"$\Delta$ Turbidity [{unit}]")
    plt.grid(ls="--")
    plt.axhline(0, ls="--", c="k")
    plt.ylim(-10, 10)
    plt.tick_params(axis="x", rotation=90)
    plt.title(f"Convolution error in HydroColor turbidity estimates")
    plt.savefig(f"results/chla/HydroColor_{absrel}_zoom.pdf")
    plt.show()
