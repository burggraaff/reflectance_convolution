from sba.bandaveraging import calculate_differences
from sba.io import load_data_file
from sba.chla import GM09
from matplotlib import pyplot as plt
from pathlib import Path

data_files = Path("data").glob("*processed.tab")
labels, wavelengths, Eds, Lws, R_rss = zip(*[load_data_file(file) for file in data_files])

difference_absolute, difference_relative = zip(*[calculate_differences(*GM09(wavelengths_data, Ed, Lw, R_rs)) for wavelengths_data, Ed, Lw, R_rs in zip(wavelengths, Eds, Lws, R_rss)])

for difference, absrel in zip([difference_absolute, difference_relative], ["abs", "rel"]):
    unit = "mg m$^{-3}$" if absrel == "abs" else "%"

    plt.figure(figsize=(5,3), tight_layout=True)
    bp = plt.boxplot(difference, whis=[5,95], showfliers=False, labels=labels, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("xkcd:dark green")
    plt.ylabel(f"Difference in Chl-a [{unit}]")
    plt.grid(ls="--")
    plt.axhline(0, ls="--", c="k")
    plt.tick_params(axis="x", rotation=90)
    plt.title(f"Convolution error in GM09 Chl-a estimates")
    plt.savefig(f"results/chla/GM09_{absrel}.pdf")
    plt.show()
