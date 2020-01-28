from sba.bandaveraging import calculate_differences
from sba.io import load_data
from sba.chla import satellite_algorithms, satellite_algorithm_labels
from matplotlib import pyplot as plt

label, wavelengths_data, Ed, Lw, R_rs = load_data()

difference_absolute, difference_relative = zip(*[calculate_differences(*func(wavelengths_data, Ed, Lw, R_rs)) for func in satellite_algorithms])

for difference, absrel in zip([difference_absolute, difference_relative], ["abs", "rel"]):
    unit = "mg m$^{-3}$" if absrel == "abs" else "%"

    plt.figure(figsize=(5,3), tight_layout=True)
    bp = plt.boxplot(difference, whis=[5,95], showfliers=False, labels=satellite_algorithm_labels, patch_artist=True)
    for patch in bp["boxes"]:
        patch.set_facecolor("xkcd:dark green")
    plt.ylabel(f"$\Delta$ Chl-a [{unit}]")
    plt.grid(ls="--")
    plt.axhline(0, ls="--", c="k")
    plt.tick_params(axis="x", rotation=20)
    plt.title(f"Convolution error in Chl-a estimates: {label}")
    plt.savefig(f"results/{label}/{label}_chla_{absrel}.pdf")
    plt.show()
