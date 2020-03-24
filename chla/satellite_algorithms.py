from sba.bandaveraging import calculate_differences
from sba.io import load_data
from sba.chla import satellite_algorithms, satellite_algorithm_labels, satellite_algorithm_colours
from matplotlib import pyplot as plt

label, wavelengths_data, Ed, Lw, R_rs = load_data()

difference_absolute, difference_relative = zip(*[calculate_differences(*func(wavelengths_data, Ed, Lw, R_rs)) for func in satellite_algorithms])

fig, ax = plt.subplots(figsize=(7,1))
bp = ax.boxplot(difference_relative, whis=[5,95], showfliers=False, labels=satellite_algorithm_labels, patch_artist=True)
for patch, colour in zip(bp["boxes"], satellite_algorithm_colours):
    patch.set_facecolor(colour)
ax.set_ylabel("$\Delta$ Chl-a [%]")
ax.grid(ls="--")
ax.axhline(0, ls="--", c="k")
ax.locator_params(axis="y", nbins=5)
ax.set_title(f"Convolution error in retrieval algorithms, {label} data")


ax2 = ax.twinx()
ax2.set_ylabel("$\Delta$ TSM [%]")
ax2.boxplot(difference_relative, whis=[5,95], showfliers=False, labels=["" for i in satellite_algorithms], medianprops={"alpha": 1}, boxprops={"alpha": 1})
ax2.locator_params(axis="y", nbins=5)


plt.savefig(f"results/{label}/{label}_chla_rel.pdf", bbox_inches="tight")
plt.show()
plt.close()
