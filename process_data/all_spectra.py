import numpy as np
import matplotlib.pyplot as plt
from astropy.io.ascii import read
from pathlib import Path
from sba.bandaveraging import split_spectrum

folder = Path("data")

all_data_files = sorted(folder.glob("*processed.tab"))
colours = ["xkcd:royal blue", "xkcd:orange", "xkcd:lime green", "xkcd:magenta", "xkcd:olive green", "xkcd:navy", "xkcd:maroon", "xkcd:bright pink", "xkcd:hunter green", "xkcd:vomit", "xkcd:chocolate", "xkcd:rust brown", "xkcd:crimson", "xkcd:black"]

nr_spectra = 0

fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,7))

for file in all_data_files:
    data = read(file)
    label = file.stem[:-10]
    nr_spectra += len(data)

    wavelengths, Ed = split_spectrum(data, "Ed")
    wavelengths, Lw = split_spectrum(data, "Lw")
    wavelengths, R_rs = split_spectrum(data, "R_rs")

    for ax, spectrum in zip(axs.ravel(), [Ed, Lw, R_rs]):
        ax.plot(wavelengths, spectrum.T, c="k", alpha=0.01, zorder=1)

    print(label)

for ax, label in zip(axs.ravel(), ["$E_d$", "$L_w$", "$R_{rs}$"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[-1].set_xlabel("Wavelength [nm]")
axs[-1].set_xlim(300, 1000)

axs[0].set_ylim(ymin=0)
axs[1].set_ylim(ymin=0)
axs[2].set_ylim(ymin=0)

axs[0].set_title(f"All spectra ({nr_spectra})")
plt.savefig(f"data/plots/spectra_all_data.pdf")
plt.show()
plt.close()

# Zoom
fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,7))

for file in all_data_files:
    data = read(file)
    label = file.stem[:-10]

    wavelengths, Ed = split_spectrum(data, "Ed")
    wavelengths, Lw = split_spectrum(data, "Lw")
    wavelengths, R_rs = split_spectrum(data, "R_rs")

    for ax, spectrum in zip(axs.ravel(), [Ed, Lw, R_rs]):
        ax.plot(wavelengths, spectrum.T, c="k", alpha=0.01, zorder=1)

    print(label)

for ax, label in zip(axs.ravel(), ["$E_d$", "$L_w$", "$R_{rs}$"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[-1].set_xlabel("Wavelength [nm]")
axs[-1].set_xlim(300, 1000)

axs[0].set_ylim(0, 2)
axs[1].set_ylim(0, 0.05)
axs[2].set_ylim(0, 0.01)

axs[0].set_title(f"All spectra ({nr_spectra})")
plt.savefig(f"data/plots/spectra_all_data_zoom.pdf")
plt.show()
plt.close()