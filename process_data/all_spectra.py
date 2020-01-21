import numpy as np
import matplotlib.pyplot as plt
from astropy.io.ascii import read
from pathlib import Path
from sba.data_processing import split_spectrum

folder = Path("data")

all_data_files = sorted(folder.glob("*processed.tab"))

data_all = [read(file) for file in all_data_files]
N_all = [len(data) for data in data_all]
N_total = sum(N_all)

labels = [file.stem[:-10] for file in all_data_files]

wavelength_range = np.arange(300, 1350, 5)
wavelengths_all = [split_spectrum(data, "Ed")[0] for data in data_all]
wavelength_extremes = [[wavelength_range[0], wavelength_range[-1]] for wavelength_range in wavelengths_all]
wavelengths_bincount = np.zeros_like(wavelength_range)
for N, extremes in zip(N_all, wavelength_extremes):
    wavelengths_bincount[np.where((wavelength_range >= extremes[0]-2.5) & (wavelength_range <= extremes[-1]+2.5))] += N

plt.figure(figsize=(4,2.5), tight_layout=True)
plt.step(wavelength_range, wavelengths_bincount, where="mid", c="k")
plt.xlim(wavelength_range[0], wavelength_range[-1]+5)
plt.ylim(ymin=0)
plt.xlabel("Wavelength [nm]")
plt.ylabel("No. spectra $N$")
plt.grid(ls="--")
plt.title("Spectral coverage of data set")
plt.savefig("data/plots/coverage.pdf")
plt.show()
plt.close()

fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,7))

for data in data_all:
    wavelengths, Ed = split_spectrum(data, "Ed")
    wavelengths, Lw = split_spectrum(data, "Lw")
    wavelengths, R_rs = split_spectrum(data, "R_rs")

    for ax, spectrum in zip(axs.ravel(), [Ed, Lw, R_rs]):
        ax.plot(wavelengths, spectrum.T, c="k", alpha=0.01, zorder=1)

for ax, label in zip(axs.ravel(), ["$E_d$", "$L_w$", "$R_{rs}$"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[-1].set_xlabel("Wavelength [nm]")
axs[-1].set_xlim(300, 1300)

axs[0].set_ylim(ymin=0)
axs[1].set_ylim(ymin=0)
axs[2].set_ylim(ymin=0)

axs[0].set_title(f"All spectra ({N_total})")
plt.savefig(f"data/plots/spectra_all_data.pdf")
plt.show()
plt.close()

# Zoom
fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,7))

for data in data_all:
    wavelengths, Ed = split_spectrum(data, "Ed")
    wavelengths, Lw = split_spectrum(data, "Lw")
    wavelengths, R_rs = split_spectrum(data, "R_rs")

    for ax, spectrum in zip(axs.ravel(), [Ed, Lw, R_rs]):
        ax.plot(wavelengths, spectrum.T, c="k", alpha=0.01, zorder=1)

for ax, label in zip(axs.ravel(), ["$E_d$", "$L_w$", "$R_{rs}$"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[-1].set_xlabel("Wavelength [nm]")
axs[-1].set_xlim(300, 800)

axs[0].set_ylim(0, 2)
axs[1].set_ylim(0, 0.015)
axs[2].set_ylim(0, 0.01)

axs[0].set_title(f"All spectra ({N_total})")
plt.savefig(f"data/plots/spectra_all_data_zoom.pdf")
plt.show()
plt.close()

# Combined plots
fig = plt.figure(figsize=(7,3), tight_layout=True)
gs = fig.add_gridspec(2, 3)

ax_Ed = fig.add_subplot(gs[0,0])
ax_Lw = fig.add_subplot(gs[0,1])
ax_R_rs = fig.add_subplot(gs[0,2])
ax_N = fig.add_subplot(gs[1,2])
ax_zoom = fig.add_subplot(gs[1,:2])

axs = [ax_Ed, ax_Lw, ax_R_rs, ax_N, ax_zoom]
ylabels = ["$E_d$", "$L_w$", "$R_{rs}$", "$N$", "$R_{rs}$"]

for data in data_all:
    wavelengths, Ed = split_spectrum(data, "Ed")
    wavelengths, Lw = split_spectrum(data, "Lw")
    wavelengths, R_rs = split_spectrum(data, "R_rs")

    for ax, spectra in zip([ax_Ed, ax_Lw, ax_R_rs, ax_zoom], [Ed, Lw, R_rs, R_rs]):
        ax.plot(wavelengths, spectra.T, c="k", alpha=0.01)

ax_N.step(wavelength_range, wavelengths_bincount, where="mid", c="k")
#plt.xlim(wavelength_range[0], wavelength_range[-1]+5)

for ax, ylabel in zip(axs, ylabels):
    ax.grid(ls="--")
    ax.set_xlim(320, 1300)
    ax.set_ylim(ymin=0)
    ax.set_xlabel("Wavelength [nm]")
    ax.set_ylabel(ylabel)

ax_Ed.set_yticks([0,1,2,3])
ax_Lw.set_yticks([0,0.02,0.04,0.06,0.08])
ax_R_rs.set_yticks([0,0.02,0.04,0.06,0.08])
ax_N.set_yticks([0,500,1000,1500])
ax_zoom.set_yticks([0,0.002,0.004,0.006,0.008])

ax_zoom.set_xlim(320, 1000)
ax_zoom.set_ylim(0, 0.008)
ax_N.set_xlim(300, 1340)

plt.savefig("data/plots/spectra_all_data_combined.png", dpi=300)
plt.show()
plt.close()
