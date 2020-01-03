import numpy as np
from bandaveraging import split_spectrum, load_data_full
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap

data = load_data_full()

# Plot map of observations
fig = plt.figure(figsize=(10, 7.5), tight_layout=True)

m = Basemap(projection='moll', lon_0=0, resolution="c")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(-90, 95, 15), labels=[1,1,0,0])
m.drawmeridians(np.arange(-180, 180, 30), labels=[0,0,1,1])

m.scatter(data["Longitude"], data["Latitude"], latlon=True, c="r", edgecolors="r", s=25)

plt.savefig("combined_map.pdf")
plt.show()

# Plot all Ed, Lu, Ls, R_rs spectra
fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0})

wavelengths, Ed = split_spectrum(data, "Ed")
wavelengths, Lu = split_spectrum(data, "Lu")
wavelengths, Ls = split_spectrum(data, "Ls")
wavelengths, Lw = split_spectrum(data, "Lw")
wavelengths, R_rs = split_spectrum(data, "R_rs")

for ax, spec in zip(axs, [Ed, Lu, Ls, Lw, R_rs]):
    ax.plot(wavelengths, spec.T, c="k", alpha=0.15, zorder=1)

for ax, label in zip(axs, ["$E_d$ [W m$^{-2}$ nm$^{-1}$]", "$L_u$ [W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$L_s$ [W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$L_w$ [W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$R_{rs}$ [sr$^{-1}$]"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[4].set_xlabel("Wavelength [nm]")

fig.suptitle("All spectra")
plt.savefig("combined_spectra.pdf")
plt.show()
plt.close()

axs[0,0].tick_params(bottom=False, labelbottom=False)
axs[0,1].tick_params(bottom=False, labelbottom=False)
axs[0,1].tick_params(left=False, labelleft=False, right=True, labelright=True)
axs[1,1].tick_params(left=False, labelleft=False, right=True, labelright=True)
axs[0,1].yaxis.set_label_position("right")
axs[1,1].yaxis.set_label_position("right")
axs[1,0].set_xlabel("Wavelength [nm]")
axs[1,1].set_xlabel("Wavelength [nm]")
axs[0,0].set_xlim(320, 950)

