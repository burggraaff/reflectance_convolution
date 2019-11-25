import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from pathlib import Path

folder = Path("data/Tara-O/")
files = list(folder.glob("Tara_HyperPro*.txt"))

data = table.vstack([read(file, data_start=35) for file in files])

header = read(files[0], data_start=32, data_end=33)
header["col1"][0] = "year"
header = header[0].as_void()

for key, new_key in zip(data.keys(), header):
    data.rename_column(key, new_key)

data.remove_columns([key for key in data.keys() if "sd" in key])
data.remove_columns([key for key in data.keys() if "LU" in key])

Es_keys = [key for key in data.keys() if "ES" in key]
R_rs_keys = [key for key in data.keys() if "Rrs" in key]

for Es_k, R_rs_k in zip(Es_keys, R_rs_keys):
    wavelength = float(Es_k[2:])
    Es = data[Es_k]
    R_rs = data[R_rs_k]
    Lw = R_rs * Es

    data[Es_k].unit = u.microwatt / (u.cm**2 * u.nm)
    data.rename_column(Es_k, f"Ed_{wavelength:.1f}")

    data[R_rs_k].unit = 1 / u.steradian
    data.rename_column(R_rs_k, f"R_rs_{wavelength:.1f}")

    Lw.name = f"Lw_{wavelength:.1f}"
    Lw.unit = u.microwatt / (u.cm**2 * u.nm * u.steradian)
    data.add_column(Lw)


# Plot map of observations
fig = plt.figure(figsize=(10, 6), tight_layout=True)

m = Basemap(projection='moll', lon_0=0, resolution="i")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(-90, 95, 15), labels=[1,1,0,0])
m.drawmeridians(np.arange(-180, 180, 30), labels=[0,0,1,1])

m.scatter(data["lon"], data["lat"], latlon=True, c="r", edgecolors="k", s=60, zorder=10)

plt.savefig("map_Tara-O.pdf")
plt.show()

# Plot all Es, Lw, R_rs spectra
wavelengths = [float(key[3:]) for key in data.keys() if "Lw" in key]

fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,7))

for row in data:
    spec_Es = [row[key] for key in data.keys() if "Ed" in key]
    spec_Lw = [row[key] for key in data.keys() if "Lw" in key]
    spec_R_rs = [row[key] for key in data.keys() if "R_rs" in key]

    for ax, spec in zip(axs.ravel(), [spec_Es, spec_Lw, spec_R_rs]):
        ax.plot(wavelengths, spec, c="k", alpha=0.1, zorder=1)

for ax, label in zip(axs.ravel(), ["$E_s$ [$\mu$W cm$^{-2}$ nm$^{-1}$]", "$L_w$ [$\mu$W cm$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$R_{rs}$ [sr$^{-1}$]"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[-1].set_xlabel("Wavelength [nm]")
axs[-1].set_xlim(400, 750)

axs[0].set_title("Tara-O spectra")
plt.savefig("spectra_Tara-O.pdf")
plt.show()
plt.close()

data.write("data/tara-o_processed.tab", format="ascii.fast_tab", overwrite=True)
