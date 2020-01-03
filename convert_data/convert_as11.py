import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from pathlib import Path

folder = Path("data/AS11/")
files = list(folder.glob("AS*HTSRB.csv"))

tabs = []
for file in files:
    wavelengths, Lw, Es, Rrs = np.loadtxt(file, delimiter=",", skiprows=43, unpack=True, usecols=[0,1,3,4])

    # Reject data with negative Rrs
    if len(np.where(Rrs < 0)[0]) > 0:
        continue

    with open(file, "r") as f:
        lines = f.readlines()
        lat, lon = lines[20:22]
        lat = float(lat[16:-6])
        lon = float(lon[16:-6])

        date, time = lines[16:18]
        date = int(date[10:-1])
        time = time[12:-6]

    cols = ["Date", "Time", "Latitude", "Longitude"] + [f"Lw_{wvl:.2f}" for wvl in wavelengths] + [f"Ed_{wvl:.2f}" for wvl in wavelengths] + [f"R_rs_{wvl:.2f}" for wvl in wavelengths]
    dtype = [int, "S8", float, float] + 3 * [float for wvl in wavelengths]
    tab = table.Table(rows=[[date, time, lat, lon, *Lw, *Es, *Rrs]], names=cols, dtype=dtype)
    tabs.append(tab)

data = table.vstack(tabs)

# Plot map of observations
fig = plt.figure(figsize=(10, 6), tight_layout=True)

m = Basemap(projection='gnom', lat_0=21, lon_0=68, llcrnrlon=65, urcrnrlon=71, llcrnrlat=18, urcrnrlat=24, resolution="h")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(16, 30, 2), labels=[1,1,0,0])
m.drawmeridians(np.arange(60, 75, 2), labels=[0,0,1,1])

m.scatter(data["Longitude"], data["Latitude"], latlon=True, c="r", edgecolors="k", s=60, zorder=10)

plt.savefig("data/plots/map_AS11.pdf")
plt.show()

# Plot all Es, Lw, R_rs spectra
wavelengths = [float(key[3:]) for key in data.keys() if "Lw" in key]

fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,7))

for row in data:
    spec_Ed = [row[key] for key in data.keys() if "Ed" in key]
    spec_Lw = [row[key] for key in data.keys() if "Lw" in key]
    spec_R_rs = [row[key] for key in data.keys() if "R_rs" in key]

    for ax, spec in zip(axs.ravel(), [spec_Ed, spec_Lw, spec_R_rs]):
        ax.plot(wavelengths, spec, c="k", alpha=0.7, zorder=1)

for ax, label in zip(axs.ravel(), ["$E_d$ [$\mu$W cm$^{-2}$ nm$^{-1}$]", "$L_w$ [$\mu$W cm$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$R_{rs}$ [sr$^{-1}$]"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[-1].set_xlabel("Wavelength [nm]")
axs[-1].set_xlim(400, 750)

axs[0].set_title(f"AS11 spectra ({len(data)})")
plt.savefig("data/plots/spectra_AS11.pdf")
plt.show()
plt.close()

data.write("data/as11_processed.tab", format="ascii.fast_tab", overwrite=True)
