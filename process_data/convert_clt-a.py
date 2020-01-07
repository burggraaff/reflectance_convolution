import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra

folder = Path("data/CLT/")
files = list(folder.glob("*ASD*.txt"))

tabs = []
for file in files:
    skiprows = 51 if "2007" in file.stem else 47
    wavelengths, Es, Rrs = np.loadtxt(file, delimiter="\t", skiprows=skiprows, unpack=True, usecols=[0,1,3])

    # Reject data with negative Rrs
    if len(np.where(Rrs < 0)[0]) > 0:
        continue

    with open(file, "r") as f:
        lines = f.readlines()
        lat, lon = lines[14:16]
        lat = float(lat[16:-6])
        lon = float(lon[16:-6])

        date, time = lines[18:20]
        date = int(date[10:])
        time = time[12:-6]

    cols = ["Date", "Time", "Latitude", "Longitude"] + [f"Ed_{wvl:.0f}" for wvl in wavelengths] + [f"R_rs_{wvl:.0f}" for wvl in wavelengths]
    dtype = [int, "S8", float, float] + 2 * [float for wvl in wavelengths]
    tab = table.Table(rows=[[date, time, lat, lon, *Es, *Rrs]], names=cols, dtype=dtype)
    tabs.append(tab)

data = table.vstack(tabs)

Ed_keys = [key for key in data.keys() if "Ed" in key]
R_rs_keys = [key for key in data.keys() if "R_rs" in key]

for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    wavelength = int(Ed_k[3:])

    data[Ed_k].unit = u.microwatt / (u.cm**2 * u.nm)
    data[Ed_k] = data[Ed_k].to(u.watt / (u.m**2 * u.nm))

    data[R_rs_k].unit = 1 / u.steradian

    Lw = data[f"Ed_{wavelength}"] * data[f"R_rs_{wavelength}"]
    Lw.name = f"Lw_{wavelength}"
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

# Plot map of observations
fig = plt.figure(figsize=(10, 6), tight_layout=True)

m = Basemap(projection='gnom', lat_0=36.9, lon_0=-75.8, llcrnrlon=-80, urcrnrlon=-70, llcrnrlat=32, urcrnrlat=42, resolution="h")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(30, 45, 2), labels=[1,1,0,0])
m.drawmeridians(np.arange(-80, -70, 2), labels=[0,0,1,1])

m.scatter(data["Longitude"], data["Latitude"], latlon=True, c="r", edgecolors="k", s=60, zorder=10)

plt.savefig("data/plots/map_CLT-A.pdf")
plt.show()

plot_spectra(data, data_label="CLT-A", alpha=0.2)

data.write("data/clt-a_processed.tab", format="ascii.fast_tab", overwrite=True)
