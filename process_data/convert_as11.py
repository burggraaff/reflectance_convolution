import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra

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

Ed_keys = [key for key in data.keys() if "Ed" in key]
Lw_keys = [key for key in data.keys() if "Lw" in key]
R_rs_keys = [key for key in data.keys() if "Rrs" in key]

for key in Ed_keys:
    data[key].unit = u.microwatt / (u.cm**2 * u.nm)
    data[key] = data[key].to(u.watt / (u.m**2 * u.nm))

for key in Lw_keys:
    data[key].unit = u.microwatt / (u.cm**2 * u.nm * u.steradian)
    data[key] = data[key].to(u.watt / (u.m**2 * u.nm * u.steradian))

for key in R_rs_keys:
    data[key].unit = 1 / u.steradian

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

plot_spectra(data, data_label="AS11", alpha=0.7)

data.write("data/as11_processed.tab", format="ascii.fast_tab", overwrite=True)