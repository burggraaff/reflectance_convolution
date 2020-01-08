import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data

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

map_data(data, data_label="AS11", projection="gnom", lat_0=21, lon_0=68, llcrnrlon=65, urcrnrlon=71, llcrnrlat=18, urcrnrlat=24, resolution="h", parallels=np.arange(16, 30, 2), meridians=np.arange(66, 72, 2))

plot_spectra(data, data_label="AS11", alpha=0.7)

write_data(data, label="AS11")
