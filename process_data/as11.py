import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import get_keys_with_label

folder = Path("data/AS11/")
files = list(folder.glob("AS*HTSRB.csv"))

tabs = []
for file in files:
    wavelengths, Lw, Es, Rrs = np.loadtxt(file, delimiter=",", skiprows=43, unpack=True, usecols=[0,1,3,4])

    date, time, lon, lat = find_auxiliary_information_seabass(file)

    cols = ["Date", "Time", "Latitude", "Longitude"] + [f"Lw_{wvl:.2f}" for wvl in wavelengths] + [f"Ed_{wvl:.2f}" for wvl in wavelengths] + [f"R_rs_{wvl:.2f}" for wvl in wavelengths]
    dtype = [int, "S8", float, float] + 3 * [float for wvl in wavelengths]
    tab = table.Table(rows=[[date, time, lat, lon, *Lw, *Es, *Rrs]], names=cols, dtype=dtype)
    tabs.append(tab)

data = table.vstack(tabs)

print(f"Original N = {len(data)}")

Ed_keys, Lw_keys, R_rs_keys = get_keys_with_label(data, "Ed", "Lw", "R_rs")

for key in Ed_keys:
    data[key].unit = u.microwatt / (u.cm**2 * u.nm)
    data[key] = data[key].to(u.watt / (u.m**2 * u.nm))

for key in Lw_keys:
    data[key].unit = u.microwatt / (u.cm**2 * u.nm * u.steradian)
    data[key] = data[key].to(u.watt / (u.m**2 * u.nm * u.steradian))

for key in R_rs_keys:
    data[key].unit = 1 / u.steradian

remove_indices = [i for i, row in enumerate(data) if any(row[key] < 0 for key in R_rs_keys)]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with R_rs < 0")

map_data(data, data_label="AS11", projection="gnom", lat_0=21, lon_0=68, llcrnrlon=65, urcrnrlon=71, llcrnrlat=18, urcrnrlat=24, resolution="h", parallels=np.arange(16, 30, 2), meridians=np.arange(66, 72, 2))

plot_spectra(data, data_label="AS11", alpha=0.7)

write_data(data, label="AS11")
