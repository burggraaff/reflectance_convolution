import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import remove_negative_R_rs, convert_to_unit

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

convert_to_unit(data, "Ed", u.microwatt / (u.centimeter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "Lw", u.microwatt / (u.centimeter**2 * u.nanometer * u.steradian), u.watt / (u.meter**2 * u.nanometer * u.steradian))
convert_to_unit(data, "R_rs", 1 / u.steradian)

remove_negative_R_rs(data)

map_data(data, data_label="AS11", projection="gnom", lat_0=21, lon_0=68, llcrnrlon=65, urcrnrlon=71, llcrnrlat=18, urcrnrlat=24, resolution="h", parallels=np.arange(16, 30, 2), meridians=np.arange(66, 72, 2))

plot_spectra(data, data_label="AS11", alpha=0.7)

write_data(data, label="AS11")
