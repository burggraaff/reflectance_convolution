import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import remove_negative_R_rs, remove_rows_based_on_threshold, convert_to_unit, add_Lw_from_Ed_Rrs

folder = Path("data/SMF/ASD/")
files = list(folder.glob("*ASD*"))

tabs = []
for file in files:
    wavelengths, Es, Rrs = np.loadtxt(file, delimiter="\t", skiprows=50, unpack=True, usecols=[0,1,3])

    date, time, lon, lat = find_auxiliary_information_seabass(file)

    cols = ["Date", "Time", "Latitude", "Longitude"] + [f"Ed_{wvl:.0f}" for wvl in wavelengths] + [f"R_rs_{wvl:.0f}" for wvl in wavelengths]
    dtype = [int, "S8", float, float] + 2 * [float for wvl in wavelengths]
    tab = table.Table(rows=[[date, time, lat, lon, *Es, *Rrs]], names=cols, dtype=dtype)
    tabs.append(tab)

data = table.vstack(tabs)

convert_to_unit(data, "Ed", u.microwatt / (u.centimeter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "R_rs", 1 / u.steradian)

data = add_Lw_from_Ed_Rrs(data)

# Remove bad Ed row
remove_rows_based_on_threshold(data, "Ed", ">", 10)

remove_negative_R_rs(data)

map_data(data, data_label="SMF-A", projection='gnom', lat_0=30, lon_0=-83, llcrnrlon=-89, urcrnrlon=-77, llcrnrlat=24, urcrnrlat=32, resolution="h", parallels=np.arange(24, 36, 2), meridians=np.arange(-90, -70, 2))

plot_spectra(data, data_label="SMF-A", alpha=0.2)

write_data(data, label="SMF-A")
