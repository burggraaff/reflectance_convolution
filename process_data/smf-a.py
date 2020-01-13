import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import get_keys_with_label, remove_negative_R_rs, remove_rows_based_on_threshold

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

Ed_keys, R_rs_keys = get_keys_with_label(data, "Ed", "R_rs")

for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    wavelength = int(Ed_k[3:])

    data[Ed_k].unit = u.microwatt / (u.cm**2 * u.nm)
    data[Ed_k] = data[Ed_k].to(u.watt / (u.m**2 * u.nm))

    data[R_rs_k].unit = 1 / u.steradian

    Lw = data[f"Ed_{wavelength}"] * data[f"R_rs_{wavelength}"]
    Lw.name = f"Lw_{wavelength}"
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

# Remove bad Ed row
remove_rows_based_on_threshold(data, "Ed", ">", 10)

remove_negative_R_rs(data)

map_data(data, data_label="SMF-A", projection='gnom', lat_0=30, lon_0=-83, llcrnrlon=-89, urcrnrlon=-77, llcrnrlat=24, urcrnrlat=32, resolution="h", parallels=np.arange(24, 36, 2), meridians=np.arange(-90, -70, 2))

plot_spectra(data, data_label="SMF-A", alpha=0.2)

write_data(data, label="SMF-A")
