import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import get_keys_with_label, remove_negative_R_rs

folder = Path("data/GasEx/")
files = list(folder.glob("AOP*.txt"))

tabs = []
for file in files:
    try:
        wavelengths, Es, Rrs = np.loadtxt(file, delimiter="\t", skiprows=40, unpack=True, usecols=[0,1,5])
    except:
        wavelengths, Es, Rrs = np.loadtxt(file, delimiter="\t", skiprows=41, unpack=True, usecols=[0,1,5])

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

for wavelength in np.arange(712, 722, 2, dtype=int):
    data.remove_columns(get_keys_with_label(data, f"_{wavelength}"))

map_data(data, data_label="GasEx", projection='gnom', lat_0=36.9, lon_0=-75.8, llcrnrlon=-80, urcrnrlon=-70, llcrnrlat=32, urcrnrlat=42, resolution="h", parallels=np.arange(30, 45, 2), meridians=np.arange(-80, -70, 2))

plot_spectra(data, data_label="GasEx", alpha=0.2)

write_data(data, label="GasEx")