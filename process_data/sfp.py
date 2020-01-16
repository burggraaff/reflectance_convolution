import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import get_keys_with_label, convert_to_unit

folder = Path("data/SFP/")
files = sorted(folder.glob("*.txt"))

tabs = []
for file in files:
    for skiprows in range(44, 52):
        try:
            wavelengths, Ed, Rrs = np.loadtxt(file, skiprows=skiprows, unpack=True, usecols=[0,3,4])
        except Exception:
            continue
        else:
            break

    date, time, lon, lat = find_auxiliary_information_seabass(file)

    cols = ["Date", "Time", "Latitude", "Longitude"] + [f"Ed_{wvl:.2f}" for wvl in wavelengths] + [f"R_rs_{wvl:.2f}" for wvl in wavelengths]
    dtype = [int, "S8", float, float] + 2 * [float for wvl in wavelengths]
    tab = table.Table(rows=[[date, time, lat, lon, *Ed, *Rrs]], names=cols, dtype=dtype)
    tabs.append(tab)

data = table.vstack(tabs)

convert_to_unit(data, "Ed", u.microwatt / (u.centimeter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "R_rs", 1 / u.steradian)

Ed_keys, R_rs_keys = get_keys_with_label(data, "Ed", "R_rs")
for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    Lw_k = Ed_k.replace("Ed", "Lw")
    Lw = data[Ed_k] * data[R_rs_k]
    Lw.name = Lw_k
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

map_data(data, data_label="SFP", projection='gnom', lat_0=25, lon_0=-81, llcrnrlon=-87, urcrnrlon=-75, llcrnrlat=20, urcrnrlat=29, resolution="h", parallels=np.arange(20, 35, 2), meridians=np.arange(-90, -70, 2))

plot_spectra(data, data_label="SFP", alpha=0.1)

write_data(data, label="SFP")
