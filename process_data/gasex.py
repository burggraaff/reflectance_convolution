import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import get_keys_with_label, convert_to_unit, add_Lw_from_Ed_Rrs

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

convert_to_unit(data, "Ed", u.microwatt / (u.centimeter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "R_rs", 1 / u.steradian)

data = add_Lw_from_Ed_Rrs(data)

for wavelength in np.arange(712, 722, 2, dtype=int):
    data.remove_columns(get_keys_with_label(data, f"_{wavelength}"))

map_data(data, data_label="GasEx", projection='gnom', lat_0=-52, lon_0=-38, llcrnrlon=-60, urcrnrlon=-30, llcrnrlat=-60, urcrnrlat=-35, resolution="h", parallels=np.arange(-60, -20, 5), meridians=np.arange(-60, -10, 5))

plot_spectra(data, data_label="GasEx", alpha=0.2)

write_data(data, label="GasEx")
