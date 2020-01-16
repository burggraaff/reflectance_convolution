import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import get_keys_with_label, convert_to_unit, rename_columns

folder = Path("data/TaraO/")
files = list(folder.glob("Tara_HyperPro*.txt"))

data = table.vstack([read(file, data_start=35) for file in files])

header = read(files[0], data_start=32, data_end=33)
header["col1"][0] = "year"
header = header[0].as_void()

for key, new_key in zip(data.keys(), header):
    data.rename_column(key, new_key)

data.remove_columns(get_keys_with_label(data, "LU"))

data.rename_column("lat", "Latitude")
data.rename_column("lon", "Longitude")

rename_columns(data, "ES", "Ed_", exclude="None")
rename_columns(data, "Rrs", "R_rs_", exclude="None")

convert_to_unit(data, "Ed", u.microwatt / (u.centimeter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "R_rs", 1 / u.steradian)

Ed_keys, R_rs_keys = get_keys_with_label(data, "Ed", "R_rs")
for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    wavelength = float(Ed_k[3:])
    Ed_sd = Ed_k + "_sd"
    R_rs_sd = R_rs_k + "_sd"

    Lw = data[Ed_k] * data[R_rs_k]
    Lw_key = f"Lw_{wavelength:.1f}"
    Lw.name = Lw_key
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

map_data(data, data_label="TaraO", lon_0=0, resolution="i")

plot_spectra(data, data_label="TaraO", alpha=0.1)

write_data(data, label="TaraO")
