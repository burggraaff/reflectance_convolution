import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import get_keys_with_label, convert_to_unit, rename_columns, add_Lw_from_Ed_Rrs

folder = Path("data/RSP/")
files = list(folder.glob("*.txt"))

data_tables = []

for file in files:
    data = read(file, data_start=37)
    header = read(file, data_start=32, data_end=33)
    header["col1"][0] = "date"
    header = header[0].as_void()
    for key, new_key in zip(data.keys(), header):
        data.rename_column(key, new_key)

    data.rename_column("lat", "Latitude")
    data.rename_column("lon", "Longitude")

    data.remove_columns(get_keys_with_label(data, "sd"))
    data.remove_columns(get_keys_with_label(data, "Lu"))

    data_tables.append(data)

data = table.vstack(data_tables)

rename_columns(data, "Ed", "Ed_")
rename_columns(data, "Rrs", "R_rs_")

convert_to_unit(data, "Ed", u.microwatt / (u.centimeter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "R_rs", 1 / u.steradian)

data = add_Lw_from_Ed_Rrs(data)

map_data(data, data_label="RSP", lon_0=0, resolution="i")

plot_spectra(data, data_label="RSP", alpha=0.3)

write_data(data, label="RSP")
