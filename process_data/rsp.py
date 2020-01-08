import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data

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

    data.remove_columns([key for key in data.keys() if "sd" in key])
    data.remove_columns([key for key in data.keys() if "Lu" in key])

    data_tables.append(data)

data = table.vstack(data_tables)

Ed_keys = [key for key in data.keys() if "Ed" in key]
R_rs_keys = [key for key in data.keys() if "Rrs" in key]

for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    wavelength = float(Ed_k[2:])

    data[Ed_k].unit = u.microwatt / (u.cm**2 * u.nm)
    data[Ed_k] = data[Ed_k].to(u.watt / (u.m**2 * u.nm))
    data.rename_column(Ed_k, f"Ed_{wavelength:.1f}")

    data[R_rs_k].unit = 1 / u.steradian
    data.rename_column(R_rs_k, f"R_rs_{wavelength:.1f}")

    Lw = data[f"Ed_{wavelength:.1f}"] * data[f"R_rs_{wavelength:.1f}"]
    Lw.name = f"Lw_{wavelength:.1f}"
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

map_data(data, data_label="RSP", lon_0=0, resolution="i")

plot_spectra(data, data_label="RSP", alpha=0.3)

write_data(data, label="RSP")
