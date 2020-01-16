import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import get_keys_with_label

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

Es_keys, R_rs_keys = get_keys_with_label(data, "ES", "Rrs")

for Es_k, R_rs_k in zip(Es_keys, R_rs_keys):
    wavelength = float(Es_k[2:])

    data[Es_k].unit = u.microwatt / (u.cm**2 * u.nm)
    data[Es_k] = data[Es_k].to(u.watt / (u.m**2 * u.nm))
    data.rename_column(Es_k, f"Ed_{wavelength:.1f}")

    data[R_rs_k].unit = 1 / u.steradian
    data.rename_column(R_rs_k, f"R_rs_{wavelength:.1f}")

    Lw = data[f"Ed_{wavelength:.1f}"] * data[f"R_rs_{wavelength:.1f}"]
    Lw.name = f"Lw_{wavelength:.1f}"
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

map_data(data, data_label="TaraO", lon_0=0, resolution="i")

plot_spectra(data, data_label="TaraO", alpha=0.1)

write_data(data, label="TaraO")
