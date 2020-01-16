import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import get_keys_with_label, convert_to_unit

folder = Path("data/TaraM/")
files = list(folder.glob("Tara_HyperPro*.txt"))

data = table.vstack([read(file, data_start=35) for file in files])

header = read(files[0], data_start=32, data_end=33)
header["col1"][0] = "year"
header = header[0].as_void()

for key, new_key in zip(data.keys(), header):
    data.rename_column(key, new_key)

data.remove_columns(get_keys_with_label(data, "sd"))
data.remove_columns(get_keys_with_label(data, "LU"))

data.rename_column("lat", "Latitude")
data.rename_column("lon", "Longitude")

Es_keys, R_rs_keys = get_keys_with_label(data, "ES", "Rrs")

for Es_k, R_rs_k in zip(Es_keys, R_rs_keys):
    wavelength = float(Es_k[2:])

    convert_to_unit(data, Es_k, u.microwatt / (u.centimeter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
    convert_to_unit(data, R_rs_k, 1 / u.steradian)

    data.rename_column(Es_k, f"Ed_{wavelength:.1f}")
    data.rename_column(R_rs_k, f"R_rs_{wavelength:.1f}")

    Lw = data[f"Ed_{wavelength:.1f}"] * data[f"R_rs_{wavelength:.1f}"]
    Lw.name = f"Lw_{wavelength:.1f}"
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

map_data(data, data_label="TaraM", projection="gnom", lat_0=43.5, lon_0=19, llcrnrlon=0, urcrnrlon=40, llcrnrlat=30, urcrnrlat=44, resolution="h", parallels=np.arange(30, 50, 2.5), meridians=np.arange(0, 45, 5))

plot_spectra(data, data_label="TaraM", alpha=0.1)

write_data(data, label="TaraM")
