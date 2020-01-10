import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import get_keys_with_label
import csv

csv.field_size_limit(1000000)  # Increase to allow large number of columns

folder = Path("data/SABOR/")
files = list(folder.glob("CCNY*.sb"))

data_tables = []
for file in files:
    data_table = read(file, data_start=35)
    header = read(files[0], data_start=32, data_end=33)
    header["col1"][0] = "year"
    header = header[0].as_void()

    for key, new_key in zip(data_table.keys(), header):
        data_table.rename_column(key, new_key)

    data_table.remove_columns(get_keys_with_label(data_table, "stokes"))
    data_table.remove_columns(get_keys_with_label(data_table, "sd"))
    data_table.remove_columns(get_keys_with_label(data_table, "sky"))
    data_table.remove_columns(get_keys_with_label(data_table, "Lt"))
    # These data are normalised to R_rs(750 nm), so we must re-calculate Lw to get a fair comparison
    data_table.remove_columns(get_keys_with_label(data_table, "Lw"))
    data_table.remove_columns(get_keys_with_label(data_table, "AOT"))

    data_tables.append(data_table)
    print(file)

data = table.vstack(data_tables)

data.remove_column("col15346")
data.rename_column("lat", "Latitude")
data.rename_column("lon", "Longitude")

Es_keys, R_rs_keys = get_keys_with_label(data, "Es", "Rrs")

for Es_k, R_rs_k in zip(Es_keys, R_rs_keys):
    wavelength = int(Es_k[2:])

    data[Es_k].unit = u.microwatt / (u.cm**2 * u.nm)
    data[Es_k] = data[Es_k].to(u.watt / (u.m**2 * u.nm))
    data.rename_column(Es_k, f"Ed_{wavelength}")

    data[R_rs_k].unit = 1 / u.steradian
    data.rename_column(R_rs_k, f"R_rs_{wavelength}")

    Lw = data[f"Ed_{wavelength}"] * data[f"R_rs_{wavelength}"]
    Lw.name = f"Lw_{wavelength}"
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

# Remove NaN (-999) values and columns with negative values (wavelength crop)
remove_wavelengths = [*np.arange(350, 358), *np.arange(750, 801)]
for wvl in remove_wavelengths:
    try:
        data.remove_columns([f"{s}_{wvl}" for s in ["Ed", "Lw", "R_rs"]])
    except KeyError:
        continue

# Clip small negative values (0 > R_rs > -1e-4) to 0
Lw_keys, R_rs_keys = get_keys_with_label(data, "Lw", "R_rs")
for Lw_k, R_rs_k in zip(Lw_keys, R_rs_keys):
    ind = np.where((data[R_rs_k] < 0) & (data[R_rs_k] > -1e-4))
    data[R_rs_k][ind] = 0
    data[Lw_k][ind] = 0

map_data(data, data_label="SABOR-S", projection="gnom", lat_0=37, lon_0=-70, llcrnrlon=-77, urcrnrlon=-64, llcrnrlat=35, urcrnrlat=43, resolution="h", parallels=np.arange(32, 45, 2), meridians=np.arange(-80, -60, 2))

plot_spectra(data, data_label="SABOR-S", alpha=0.5)

write_data(data, label="SABOR-S")
