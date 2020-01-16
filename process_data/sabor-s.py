import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import get_keys_with_label, remove_negative_R_rs, convert_to_unit, rename_columns
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

rename_columns(data, "Es", "Ed_")
rename_columns(data, "Rrs", "R_rs_")

convert_to_unit(data, "Ed", u.microwatt / (u.centimeter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "R_rs", 1 / u.steradian)

Ed_keys, R_rs_keys = get_keys_with_label(data, "Ed", "R_rs")
for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    Lw = data[Ed_k] * data[R_rs_k]
    Lw.name = Ed_k.replace("Ed", "Lw")
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

# Remove NaN (-999) values and columns with negative values (wavelength crop)
remove_wavelengths = [*np.arange(350, 358), *np.arange(750, 801)]
for wvl in remove_wavelengths:
    remove_keys = get_keys_with_label(data, f"{wvl}")
    data.remove_columns(remove_keys)

remove_negative_R_rs(data)

map_data(data, data_label="SABOR-S", projection="gnom", lat_0=37, lon_0=-70, llcrnrlon=-77, urcrnrlon=-64, llcrnrlat=35, urcrnrlat=43, resolution="h", parallels=np.arange(32, 45, 2), meridians=np.arange(-80, -60, 2))

plot_spectra(data, data_label="SABOR-S", alpha=0.5)

write_data(data, label="SABOR-S")
