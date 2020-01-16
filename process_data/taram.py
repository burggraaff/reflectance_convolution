import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import get_keys_with_label, convert_to_unit, rename_columns, add_Lw_from_Ed_Rrs

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

rename_columns(data, "ES", "Ed_", exclude="None")
rename_columns(data, "Rrs", "R_rs_", exclude="None")

convert_to_unit(data, "Ed", u.microwatt / (u.centimeter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "R_rs", 1 / u.steradian)

data = add_Lw_from_Ed_Rrs(data)

map_data(data, data_label="TaraM", projection="gnom", lat_0=43.5, lon_0=19, llcrnrlon=0, urcrnrlon=40, llcrnrlat=30, urcrnrlat=44, resolution="h", parallels=np.arange(30, 50, 2.5), meridians=np.arange(0, 45, 5))

plot_spectra(data, data_label="TaraM", alpha=0.1)

write_data(data, label="TaraM")
