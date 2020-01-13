import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import get_keys_with_label, split_spectrum, remove_rows_based_on_threshold

data = read("data/TAOM/ep1_hr3.avg.prod_1_1001.ftp", data_start=30)
header = read("data/TAOM/ep1_hr3.avg.prod_1_1001.ftp", data_start=27, data_end=28)
header["col1"][0] = "year"
header = header[0].as_void()

for key, new_key in zip(data.keys(), header):
    data.rename_column(key, new_key)

date, time, lon, lat = find_auxiliary_information_seabass("data/TAOM/ep1_hr3.avg.prod_1_1001.ftp")
data.add_column(table.Column(name="Latitude", data=[lat]*len(data)))
data.add_column(table.Column(name="Longitude", data=[lon]*len(data)))

data.remove_columns(get_keys_with_label(data, "Lwn"))

Lw_keys, R_rs_keys = get_keys_with_label(data, "Lw", "Rrs")

for Lw_k, R_rs_k in zip(Lw_keys, R_rs_keys):
    wavelength = float(Lw_k[2:])

    data[Lw_k].unit = u.microwatt / (u.cm**2 * u.nm * u.steradian)
    data[Lw_k] = data[Lw_k].to(u.watt / (u.m**2 * u.nm * u.steradian))
    data.rename_column(Lw_k, f"Lw_{wavelength:.4f}")

    data[R_rs_k].unit = 1 / u.steradian
    data.rename_column(R_rs_k, f"R_rs_{wavelength:.4f}")

    Ed = data[f"Lw_{wavelength:.4f}"] / data[f"R_rs_{wavelength:.4f}"]
    Ed.name = f"Ed_{wavelength:.4f}"
    Ed.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Ed)

# Remove columns commonly missing data
wavelengths, Ed = split_spectrum(data, "Ed")
for wavelength in wavelengths[np.where((wavelengths < 371) | (wavelengths > 734))]:
    data.remove_columns(get_keys_with_label(data, f"_{wavelength:.4f}"))

# Remove columns with noisy data
wavelengths, Ed = split_spectrum(data, "Ed")
for wavelength in wavelengths[np.where(wavelengths > 685)]:
    data.remove_columns(get_keys_with_label(data, f"_{wavelength:.4f}"))

# Remove rows that still miss data
remove_rows_based_on_threshold(data, "R_rs", "<", -900)

# Remove outliers with R_rs > 0.02
remove_rows_based_on_threshold(data, "R_rs", ">", 0.02)

map_data(data, data_label="TAOM", projection="gnom", lat_0=0, lon_0=-155, llcrnrlon=-185, urcrnrlon=-125, llcrnrlat=-30, urcrnrlat=30, resolution="h", parallels=np.arange(-30, 40, 10), meridians=np.arange(-190, -120, 10))

plot_spectra(data, data_label="TAOM", alpha=0.5)

write_data(data, label="TAOM")
