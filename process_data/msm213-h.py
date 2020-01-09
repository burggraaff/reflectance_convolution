import numpy as np
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.bandaveraging import split_spectrum

Lw = read("data/MSM21_3/MSM21_3_Lw-5nm.tab", data_start=132, header_start=131)
Rrs = read("data/MSM21_3/MSM21_3_Rrs-5nm.tab", data_start=133, header_start=132)

wavelengths = np.arange(360, 805, 5)
for wvl in wavelengths:
    try:  # mu gets properly loaded on Linux
        Lw.rename_column(f"Lw_{wvl} [µW/cm**2/nm/sr]", f"Lw_{wvl}")
    except KeyError: # but not on Windows
        Lw.rename_column(f"Lw_{wvl} [ÂµW/cm**2/nm/sr]", f"Lw_{wvl}")
    Lw[f"Lw_{wvl}"].unit = u.microwatt / (u.centimeter**2 * u.nanometer * u.steradian)
    Lw[f"Lw_{wvl}"] = Lw[f"Lw_{wvl}"].to(u.watt / (u.meter**2 * u.nanometer * u.steradian))

    Rrs.rename_column(f"Rrs_{wvl} [1/sr]", f"R_rs_{wvl}")
    Rrs[f"R_rs_{wvl}"].unit = 1 / u.steradian

data = table.join(Lw, Rrs, keys=["Date/Time"])
data.rename_column("Event_1", "Event")
data.rename_column("Latitude_1", "Latitude")
data.rename_column("Longitude_1", "Longitude")

for wvl in wavelengths:
    Ed = data[f"Lw_{wvl}"] / data[f"R_rs_{wvl}"]
    Ed.name = f"Ed_{wvl}"
    Ed.unit = u.watt / (u.meter**2 * u.nanometer)
    data.add_column(Ed)

# Remove rows with NaN values
R_rs_keys = [key for key in data.keys() if "R_rs" in key]
remove_indices = [i for i, row_mask in enumerate(data.mask) if any(row_mask[key] for key in R_rs_keys)]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with NaN values")

# Remove rows with consecutive jumps in Ed >= threshold between wavelengths
threshold = 0.2
wavelengths, Ed = split_spectrum(data, "Ed")
diffs = np.diff(Ed, axis=1)
diffs_abs = np.abs(diffs)
all_rows = [np.where((col1 >= threshold) & (col2 >= threshold))[0] for col1, col2 in zip(diffs_abs.T, diffs_abs.T[1:])]
all_rows = [row for sub in all_rows for row in sub]
remove_indices = np.unique(all_rows)
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with consecutive Ed jumps >= {threshold}")

# Remove rows with singular jumps in Ed > 0.4 between wavelengths
wavelengths, Ed = split_spectrum(data, "Ed")
diffs = np.diff(Ed, axis=1)
remove_indices = np.unique(np.where(np.abs(diffs) >= 0.35)[0])
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with Ed jumps >= 0.35")

# Remove rows with negative R_rs
remove_indices = [i for i, row in enumerate(data) if any(row[key] <= -0.001 for key in R_rs_keys)]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with negative values")

map_data(data, data_label="MSM213-H", projection='gnom', lat_0=66, lon_0=-40.5, llcrnrlon=-53, urcrnrlon=-12, llcrnrlat=58, urcrnrlat=70.5, resolution="h", parallels=np.arange(55, 75, 5), meridians=np.arange(-60, -5, 5))

plot_spectra(data, data_label="MSM213-H", alpha=0.05)

data.remove_columns(["Event_2", "Sample label_1", "Sample label_2", "Altitude [m]", "Latitude_2", "Longitude_2"])
write_data(data, label="MSM213-H")
