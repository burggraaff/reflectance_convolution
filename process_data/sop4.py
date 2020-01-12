import numpy as np
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import split_spectrum, get_keys_with_label, remove_negative_R_rs

Ed = read("data/SOP4/SO-P4_irrad.tab", data_start=142, header_start=141)
Lu = read("data/SOP4/SO-P4_rad_up_40deg.tab", data_start=142, header_start=141)
Ls = read("data/SOP4/SO-P4_sky_rad_40deg.tab", data_start=142, header_start=141)

wavelengths = np.arange(320, 955, 5)
for wvl in wavelengths:
    Ed.rename_column(f"Ed_{wvl} [W/m**2/nm]", f"Ed_{wvl}")
    Ed[f"Ed_{wvl}"].unit = u.watt / (u.meter**2 * u.nanometer)

    try:  # mu gets properly loaded on Linux
        Lu.rename_column(f"Lu_{wvl} [µW/cm**2/nm/sr]", f"Lu_{wvl}")
    except KeyError: # but not on Windows
        Lu.rename_column(f"Lu_{wvl} [ÂµW/cm**2/nm/sr]", f"Lu_{wvl}")
    Lu[f"Lu_{wvl}"].unit = u.microwatt / (u.centimeter**2 * u.nanometer * u.steradian)
    Lu[f"Lu_{wvl}"] = Lu[f"Lu_{wvl}"].to(u.watt / (u.meter**2 * u.nanometer * u.steradian))

    Ls.rename_column(f"Ls_{wvl} [W/m**2/nm/sr]", f"Ls_{wvl}")
    Ls[f"Ls_{wvl}"].unit = u.watt / (u.meter**2 * u.nanometer * u.steradian)

data = table.join(Ed, Lu, keys=["Date/Time"])
data = table.join(data, Ls, keys=["Date/Time"])

for wvl in wavelengths:
    Lw = data[f"Lu_{wvl}"] - 0.028 * data[f"Ls_{wvl}"]
    Lw.name = f"Lw_{wvl}"
    data.add_column(Lw)
    R_rs = data[f"Lw_{wvl}"] / data[f"Ed_{wvl}"]
    R_rs.name = f"R_rs_{wvl}"
    R_rs.unit = 1 / u.steradian
    data.add_column(R_rs)

# Normalise by R_rs(750 nm), re-calculate Lw
normalisation = data["R_rs_750"].copy()
for wvl in wavelengths:
    data[f"R_rs_{wvl}"] = data[f"R_rs_{wvl}"] - normalisation
    data[f"Lw_{wvl}"] = data[f"R_rs_{wvl}"] * data[f"Ed_{wvl}"]

# Remove rows where the maximum Ed is unphysically small (< threshold)
threshold = 0.01
wavelengths, Ed = split_spectrum(data, "Ed")
max_values = np.max(Ed, axis=1)
remove_indices = np.where(max_values < threshold)[0]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with max(Ed) < {threshold}")

# Remove columns that are consistently negative in R_rs
for wvl in wavelengths[(wavelengths < 360) | (wavelengths > 750)]:
    remove_keys = get_keys_with_label(data, f"{wvl:.0f}")
    data.remove_columns(remove_keys)
print("Removed columns with wavelengths <360 and >750")

# Remove rows with consecutive jumps in Ed >= threshold between wavelengths
threshold = 0.005
wavelengths, R_rs = split_spectrum(data, "R_rs")
diffs = np.diff(R_rs, axis=1)
diffs_abs = np.abs(diffs)
all_rows = [np.where((col1 >= threshold) & (col2 >= threshold))[0] for col1, col2 in zip(diffs_abs.T, diffs_abs.T[1:])]
all_rows = [row for sub in all_rows for row in sub]
remove_indices = np.unique(all_rows)
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with consecutive R_rs jumps >= {threshold}")

# Remove rows where Ed_405 is abnormally low compared to Ed_400
diff = data["Ed_400"] - data["Ed_405"]
remove_indices = np.where(diff > 0.01)[0]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows where Ed(400 nm) - Ed(405 nm) > 0.01")

remove_negative_R_rs(data)

map_data(data, data_label="SOP4", projection='gnom', lat_0=56, lon_0=5, llcrnrlon=-2, urcrnrlon=11, llcrnrlat=52, urcrnrlat=59, resolution="h", parallels=np.arange(40, 70, 2), meridians=np.arange(-20, 20, 2))

plot_spectra(data, data_label="SOP4", alpha=0.05)

data.remove_columns(["Latitude_1", "Longitude_1", "Altitude [m]_1", "Latitude_2", "Longitude_2", "Altitude [m]_2"])
write_data(data, "SOP4")
