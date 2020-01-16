import numpy as np
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import split_spectrum, get_keys_with_label, remove_negative_R_rs, convert_to_unit, rename_columns

Ed = read("data/SOP4/SO-P4_irrad.tab", data_start=142, header_start=141)
Lu = read("data/SOP4/SO-P4_rad_up_40deg.tab", data_start=142, header_start=141)
Ls = read("data/SOP4/SO-P4_sky_rad_40deg.tab", data_start=142, header_start=141)

data = table.join(Ed, Lu, keys=["Date/Time"])
data = table.join(data, Ls, keys=["Date/Time"])

rename_columns(data, "Ed", "Ed", strip=True)
rename_columns(data, "Lu", "Lu", strip=True)
rename_columns(data, "Ls", "Ls", strip=True)

Ed_keys, Lu_keys, Ls_keys = get_keys_with_label(data, "Ed", "Lu", "Ls")
for Ed_k, Lu_k, Ls_k in zip(Ed_keys, Lu_keys, Ls_keys):
    convert_to_unit(data, Ed_k, u.watt / (u.meter**2 * u.nanometer))
    convert_to_unit(data, Lu_k, u.microwatt / (u.centimeter**2 * u.nanometer * u.steradian), u.watt / (u.meter**2 * u.nanometer * u.steradian))
    convert_to_unit(data, Ls_k, u.watt / (u.meter**2 * u.nanometer * u.steradian))

    Lw_k = Lu_k.replace("Lu", "Lw")
    Lw = data[Lu_k] - 0.028 * data[Ls_k]
    Lw.name = Lw_k
    data.add_column(Lw)

    R_rs = data[Lw_k] / data[Ed_k]
    R_rs.name = Lw_k.replace("Lw", "R_rs")
    R_rs.unit = 1 / u.steradian
    data.add_column(R_rs)

# Normalise by R_rs(750 nm), re-calculate Lw
normalisation = data["R_rs_750"].copy()
Ed_keys, Lw_keys, R_rs_keys = get_keys_with_label(data, "Ed", "Lw", "R_rs")
for Ed_k, Lw_k, R_rs_k in zip(Ed_keys, Lw_keys, R_rs_keys):
    data[R_rs_k] = data[R_rs_k] - normalisation
    data[Lw_k] = data[R_rs_k] * data[Ed_k]

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
