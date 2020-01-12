import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import get_keys_with_label, remove_negative_R_rs, remove_rows_based_on_threshold

folder = Path("data/CARIACO/")
files = sorted(folder.glob("*.txt"))

tabs = []
for file in files:
    for skiprows in range(40, 70):
        try:
            wavelengths, Ed, Rrs = np.loadtxt(file, skiprows=skiprows, unpack=True, usecols=[0,3,5])
        except Exception:
            continue
        else:
            break

    date, time, lon, lat = find_auxiliary_information_seabass(file)

    cols = ["Date", "Time", "Latitude", "Longitude"] + [f"Ed_{wvl:.0f}" for wvl in wavelengths] + [f"R_rs_{wvl:.0f}" for wvl in wavelengths]
    dtype = [int, "S8", float, float] + 2 * [float for wvl in wavelengths]
    tab = table.Table(rows=[[date, time, lat, lon, *Ed, *Rrs]], names=cols, dtype=dtype)
    tabs.append(tab)

data = table.vstack(tabs)

Ed_keys, R_rs_keys = get_keys_with_label(data, "Ed", "R_rs")

for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    wavelength = int(Ed_k[3:])

    data[Ed_k].unit = u.microwatt / (u.cm**2 * u.nm)
    data[Ed_k] = data[Ed_k].to(u.watt / (u.m**2 * u.nm))

    data[R_rs_k].unit = 1 / u.steradian

    Lw = data[f"Ed_{wavelength}"] * data[f"R_rs_{wavelength}"]
    Lw.name = f"Lw_{wavelength}"
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

# Remove rows with NaN values
remove_indices = [i for i, row_mask in enumerate(data.mask) if any(row_mask[key] for key in R_rs_keys)]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with NaN values")

# Remove rows with missing R_rs values (< -90)
remove_rows_based_on_threshold(data, "R_rs", "<", -90)

# Remove rows with missing Ed values (<)
remove_indices = [i for i, row in enumerate(data) if row["Ed_600"] <= 0.1]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with missing Ed values")

# Remove rows with unphysically high R_rs values (> 0.02)
remove_indices = [i for i, row in enumerate(data) if any(row[key] > 0.02 for key in R_rs_keys)]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with R_rs > 0.02")

remove_negative_R_rs(data)

map_data(data, data_label="CARIACO", projection='gnom', lat_0=10.5, lon_0=-64.67, llcrnrlon=-70, urcrnrlon=-59, llcrnrlat=5, urcrnrlat=15, resolution="h", parallels=np.arange(4, 16, 2), meridians=np.arange(-70, -56, 2))

plot_spectra(data, data_label="CARIACO", alpha=0.1)

write_data(data, label="CARIACO")
