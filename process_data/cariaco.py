import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass

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

    # reject data with negative Rrs
    if len(np.where(Rrs < 0)[0]) > 0:
        continue

    # reject data that do not cover the full spectral range
    if 380 not in wavelengths:
        continue

    date, time, lon, lat = find_auxiliary_information_seabass(file)

    cols = ["Date", "Time", "Latitude", "Longitude"] + [f"Ed_{wvl:.0f}" for wvl in wavelengths] + [f"R_rs_{wvl:.0f}" for wvl in wavelengths]
    dtype = [int, "S8", float, float] + 2 * [float for wvl in wavelengths]
    tab = table.Table(rows=[[date, time, lat, lon, *Ed, *Rrs]], names=cols, dtype=dtype)
    tabs.append(tab)

data = table.vstack(tabs)

Ed_keys = [key for key in data.keys() if "Ed" in key]
R_rs_keys = [key for key in data.keys() if "R_rs" in key]

normalisation = data["R_rs_748"].copy()
for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    wavelength = int(Ed_k[3:])

    data[Ed_k].unit = u.microwatt / (u.cm**2 * u.nm)
    data[Ed_k] = data[Ed_k].to(u.watt / (u.m**2 * u.nm))

    data[R_rs_k].unit = 1 / u.steradian
    # Normalise to R_rs at 748 nm
    data[R_rs_k] = data[R_rs_k] - normalisation

    Lw = data[f"Ed_{wavelength}"] * data[f"R_rs_{wavelength}"]
    Lw.name = f"Lw_{wavelength}"
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

# Remove columns that are consistently negative in R_rs
for wvl in wavelengths[wavelengths >= 680]:
    remove_keys = [f"{s}_{wvl:.0f}" for s in ("Ed", "Lw", "R_rs")]
    data.remove_columns(remove_keys)

R_rs_keys = [key for key in data.keys() if "R_rs" in key]
remove_indices = [i for i, row in enumerate(data) if any(row[key] > 0.2 for key in R_rs_keys)]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with R_rs > 0.2")

remove_indices = [i for i, row in enumerate(data) if row["Ed_600"] <= 0.1]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with missing values")

map_data(data, data_label="CARIACO", projection='gnom', lat_0=10.5, lon_0=-64.67, llcrnrlon=-70, urcrnrlon=-59, llcrnrlat=5, urcrnrlat=15, resolution="h", parallels=np.arange(4, 16, 2), meridians=np.arange(-70, -56, 2))

plot_spectra(data, data_label="CARIACO", alpha=0.1)

write_data(data, label="CARIACO")
