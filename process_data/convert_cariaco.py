import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra

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

    with open(file, "r") as f:
        lines = f.readlines()
        lat_line = [line for line in lines if "north_latitude" in line][0]
        lat = float(lat_line[16:-6])
        lon_line = [line for line in lines if "west_longitude" in line][0]
        lon = float(lon_line[16:-6])

        date_line = [line for line in lines if "end_date" in line][0]
        date = int(date_line[10:])
        time_line = [line for line in lines if "start_time" in line][0]
        time = time_line[12:-6]

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

# Plot map of observations
fig = plt.figure(figsize=(10, 6), tight_layout=True)

m = Basemap(projection='gnom', lat_0=10.5, lon_0=-64.67, llcrnrlon=-70, urcrnrlon=-59, llcrnrlat=5, urcrnrlat=15, resolution="h")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(4, 16, 2), labels=[1,1,0,0])
m.drawmeridians(np.arange(-70, -56, 2), labels=[0,0,1,1])

m.scatter(data["Longitude"], data["Latitude"], latlon=True, c="r", edgecolors="k", s=60, zorder=10)

plt.savefig("data/plots/map_CARIACO.pdf")
plt.show()

plot_spectra(data, data_label="CARIACO", alpha=0.1)

data.write("data/cariaco_processed.tab", format="ascii.fast_tab", overwrite=True)
