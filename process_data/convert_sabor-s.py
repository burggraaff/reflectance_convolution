import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra
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

    data_table.remove_columns([key for key in data_table.keys() if "stokes" in key])
    data_table.remove_columns([key for key in data_table.keys() if "sd" in key])
    data_table.remove_columns([key for key in data_table.keys() if "sky" in key])
    data_table.remove_columns([key for key in data_table.keys() if "Lt" in key])
    # These data are normalised to R_rs(750 nm), so we must re-calculate Lw to get a fair comparison
    data_table.remove_columns([key for key in data_table.keys() if "Lw" in key])
    data_table.remove_columns([key for key in data_table.keys() if "AOT" in key])

    data_tables.append(data_table)
    print(file)

data = table.vstack(data_tables)

data.remove_column("col15346")
data.rename_column("lat", "Latitude")
data.rename_column("lon", "Longitude")

Es_keys = [key for key in data.keys() if "Es" in key]
R_rs_keys = [key for key in data.keys() if "Rrs" in key]

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

# Plot map of observations
fig = plt.figure(figsize=(10, 6), tight_layout=True)

m = Basemap(projection='gnom', lat_0=37, lon_0=-70, llcrnrlon=-77, urcrnrlon=-64, llcrnrlat=35, urcrnrlat=43, resolution="h")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(32, 45, 2), labels=[1,1,0,0])
m.drawmeridians(np.arange(-80, -60, 2), labels=[0,0,1,1])

m.scatter(data["Longitude"], data["Latitude"], latlon=True, c="r", edgecolors="k", s=60, zorder=10)

plt.savefig("data/plots/map_SABOR-S.pdf")
plt.show()

plot_spectra(data, data_label="SABOR-S", alpha=0.5)

data.write("data/sabor-s_processed.tab", format="ascii.fast_tab", overwrite=True)
