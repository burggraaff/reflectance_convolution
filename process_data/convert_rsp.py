import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra

folder = Path("data/RSP/")
files = list(folder.glob("*.txt"))

data_tables = []

for file in files:
    data = read(file, data_start=37)
    header = read(file, data_start=32, data_end=33)
    header["col1"][0] = "date"
    header = header[0].as_void()
    for key, new_key in zip(data.keys(), header):
        data.rename_column(key, new_key)

    data.rename_column("lat", "Latitude")
    data.rename_column("lon", "Longitude")

    data.remove_columns([key for key in data.keys() if "sd" in key])
    data.remove_columns([key for key in data.keys() if "Lu" in key])

    data_tables.append(data)

data = table.vstack(data_tables)

Ed_keys = [key for key in data.keys() if "Ed" in key]
R_rs_keys = [key for key in data.keys() if "Rrs" in key]

for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    wavelength = float(Ed_k[2:])

    data[Ed_k].unit = u.microwatt / (u.cm**2 * u.nm)
    data[Ed_k] = data[Ed_k].to(u.watt / (u.m**2 * u.nm))
    data.rename_column(Ed_k, f"Ed_{wavelength:.1f}")

    data[R_rs_k].unit = 1 / u.steradian
    data.rename_column(R_rs_k, f"R_rs_{wavelength:.1f}")

    Lw = data[f"Ed_{wavelength:.1f}"] * data[f"R_rs_{wavelength:.1f}"]
    Lw.name = f"Lw_{wavelength:.1f}"
    Lw.unit = u.watt / (u.m**2 * u.nm * u.steradian)
    data.add_column(Lw)

# Plot map of observations
fig = plt.figure(figsize=(10, 6), tight_layout=True)

m = Basemap(projection='moll', lon_0=0, resolution="i")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(-90, 95, 15), labels=[1,1,0,0])
m.drawmeridians(np.arange(-180, 180, 30), labels=[0,0,1,1])

m.scatter(data["Longitude"], data["Latitude"], latlon=True, c="r", edgecolors="k", s=60, zorder=10)

plt.savefig("data/plots/map_RSP.pdf")
plt.show()

plot_spectra(data, data_label="RSP", alpha=0.3)

data.write("data/rsp_processed.tab", format="ascii.fast_tab", overwrite=True)