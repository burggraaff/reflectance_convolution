import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from pathlib import Path

folder = Path("data")

all_data_files = list(folder.glob("*processed.tab"))

data = [read(file, include_names=["Latitude", "Longitude"]) for file in all_data_files]
labels = [file.stem[:-10] for file in all_data_files]

for label, tab in zip(labels, data):
    tab.add_column(table.Column(name="Label", data=[label for row in tab]))

data = table.vstack(data)

# Plot map of observations
fig = plt.figure(figsize=(8, 5), tight_layout=True)

m = Basemap(projection='moll', lon_0=0, resolution="i")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(-90, 95, 15), labels=[1,1,0,0])
m.drawmeridians(np.arange(-180, 180, 30), labels=[0,0,1,1])

m.scatter(data["Longitude"], data["Latitude"], latlon=True, c="r", edgecolors="", s=30, zorder=10)

plt.savefig("data/plots/map_all_data.pdf")
plt.show()
