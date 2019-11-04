import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u

Ed = read("data/NORCOHAB/HE302_irrad.tab", data_start=186, header_start=185)
Lu = read("data/NORCOHAB/HE302_rad.tab", data_start=186, header_start=185)
Ls = read("data/NORCOHAB/HE302_ssr.tab", data_start=186, header_start=185)

wavelengths = np.arange(320, 955, 5)
for wvl in wavelengths:
    Ed.rename_column(f"Ed_{wvl} [W/m**2/nm]", f"Ed_{wvl}")
    Ed[f"Ed_{wvl}"].unit = u.watt / (u.meter**2 * u.nanometer)

    Lu.rename_column(f"Lu_{wvl} [µW/cm**2/nm/sr]", f"Lu_{wvl}")
    Lu[f"Lu_{wvl}"].unit = u.microwatt / (u.centimeter**2 * u.nanometer * u.steradian)
    Lu[f"Lu_{wvl}"] = Lu[f"Lu_{wvl}"].to(u.watt / (u.meter**2 * u.nanometer * u.steradian))

    Ls.rename_column(f"Ls_{wvl} [W/m**2/nm/sr]", f"Ls_{wvl}")
    Ls[f"Ls_{wvl}"].unit = u.watt / (u.meter**2 * u.nanometer * u.steradian)

combined_table = table.join(Ed, Lu, keys=["Event"])
combined_table = table.join(combined_table, Ls, keys=["Event"])

for wvl in wavelengths:
    Lw = combined_table[f"Lu_{wvl}"] - 0.028 * combined_table[f"Ls_{wvl}"]
    Lw.name = f"Lw_{wvl}"
    combined_table.add_column(Lw)
    R_rs = combined_table[f"Lw_{wvl}"] / combined_table[f"Ed_{wvl}"]
    R_rs.name = f"R_rs_{wvl}"
    R_rs.unit = 1 / u.steradian
    combined_table.add_column(R_rs)

# Plot map of observations
fig = plt.figure(figsize=(10, 7.5), tight_layout=True)

m = Basemap(projection='gnom', lat_0=55, lon_0=0, llcrnrlon=-10, urcrnrlon=11, llcrnrlat=50.5, urcrnrlat=59.5, resolution="h")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(40, 70, 2), labels=[1,1,0,0])
m.drawmeridians(np.arange(-20, 20, 2), labels=[0,0,1,1])

m.scatter(combined_table["Longitude"], combined_table["Latitude"], latlon=True, c="r", edgecolors="k", s=60)

plt.savefig("NORCOHAB.pdf")
plt.show()

