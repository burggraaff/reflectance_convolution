import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u

Ed = read("data/ARCHEMHAB/MSM21_3_Ed-5nm.tab", data_start=142, header_start=141)
Lu = read("data/ARCHEMHAB/MSM21_3_Lsfc-5nm.tab", data_start=142, header_start=141)
Ls = read("data/ARCHEMHAB/MSM21_3_Lsky-5nm.tab", data_start=141, header_start=140)

wavelengths = np.arange(320, 955, 5)
for wvl in wavelengths:
    Ed.rename_column(f"Ed_{wvl} [W/m**2/nm]", f"Ed_{wvl}")
    Ed[f"Ed_{wvl}"].unit = u.watt / (u.meter**2 * u.nanometer)

    Lu.rename_column(f"Lu_{wvl} [µW/cm**2/nm/sr]", f"Lu_{wvl}")
    Lu[f"Lu_{wvl}"].unit = u.microwatt / (u.centimeter**2 * u.nanometer * u.steradian)
    Lu[f"Lu_{wvl}"] = Lu[f"Lu_{wvl}"].to(u.watt / (u.meter**2 * u.nanometer * u.steradian))

    Ls.rename_column(f"Ls_{wvl} [W/m**2/nm/sr]", f"Ls_{wvl}")
    Ls[f"Ls_{wvl}"].unit = u.watt / (u.meter**2 * u.nanometer * u.steradian)

combined_table = table.join(Ed, Lu, keys=["Date/Time"])
combined_table = table.join(combined_table, Ls, keys=["Date/Time"])

for wvl in wavelengths:
    Lw = combined_table[f"Lu_{wvl}"] - 0.028 * combined_table[f"Ls_{wvl}"]
    Lw.name = f"Lw_{wvl}"
    combined_table.add_column(Lw)
    R_rs = combined_table[f"Lw_{wvl}"] / combined_table[f"Ed_{wvl}"]
    R_rs.name = f"R_rs_{wvl}"
    R_rs.unit = 1 / u.steradian
    combined_table.add_column(R_rs)

# Plot map of observations
fig = plt.figure(figsize=(10, 10), tight_layout=True)

m = Basemap(projection='gnom', lat_0=66, lon_0=-40.5, llcrnrlon=-53, urcrnrlon=-12, llcrnrlat=58, urcrnrlat=70.5, resolution="h")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(55, 75, 5), labels=[1,1,0,0])
m.drawmeridians(np.arange(-60, -5, 5), labels=[0,0,1,1])

m.scatter(combined_table["Longitude"], combined_table["Latitude"], latlon=True, c="r", edgecolors="k", s=60)

plt.savefig("ARCHEMHAB_map.pdf")
plt.show()

# Plot all Ed, Lu, Ls, R_rs spectra
fig, axs = plt.subplots(nrows=2, ncols=2, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0})

for row in combined_table:
    spec_Ed = [row[f"Ed_{wvl}"] for wvl in wavelengths]
    spec_Lu = [row[f"Lu_{wvl}"] for wvl in wavelengths]
    spec_Ls = [row[f"Ls_{wvl}"] for wvl in wavelengths]
    spec_R_rs = [row[f"R_rs_{wvl}"] for wvl in wavelengths]

    for ax, spec in zip(axs.ravel(), [spec_Ed, spec_Lu, spec_Ls, spec_R_rs]):
        ax.plot(wavelengths, spec, c="k", alpha=0.02, zorder=1)

for ax, label in zip(axs.ravel(), ["$E_d$ [W m$^{-2}$ nm$^{-1}$]", "$L_u$ [W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$L_s$ [W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$R_{rs}$ [sr$^{-1}$]"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[0,0].tick_params(bottom=False, labelbottom=False)
axs[0,1].tick_params(bottom=False, labelbottom=False)
axs[0,1].tick_params(left=False, labelleft=False, right=True, labelright=True)
axs[1,1].tick_params(left=False, labelleft=False, right=True, labelright=True)
axs[0,1].yaxis.set_label_position("right")
axs[1,1].yaxis.set_label_position("right")
axs[1,0].set_xlabel("Wavelength [nm]")
axs[1,1].set_xlabel("Wavelength [nm]")
axs[0,0].set_xlim(320, 950)

fig.suptitle("ARCHEMHAB spectra")
plt.savefig("ARCHEMHAB_spectra.pdf")
plt.show()
plt.close()
