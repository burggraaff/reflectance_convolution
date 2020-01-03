import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u

Lw = read("data/MSM21_3/MSM21_3_Lw-5nm.tab", data_start=132, header_start=131)
Rrs = read("data/MSM21_3/MSM21_3_Rrs-5nm.tab", data_start=133, header_start=132)

wavelengths = np.arange(360, 805, 5)
for wvl in wavelengths:
    try:  # mu gets properly loaded on Linux
        Lw.rename_column(f"Lw_{wvl} [µW/cm**2/nm/sr]", f"Lw_{wvl}")
    except KeyError: # but not on Windows
        Lw.rename_column(f"Lw_{wvl} [ÂµW/cm**2/nm/sr]", f"Lw_{wvl}")
    Lw[f"Lw_{wvl}"].unit = u.microwatt / (u.centimeter**2 * u.nanometer * u.steradian)
    Lw[f"Lw_{wvl}"] = Lw[f"Lw_{wvl}"].to(u.watt / (u.meter**2 * u.nanometer * u.steradian))

    Rrs.rename_column(f"Rrs_{wvl} [1/sr]", f"R_rs_{wvl}")
    Rrs[f"R_rs_{wvl}"].unit = 1 / u.steradian

combined_table = table.join(Lw, Rrs, keys=["Date/Time"])
combined_table.rename_column("Event_1", "Event")
combined_table.rename_column("Latitude_1", "Latitude")
combined_table.rename_column("Longitude_1", "Longitude")

for wvl in wavelengths:
    Ed = combined_table[f"Lw_{wvl}"] / combined_table[f"R_rs_{wvl}"]
    Ed.name = f"Ed_{wvl}"
    Ed.unit = u.watt / (u.meter**2 * u.nanometer)
    combined_table.add_column(Ed)

R_rs_keys = [key for key in combined_table.keys() if "R_rs" in key]
remove_indices = [i for i, row in enumerate(combined_table) if any(row[key] <= -0.001 for key in R_rs_keys)]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with negative values")

remove_indices = [i for i, row in enumerate(combined_table) if row["R_rs_400"] < 0 or row["R_rs_800"] >= 0.003]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with values of R_rs(400 nm) < 0 or R_rs(800 nm) >= 0.003")

remove_indices = [i for i, row_mask in enumerate(combined_table.mask) if any(row_mask[key] for key in R_rs_keys)]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with NaN values")
# Plot map of observations
fig = plt.figure(figsize=(10, 10), tight_layout=True)

m = Basemap(projection='gnom', lat_0=66, lon_0=-40.5, llcrnrlon=-53, urcrnrlon=-12, llcrnrlat=58, urcrnrlat=70.5, resolution="h")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(55, 75, 5), labels=[1,1,0,0])
m.drawmeridians(np.arange(-60, -5, 5), labels=[0,0,1,1])

m.scatter(combined_table["Longitude"], combined_table["Latitude"], latlon=True, c="r", edgecolors="k", s=60, zorder=10)

plt.savefig("map_MSM21_3H.pdf")
plt.show()


# Plot all Ed, Lw, R_rs spectra
fig, axs = plt.subplots(nrows=3, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0})

for row in combined_table:
    spec_Ed = [row[f"Ed_{wvl}"] for wvl in wavelengths]
    spec_Lw = [row[f"Lw_{wvl}"] for wvl in wavelengths]
    spec_R_rs = [row[f"R_rs_{wvl}"] for wvl in wavelengths]

    for ax, spec in zip(axs.ravel(), [spec_Ed, spec_Lw, spec_R_rs]):
        ax.plot(wavelengths, spec, c="k", alpha=0.05, zorder=1)

for ax, label in zip(axs.ravel(), ["$E_d$ [W m$^{-2}$ nm$^{-1}$]", "$L_w$ [W m$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$R_{rs}$ [sr$^{-1}$]"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[0].tick_params(bottom=False, labelbottom=False)
axs[1].tick_params(bottom=False, labelbottom=False)
axs[2].set_xlabel("Wavelength [nm]")
axs[2].set_xlabel("Wavelength [nm]")
axs[0].set_xlim(360, 800)

axs[0].set_title(f"MSM21_3H spectra ({len(combined_table)})")
plt.savefig("spectra_MSM21_3H.pdf")
plt.show()
plt.close()

combined_table.remove_columns(["Event_2", "Sample label_1", "Sample label_2", "Altitude [m]", "Latitude_2", "Longitude_2"])
combined_table.write("data/msm21_3h_processed.tab", format="ascii.fast_tab", overwrite=True)
