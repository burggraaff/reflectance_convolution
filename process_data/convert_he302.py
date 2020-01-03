import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
from astropy.io.ascii import read
from astropy import table
from astropy import units as u

Ed = read("data/HE302/HE302_irrad.tab", data_start=186, header_start=185)
Lu = read("data/HE302/HE302_rad.tab", data_start=186, header_start=185)
Ls = read("data/HE302/HE302_ssr.tab", data_start=186, header_start=185)
Rrs = read("data/HE302/HE302_rrs.tab", data_start=186, header_start=185)

wavelengths = np.arange(320, 955, 5)
for wvl in wavelengths:
    Ed.rename_column(f"Ed_{wvl} [W/m**2/nm]", f"Ed_{wvl}")
    Ed[f"Ed_{wvl}"].unit = u.watt / (u.meter**2 * u.nanometer)

    Rrs.rename_column(f"Rrs_{wvl} [1/sr]", f"R_rs_{wvl}")
    Rrs[f"R_rs_{wvl}"].unit = 1 / u.steradian

combined_table = table.join(Ed, Rrs, keys=["Event"])

for wvl in wavelengths:
    Lw = combined_table[f"Ed_{wvl}"] * combined_table[f"R_rs_{wvl}"]
    Lw.name = f"Lw_{wvl}"
    Lw.unit = u.watt / (u.meter**2 * u.nanometer)
    combined_table.add_column(Lw)

remove_indices = [i for i, row in enumerate(combined_table) if row["R_rs_400"] < 0 or row["R_rs_800"] >= 0.003]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with values of R_rs(400 nm) < 0 or R_rs(800 nm) >= 0.003")

for key in ["Date/Time", "Latitude", "Longitude", "Altitude [m]"]:
    combined_table.rename_column(f"{key}_1", key)
    combined_table.remove_column(f"{key}_2")

# Plot map of observations
fig = plt.figure(figsize=(10, 7.5), tight_layout=True)

m = Basemap(projection='gnom', lat_0=55, lon_0=0, llcrnrlon=-10, urcrnrlon=11, llcrnrlat=50.5, urcrnrlat=59.5, resolution="h")
m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
m.drawmapboundary(fill_color="#DDEEFF")
m.drawcoastlines()

m.drawparallels(np.arange(40, 70, 2), labels=[1,1,0,0])
m.drawmeridians(np.arange(-20, 20, 2), labels=[0,0,1,1])

m.scatter(combined_table["Longitude"], combined_table["Latitude"], latlon=True, c="r", edgecolors="k", s=60, zorder=10)

plt.savefig("data/plots/map_HE302.pdf")
plt.show()

# Plot all Ed, Lu, Lsky, Lw, R_rs spectra
fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,10))

for row in combined_table:
    spec_Lw = [row[f"Lw_{wvl}"] for wvl in wavelengths]
    spec_Ed = [row[f"Ed_{wvl}"] for wvl in wavelengths]
    spec_R_rs = [row[f"R_rs_{wvl}"] for wvl in wavelengths]

    for ax, spec in zip(axs.ravel(), [spec_Lw, spec_Ed, spec_R_rs]):
        ax.plot(wavelengths, spec, c="k", alpha=0.15, zorder=1)

for ax, label in zip(axs.ravel(), ["$L_w$ [W cm$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$E_d$ [W cm$^{-2}$ nm$^{-1}$]", "$R_{rs}$ [sr$^{-1}$]"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[-1].set_xlabel("Wavelength [nm]")
axs[-1].set_xlim(320, 950)

axs[0].set_title(f"HE302 spectra ({len(combined_table)})")
plt.savefig("data/plots/spectra_HE302.pdf")
plt.show()
plt.close()

combined_table.write("data/he302_processed.tab", format="ascii.fast_tab", overwrite=True)

wavelengths_interp = np.arange(380, 800.5, 0.5)

def interpolate_row(row, spectrum):
    spectrum_data = np.array([row[f"{spectrum}_{wvl}"] for wvl in wavelengths])
    spectrum_interpolated = np.interp(wavelengths_interp, wavelengths, spectrum_data, left=0, right=0)
    return spectrum_interpolated

def interpolate_table(data_table, spectrum):
    interpolated_data = np.array([interpolate_row(row, spectrum) for row in data_table])
    table_names = [f"{spectrum}_{wvl}" for wvl in wavelengths_interp]
    interpolated_table = table.Table(data=interpolated_data, names=table_names)
    return interpolated_table
