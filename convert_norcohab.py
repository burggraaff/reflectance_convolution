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

    try:  # mu gets properly loaded on Linux
        Lu.rename_column(f"Lu_{wvl} [µW/cm**2/nm/sr]", f"Lu_{wvl}")
    except KeyError: # but not on Windows
        Lu.rename_column(f"Lu_{wvl} [ÂµW/cm**2/nm/sr]", f"Lu_{wvl}")
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

Lw_keys_NIR = [key for key in combined_table.keys() if "Lw_8" in key]
NIR_mean = np.mean([combined_table[key] for key in Lw_keys_NIR], axis=0)
combined_table.add_column(table.Column(data=NIR_mean, name="offset"))
for key in combined_table.keys():
    if "Lw" in key:
        combined_table[key] = combined_table[key] - combined_table["offset"]

remove_indices = [i for i, row in enumerate(combined_table) if row["Lw_400"] < 0]
combined_table.remove_rows(remove_indices)

for wvl in wavelengths:
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

plt.savefig("NORCOHAB_map.pdf")
plt.show()

# Plot all Ed, Lu, Lsky, Lw, R_rs spectra
fig, axs = plt.subplots(nrows=5, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,10))

for row in combined_table:
    spec_Ed = [row[f"Ed_{wvl}"] for wvl in wavelengths]
    spec_Lu = [row[f"Lu_{wvl}"] for wvl in wavelengths]
    spec_Ls = [row[f"Ls_{wvl}"] for wvl in wavelengths]
    spec_Lw = [row[f"Lw_{wvl}"] for wvl in wavelengths]
    spec_R_rs = [row[f"R_rs_{wvl}"] for wvl in wavelengths]

    for ax, spec in zip(axs.ravel(), [spec_Lu, spec_Ls, spec_Ed, spec_Lw, spec_R_rs]):
        ax.plot(wavelengths, spec, c="k", alpha=0.15, zorder=1)

for ax, label in zip(axs.ravel(), ["$L_u$ [$\mu$W cm$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$L_{sky}$ [$\mu$W cm$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$E_d$ [$\mu$W cm$^{-2}$ nm$^{-1}$]", "$L_w$ [$\mu$W cm$^{-2}$ nm$^{-1}$ sr$^{-1}$]", "$R_{rs}$ [sr$^{-1}$]"]):
    ax.set_ylabel(label)
    ax.grid(ls="--", zorder=0)

axs[-1].set_xlabel("Wavelength [nm]")
axs[-1].set_xlim(320, 950)

axs[0].set_title(f"NORCOHAB spectra ({len(combined_table)})")
plt.savefig("NORCOHAB_spectra.pdf")
plt.show()
plt.close()

combined_table.remove_columns(["Date/Time_1", "Latitude_1", "Longitude_1", "Altitude [m]_1", "Date/Time_2", "Latitude_2", "Longitude_2", "Altitude [m]_2"])
combined_table.write("data/norcohab_processed.tab", format="ascii.fast_tab", overwrite=True)

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
