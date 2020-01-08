import numpy as np
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data

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

# Remove rows with jumps in Ed >0.5 between wavelengths
Ed_keys = [key for key in combined_table.keys() if "Ed" in key]
remove_indices = [i for i, row in enumerate(combined_table) if any([np.abs(row[key1]-row[key2]) >= 0.5 for key1, key2 in zip(Ed_keys, Ed_keys[1:])])]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with Ed jumps > 0.5")

# Remove rows with negative R_rs
R_rs_keys = [key for key in combined_table.keys() if "R_rs" in key]
remove_indices = [i for i, row in enumerate(combined_table) if any(row[key] <= -0.001 for key in R_rs_keys)]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with negative values")

# Remove rows with NaN values
remove_indices = [i for i, row_mask in enumerate(combined_table.mask) if any(row_mask[key] for key in R_rs_keys)]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with NaN values")

map_data(combined_table, data_label="MSM21_3H", projection='gnom', lat_0=66, lon_0=-40.5, llcrnrlon=-53, urcrnrlon=-12, llcrnrlat=58, urcrnrlat=70.5, resolution="h", parallels=np.arange(55, 75, 5), meridians=np.arange(-60, -5, 5))

plot_spectra(combined_table, data_label="MSM21_3H", alpha=0.05)

combined_table.remove_columns(["Event_2", "Sample label_1", "Sample label_2", "Altitude [m]", "Latitude_2", "Longitude_2"])
combined_table.write("data/msm21_3h_processed.tab", format="ascii.fast_tab", overwrite=True)
