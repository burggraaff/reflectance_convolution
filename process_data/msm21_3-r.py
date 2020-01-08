import numpy as np
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data

Ed = read("data/MSM21_3/MSM21_3_Ed-5nm.tab", data_start=142, header_start=141)
Lu = read("data/MSM21_3/MSM21_3_Lsfc-5nm.tab", data_start=142, header_start=141)
Ls = read("data/MSM21_3/MSM21_3_Lsky-5nm.tab", data_start=141, header_start=140)

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

R_rs_keys = [key for key in combined_table.keys() if "R_rs" in key]
remove_indices = [i for i, row in enumerate(combined_table) if any(row[key] <= -0.001 for key in R_rs_keys)]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with negative values")

remove_indices = [i for i, row in enumerate(combined_table) if row["R_rs_400"] < 0 or row["R_rs_800"] >= 0.003]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with values of R_rs(400 nm) < 0 or R_rs(800 nm) >= 0.003")

map_data(combined_table, data_label="MSM21_3R", projection='gnom', lat_0=66, lon_0=-40.5, llcrnrlon=-53, urcrnrlon=-12, llcrnrlat=58, urcrnrlat=70.5, resolution="h", parallels=np.arange(55, 75, 5), meridians=np.arange(-60, -5, 5))

plot_spectra(combined_table, data_label="MSM21_3R", alpha=0.05)

combined_table.remove_columns(["Latitude_1", "Longitude_1", "Altitude [m]_1", "Latitude_2", "Longitude_2", "Altitude [m]_2"])
write_data(combined_table, label="MSM21_3-R")
