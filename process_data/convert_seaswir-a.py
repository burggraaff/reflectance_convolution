import numpy as np
from astropy.io.ascii import read
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data

wavelengths = np.arange(350, 1301, 1)

Ld = read("data/SeaSWIR/SeaSWIR_ASD_Ldspec.tab", data_start=974, header_start=973)
Ldkeys = [key for key in Ld.keys() if "Ld" in key]
for Ldkey, wvl in zip(Ldkeys, wavelengths):
    # multiply by pi to convert L to E (integrate over hemisphere)
    # divide by 1e5 for normalisation to W/m^2/nm (empirical)
    Ed = Ld[Ldkey] * np.pi / 1e5
    Ed.name = f"Ed_{wvl}"
    Ed.unit = u.watt / (u.m**2 * u.nm)
    Ld.add_column(Ed)
    Ld.remove_column(Ldkey)
Ed = Ld
# Units of Ld are not provided

Rrs = read("data/SeaSWIR/SeaSWIR_ASD_Rw.tab", data_start=974, header_start=973)
Rrskeys = [key for key in Rrs.keys() if "Refl" in key]
for key, wvl in zip(Rrskeys, wavelengths):
    Rrs.rename_column(key, f"R_rs_{wvl}")

combined_table = table.join(Ed, Rrs, keys=["Station"])

for wvl in wavelengths:
    Lw = combined_table[f"R_rs_{wvl}"] * combined_table[f"Ed_{wvl}"]
    Lw.name = f"Lw_{wvl}"
    Lw.unit = u.watt / (u.meter**2 * u.nanometer * u.steradian)
    combined_table.add_column(Lw)

for key in combined_table.keys():
    if key[-2:] == "_1":
        combined_table.rename_column(key, key[:-2])
    elif key[-2:] == "_2":
        combined_table.remove_column(key)

# Check for rows where no Ed data were provided
remove_indices = [i for i, mask_row in enumerate(combined_table.mask) if mask_row["Ed_500"]]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows without Ed data")

R_rs_keys = [key for key in combined_table.keys() if "R_rs" in key]
remove_indices = [i for i, row in enumerate(combined_table) if any(row[key] <= -0.001 for key in R_rs_keys)]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with negative values")

#remove_indices = [i for i, row in enumerate(combined_table) if any(row[key] > 0.1 for key in R_rs_keys)]
#combined_table.remove_rows(remove_indices)
#print(f"Removed {len(remove_indices)} rows with values > 0.1")

remove_indices = [i for i, row in enumerate(combined_table) if row["R_rs_400"] < 0]
combined_table.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with values of R_rs(400 nm) < 0 or R_rs(800 nm) >= 0.003")

map_data(combined_table, data_label="SeaSWIR-A", projection='merc', lat_0=10, lon_0=-30, llcrnrlon=-60, urcrnrlon=7, llcrnrlat=-38, urcrnrlat=55, resolution="i", figsize=(5,10), parallels=np.arange(-40, 60, 10), meridians=np.arange(-60, 20, 10))

plot_spectra(combined_table, data_label="SeaSWIR-A", alpha=0.05)

combined_table.remove_columns(["Date/Time (end, UTC)", "Date/Time (end, local time)", "Date/Time (start, local time)"])
combined_table.write("data/seaswir-a_processed.tab", format="ascii.fast_tab", overwrite=True)
