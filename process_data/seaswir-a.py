import numpy as np
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data

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

data = table.join(Ed, Rrs, keys=["Station"])

for wvl in wavelengths:
    Lw = data[f"R_rs_{wvl}"] * data[f"Ed_{wvl}"]
    Lw.name = f"Lw_{wvl}"
    Lw.unit = u.watt / (u.meter**2 * u.nanometer * u.steradian)
    data.add_column(Lw)

for key in data.keys():
    if key[-2:] == "_1":
        data.rename_column(key, key[:-2])
    elif key[-2:] == "_2":
        data.remove_column(key)

# Check for rows where no Ed data were provided
remove_indices = [i for i, mask_row in enumerate(data.mask) if mask_row["Ed_500"]]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows without Ed data")

R_rs_keys = [key for key in data.keys() if "R_rs" in key]
remove_indices = [i for i, row in enumerate(data) if any(row[key] <= -0.001 for key in R_rs_keys)]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with negative values")

#remove_indices = [i for i, row in enumerate(data) if any(row[key] > 0.1 for key in R_rs_keys)]
#data.remove_rows(remove_indices)
#print(f"Removed {len(remove_indices)} rows with values > 0.1")

remove_indices = [i for i, row in enumerate(data) if row["R_rs_400"] < 0]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with values of R_rs(400 nm) < 0 or R_rs(800 nm) >= 0.003")

map_data(data, data_label="SeaSWIR-A", projection='merc', lat_0=10, lon_0=-30, llcrnrlon=-60, urcrnrlon=7, llcrnrlat=-38, urcrnrlat=55, resolution="i", figsize=(5,10), parallels=np.arange(-40, 60, 10), meridians=np.arange(-60, 20, 10))

plot_spectra(data, data_label="SeaSWIR-A", alpha=0.05)

data.remove_columns(["Date/Time (end, UTC)", "Date/Time (end, local time)", "Date/Time (start, local time)"])
write_data(data, label="SeaSWIR-A")
