import numpy as np
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import get_keys_with_label, remove_negative_R_rs, convert_to_unit, rename_columns, add_Lw_from_Ed_Rrs

wavelengths = np.arange(350, 1301, 1)

Ld = read("data/SeaSWIR/SeaSWIR_ASD_Ldspec.tab", data_start=974, header_start=973)
Ldkeys = get_keys_with_label(Ld, "Ld")
for Ldkey, wvl in zip(Ldkeys, wavelengths):
    # multiply by pi to convert L to E (Mobley99)
    # divide by 1e5 for normalisation to W/m^2/nm (empirical)
    Ed = Ld[Ldkey] * np.pi / 1e5
    Ed.name = f"Ed_{wvl}"
    Ld.add_column(Ed)
    Ld.remove_column(Ldkey)
Ed = Ld
# Units of Ld are not provided

Rrs = read("data/SeaSWIR/SeaSWIR_ASD_Rw.tab", data_start=974, header_start=973)
rename_columns(Rrs, "Refl (", "R_rs_", strip=True)
R_rs_keys = get_keys_with_label(Rrs, "R_rs")
for R_rs_k in R_rs_keys:
    # Convert R_w to R_rs
    Rrs[R_rs_k] = Rrs[R_rs_k] / np.pi

data = table.join(Ed, Rrs, keys=["Station"])

convert_to_unit(data, "Ed", u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "R_rs", 1 / u.steradian)

data = add_Lw_from_Ed_Rrs(data)

for key in data.keys():
    if key[-2:] == "_1":
        data.rename_column(key, key[:-2])
    elif key[-2:] == "_2":
        data.remove_column(key)

# Check for rows where no Ed data were provided
remove_indices = [i for i, mask_row in enumerate(data.mask) if mask_row["Ed_500"]]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows without Ed data")

remove_negative_R_rs(data)

map_data(data, data_label="SeaSWIR-A", projection='merc', lat_0=10, lon_0=-30, llcrnrlon=-60, urcrnrlon=7, llcrnrlat=-38, urcrnrlat=55, resolution="i", figsize=(5,10), parallels=np.arange(-40, 60, 10), meridians=np.arange(-60, 20, 10))

plot_spectra(data, data_label="SeaSWIR-A", alpha=0.05)

data.remove_columns(["Date/Time (end, UTC)", "Date/Time (end, local time)", "Date/Time (start, local time)"])
write_data(data, label="SeaSWIR-A")
