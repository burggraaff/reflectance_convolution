import numpy as np
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import get_keys_with_label, convert_to_unit, rename_columns

Ed = read("data/HE302/HE302_irrad.tab", data_start=186, header_start=185)
Rrs = read("data/HE302/HE302_rrs.tab", data_start=186, header_start=185)

data = table.join(Ed, Rrs, keys=["Event"])

rename_columns(data, "Ed", "Ed", strip=True)
rename_columns(data, "Rrs", "R_rs", strip=True)

Ed_keys, R_rs_keys = get_keys_with_label(data, "Ed", "R_rs")
for Ed_k, R_rs_k in zip(Ed_keys, R_rs_keys):
    convert_to_unit(data, Ed_k, u.watt / (u.meter**2 * u.nanometer))
    convert_to_unit(data, R_rs_k, 1 / u.steradian)

    Lw = data[Ed_k] * data[R_rs_k]
    Lw.name = Ed_k.replace("Ed", "Lw")
    Lw.unit = u.watt / (u.meter**2 * u.nanometer * u.steradian)
    data.add_column(Lw)

remove_indices = [i for i, row in enumerate(data) if row["R_rs_800"] >= 0.003]
data.remove_rows(remove_indices)
print(f"Removed {len(remove_indices)} rows with values of R_rs(800 nm) >= 0.003")

for key in ["Date/Time", "Latitude", "Longitude", "Altitude [m]"]:
    data.rename_column(f"{key}_1", key)
    data.remove_column(f"{key}_2")

map_data(data, data_label="HE302", projection='gnom', lat_0=55, lon_0=0, llcrnrlon=-10, urcrnrlon=11, llcrnrlat=50.5, urcrnrlat=59.5, resolution="h", parallels=np.arange(40, 70, 2), meridians=np.arange(-20, 20, 2))

plot_spectra(data, data_label="HE302", alpha=0.15)

write_data(data, label="HE302")
