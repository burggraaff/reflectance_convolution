import numpy as np
from astropy import table
from astropy import units as u
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data
from sba.data_processing import remove_negative_R_rs, get_keys_with_label, convert_to_unit, rename_columns, add_Lw_from_Ed_Rrs

Ed = read("data/SeaSWIR/SeaSWIR_TRIOS_Ed.tab", data_start=238, header_start=237)

wavelengths = np.arange(350, 902.5, 2.5)

Rrs = read("data/SeaSWIR/SeaSWIR_TRIOS_Rw.tab", data_start=473, format="no_header", delimiter="\t")
colnames = ["Event", "Campaign", "Station", "Location", "Comment (TRIOS missing?)", "Comment (ASD missing?)", "ID", "Latitude", "Longitude", "Date/Time (water sample, UTC)", "Date/Time (TRIOS start, UTC)", "Date/Time (TRIOS end, UTC)", "Ratio (drho/rho, 750nm)", "Ratio (dLsky/Lsky, 750nm)", "Ratio (dLw/Lw, 750nm)", "Ratio (dEd/Ed, 750nm)", "Ratio (d(Lsk/Ed)/(Lsk/Ed), 750nm)"]
colnames_rrs = [f"R_rs_{wvl:.1f}" for wvl in wavelengths]
colnames_rrs2 = [f"R_rs_err_{wvl:.1f}" for wvl in wavelengths]
colnames = colnames + colnames_rrs + colnames_rrs2

for j, newname in enumerate(colnames, 1):
    Rrs.rename_column(f"col{j}", newname)
Rrs.remove_columns(colnames_rrs2)
Rrs.remove_columns(['Ratio (drho/rho, 750nm)', 'Ratio (dLsky/Lsky, 750nm)', 'Ratio (dLw/Lw, 750nm)', 'Ratio (dEd/Ed, 750nm)', 'Ratio (d(Lsk/Ed)/(Lsk/Ed), 750nm)'])

data = table.join(Ed, Rrs, keys=["ID"])

rename_columns(data, "Ed [mW/m**2/nm] (", "Ed_", strip=True)

convert_to_unit(data, "Ed", u.milliwatt / (u.meter**2 * u.nanometer), u.watt / (u.meter**2 * u.nanometer))
convert_to_unit(data, "R_rs", 1 / u.steradian)

R_rs_keys = get_keys_with_label(data, "R_rs")
for R_rs_k in zip(R_rs_keys):
    # Convert R_w to R_rs
    data[R_rs_k] = data[R_rs_k] / np.pi

data = add_Lw_from_Ed_Rrs(data)

remove_negative_R_rs(data)

map_data(data, data_label="SeaSWIR-R", projection='merc', lat_0=10, lon_0=-30, llcrnrlon=-60, urcrnrlon=7, llcrnrlat=-38, urcrnrlat=55, resolution="i", figsize=(5,10), parallels=np.arange(-40, 60, 10), meridians=np.arange(-60, 20, 10))

plot_spectra(data, data_label="SeaSWIR-R", alpha=0.05)

data.remove_columns(["Event_1", "Event_2", "Campaign_1", "Campaign_2", "Station_1", "Station_2"])
write_data(data, label="SeaSWIR-R")
