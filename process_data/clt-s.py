import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import split_spectrum, remove_rows_based_on_threshold, get_keys_with_label
from datetime import datetime

folder = Path("data/CLT/HyperSAS/")
master_files = list(folder.glob("CLT*.txt"))

master_table = table.vstack([read(file, data_start=48, header_start=46, include_names=["!Station", "time_GMT"]) for file in master_files])
master_table.rename_column("!Station", "Station")

# Remove data without a given 3-minute window
master_table.remove_rows(np.where(master_table["time_GMT"] == "-999")[0])

def read_data_file(filename):
    header = read(filename, data_start=26, data_end=27)
    header["col1"][0] = "year"
    header = header[0].as_void()

    for skiprows in range(85, 100):
        try:
            data_array = np.genfromtxt(filename, skip_header=skiprows, dtype=[int, int, "S10"] + [float]*(len(header)-3))
        except Exception:
            continue
        else:
            break

    data_table = table.Table(data=data_array, names=header)

    return data_table

def convert_timestamp(timestamp):
    try:
        dt = datetime.strptime(timestamp, "%H:%M:%S")
    except ValueError as e:
        if timestamp[-2:] == "60":
            if timestamp[3:5] == "59":
                timestamp_new = f"{int(timestamp[:2]) + 1}:00:00"
            else:
                timestamp_new = timestamp[:2] + f":{int(timestamp[3:5]) + 1}:00"
            dt = datetime.strptime(timestamp_new, "%H:%M:%S")
        else:
            raise e
    time = 60*dt.hour + dt.minute + dt.second/60
    return time

data = []
for j, row in enumerate(master_table):
    # Read each data file and combine them into one big table
    data_files = [f"data/CLT/HyperSAS/{row['Station']}_SAS-H_L3a_{quantity}.txt" for quantity in ["Es", "Lsky", "Lt"]]
    try:
        Es, Lsky, Lt = [read_data_file(filename) for filename in data_files]
    except FileNotFoundError:
        # If a file is missing, go on to the next
        continue
    data_table = table.join(Es, Lsky)
    data_table = table.join(data_table, Lt)

    # Get only the data points within the suggested 3-minute range
    time_start = convert_timestamp(row["time_GMT"])
    time_table = np.array([convert_timestamp(timestamp) for timestamp in data_table["time"]])
    remove_indices = np.where((time_table < time_start) | (time_table > time_start+3))[0]
    data_table.remove_rows(remove_indices)

    # If no data are available, go to the next file
    if len(data_table) == 0:
        continue

    # Calculate median values and replace table with only these
    means = [np.median(data_table[key]) for key in data_table.keys()[3:]]
    mean_data = [*data_table[0]["year", "jd"], row["time_GMT"], *means]
    data_table.add_row(mean_data)
    data_table.remove_rows(np.arange(len(data_table)-1))

    # Finally, load lat/lon
    *_, lon, lat = find_auxiliary_information_seabass(data_files[0])
    data_table.add_column(table.Column(name="Longitude", data=[lon]))
    data_table.add_column(table.Column(name="Latitude", data=[lat]))

    data.append(data_table)

    print(j, row["Station"])

data = table.vstack(data)

Es_keys, Lsky_keys, Lt_keys = get_keys_with_label(data, "Es", "Lsky", "Lt")
for Es_k, Lsky_k, Lt_k in zip(Es_keys, Lsky_keys, Lt_keys):
    wavelength = int(Es_k[2:])

    data[Es_k].unit = u.microwatt / (u.centimeter**2 * u.nanometer)
    data[Es_k] = data[Es_k].to(u.watt / (u.meter**2 * u.nanometer))
    data.rename_column(Es_k, f"Ed_{wavelength}")

    data[Lsky_k].unit = u.microwatt / (u.centimeter**2 * u.nanometer * u.steradian)
    data[Lsky_k] = data[Lsky_k].to(u.watt / (u.meter**2 * u.nanometer * u.steradian))
    data.rename_column(Lsky_k, f"Ls_{wavelength}")

    data[Lt_k].unit = u.microwatt / (u.centimeter**2 * u.nanometer * u.steradian)
    data[Lt_k] = data[Lt_k].to(u.watt / (u.meter**2 * u.nanometer * u.steradian))
    data.rename_column(Lt_k, f"Lt_{wavelength}")

Ed_keys, Ls_keys, Lt_keys = get_keys_with_label(data, "Ed", "Ls", "Lt")
for Ed_k, Ls_k, Lt_k in zip(Ed_keys, Ls_keys, Lt_keys):
    wavelength = int(Ed_k[3:])
    Lw = data[Lt_k] - 0.028 * data[Ls_k]
    Lw.name = f"Lw_{wavelength}"
    data.add_column(Lw)

    R_rs = data[f"Lw_{wavelength}"] / data[Ed_k]
    R_rs.name = f"R_rs_{wavelength}"
    R_rs.unit = 1 / u.steradian
    data.add_column(R_rs)

data.remove_columns(Ls_keys)
data.remove_columns(Lt_keys)

plot_spectra(data, data_label="CLT-S", alpha=0.05)

write_data(data, label="CLT-S")
