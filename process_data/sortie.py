import numpy as np
from astropy import table
from astropy import units as u
from pathlib import Path
from sba.plotting import plot_spectra, map_data
from sba.io import read, write_data, find_auxiliary_information_seabass
from sba.data_processing import get_keys_with_label, remove_negative_R_rs

folder = Path("data/SORTIE/")
files = list(folder.glob("SORTIE*.dat"))

tabs = []
for file in files:
    wavelengths, Lw, Ed = np.loadtxt(file, delimiter=",", skiprows=86, unpack=True, usecols=[0,8,11])

    date, time, lon, lat = find_auxiliary_information_seabass(file)

    cols = ["Date", "Time", "Latitude", "Longitude"] + [f"Lw_{wvl:.2f}" for wvl in wavelengths] + [f"Ed_{wvl:.2f}" for wvl in wavelengths]
    dtype = [int, "S8", float, float] + 2 * [float for wvl in wavelengths]
    tab = table.Table(rows=[[date, time, lat, lon, *Lw, *Ed]], names=cols, dtype=dtype)
    tabs.append(tab)

data = table.vstack(tabs)
