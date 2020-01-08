from astropy.io.ascii import read
from numpy import loadtxt, genfromtxt
from pathlib import Path
import sys
from .bandaveraging import split_spectrum

def load_data():
    filename = Path(sys.argv[1])
    data = read(filename)

    label = filename.stem[:-10]

    wavelengths, Ed = split_spectrum(data, "Ed")
    wavelengths, Lw = split_spectrum(data, "Lw")
    wavelengths, R_rs = split_spectrum(data, "R_rs")

    return label, wavelengths, Ed, Lw, R_rs

def write_data(data, label, **kwargs):
    label_lowercase = label.lower()
    data.write(f"data/{label_lowercase}_processed.tab", format="ascii.fast_tab", overwrite=True, **kwargs)


def find_phrase(lines, phrase):
    first_line = [line for line in lines if phrase in line][0]
    return first_line


def find_coordinates_seabass(lines):
    lat_line = find_phrase(lines, "north_latitude")
    latitude = float(lat_line[16:-6])

    lon_line = find_phrase(lines, "west_longitude")
    longitude = float(lon_line[16:-6])

    return longitude, latitude


def fine_datetime_seabass(lines):
    date_line = find_phrase(lines, "end_date")
    date = int(date_line[10:])

    time_line = find_phrase(lines, "start_time")
    time = time_line[12:-6]

    return date, time


def find_auxiliary_information_seabass(file):
    with open(file, "r") as f:
        lines = f.readlines()
        lon, lat = find_coordinates_seabass(lines)
        date, time = fine_datetime_seabass(lines)

    return date, time, lon, lat
