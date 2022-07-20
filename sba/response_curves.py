"""
Module with functions for loading and using spectral response functions (SRFs)
"""
import sys
from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
import xarray as xr
from . import bandaveraging as ba, plotting as p


class Band(object):
    def __init__(self, label, wavelengths, response, colour="k"):
        self.label = label
        self.colour = colour

        assert len(wavelengths) == len(response)
        self.wavelengths = wavelengths
        self.response = response

    def __repr__(self):
        return self.label

    def convolve(self, *args, **kwargs):
        result = ba.bandaverage_multi(self.wavelengths, self.response, *args, **kwargs)
        return result


class Sensor(object):
    def __init__(self, name, band_labels, colours, response_wavelengths, responses):
        self.name = name

        assert len(band_labels) == len(colours) == len(response_wavelengths) == len(responses)
        self.bands = [Band(label, wavelengths, response, colour) for label, wavelengths, response, colour in zip(band_labels, response_wavelengths, responses, colours)]

    def __repr__(self):
        return f"{self.name} ({len(self.bands)} bands)"

    def get_band_labels(self):
        return [band.label for band in self.bands]

    def get_band_colours(self):
        return [band.colour for band in self.bands]

    def plot(self, ax=None, saveto=None):
        if ax is None:
            fig, ax = plt.subplots(figsize=(6,2), tight_layout=True)
            if saveto is None:
                saveto = f"spectral_response/SRF_{self.name}.pdf"

            independent = True
        else:
            independent = False

        for band in self.bands:
            ax.plot(band.wavelengths, band.response, label=band.label, c=band.colour)

        ax.set_xticks(np.arange(300, 1500, 100))
        ax.set_xlim(300, 1300)
        ax.set_ylim(0, 1.01)
        if independent:
            ax.set_xlabel("Wavelength [nm]")
            ax.set_ylabel("Relative response")
            ax.set_title(self.name)
        ax.set_yticks([0,0.25,0.5,0.75,1])
        ax.grid(ls="--", color="0.5")

        if saveto is not None:
            plt.savefig(saveto)
            plt.show()
            plt.close()

    def band_average(self, *args, **kwargs):
        result = np.array([band.convolve(*args, **kwargs) for band in self.bands])
        return result

    def boxplot_relative(self, *args, **kwargs):
        p.boxplot_relative(*args, band_labels=self.get_band_labels(), sensor_label=self.name, colours=self.get_band_colours(), **kwargs)

    def boxplot_absolute(self, *args, **kwargs):
        p.boxplot_absolute(*args, band_labels=self.get_band_labels(), sensor_label=self.name, colours=self.get_band_colours(), **kwargs)


def generate_boxcar(center, fwhm, boxcar_wavelength_step=0.1):
    half_width = fwhm / 2.
    wavelengths_in_boxcar = np.arange(center-half_width, center+half_width+boxcar_wavelength_step, boxcar_wavelength_step)
    response = np.ones_like(wavelengths_in_boxcar)
    boxcar_sensor = Sensor("Boxcar", [f"{center:.1f} +- {half_width:.1f} nm"], [""], [wavelengths_in_boxcar], [response])
    return boxcar_sensor


def generate_gaussian(center, fwhm, wavelengths=np.arange(320, 800, 0.1)):
    sigma = fwhm/2.355
    response = np.exp(-(wavelengths-center)**2 / (2 * sigma**2))
    gaussian_sensor = Sensor("Gaussian", [f"{center:.1f} +- {sigma:.1f} nm"], [""], [wavelengths], [response])
    return gaussian_sensor


def read_synthetic_sensor_type():
    sensor_type = sys.argv[2]
    if sensor_type == "gauss":
        sensor_type = "gaussian"
    assert sensor_type in ["boxcar", "gaussian"], f"Unknown sensor type {sensor_type}"
    return sensor_type


def load_synthetic_sensor(sensor_type, center, fwhm, **kwargs):
    if sensor_type == "boxcar":
        func = generate_boxcar
    elif sensor_type in ["gauss", "gaussian"]:
        func = generate_gaussian
    sensor = func(center, fwhm, **kwargs)
    return sensor


def load_OLI():
    band_data_labels = ["CA", "Blue", "Green", "Red", "NIR"]
    band_labels = ["Band 1\nCoastal aerosol", "Band 2\nBlue", "Band 3\nGreen", "Band 4\nRed", "Band 5\nNIR"]
    responses = [np.loadtxt(f"spectral_response/OLI/{band}.txt", skiprows=1, unpack=True, usecols=[0,1]) for band in band_data_labels]
    response_wavelengths = [r[0] for r in responses]
    responses = [r[1] for r in responses]
    colours = ["xkcd:dark blue", "xkcd:blue", "xkcd:green", "xkcd:red", "xkcd:dark brown"]

    OLI = Sensor("OLI", band_labels, colours, response_wavelengths, responses)

    return OLI


def load_ETM_plus():
    bands = [1,2,3,4]
    band_labels = ["Blue", "Green", "Red", "NIR"]
    responses = [np.loadtxt(f"spectral_response/ETM+/spectral_b{x}.dat", skiprows=3, unpack=True) for x in bands]
    response_wavelengths = [r[0] for r in responses]
    responses = [r[1] for r in responses]
    colours = ["b", "g", "r", "k"]

    ETM_plus = Sensor("ETM+", band_labels, colours, response_wavelengths, responses)

    return ETM_plus


def load_VIIRS():
    band_labels = [f"M{j}" for j in np.arange(1,9)]
    wavelengths_viirs, *responses_raw = np.loadtxt("spectral_response/VIIRSN_IDPSv3_RSRs.txt", skiprows=5, unpack=True, usecols=np.arange(9))
    response_wavelengths = [wavelengths_viirs for band in band_labels]
    responses_raw = np.array(responses_raw)
    responses = responses_raw / responses_raw.max(axis=1)[:, np.newaxis]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:cyan", "xkcd:forest green", "xkcd:dark red", "xkcd:dark brown", "k", "k"]

    VIIRS = Sensor("VIIRS", band_labels, colours, response_wavelengths, responses)

    return VIIRS


def load_SeaWiFS():
    band_labels = [f"{wvl} nm" for wvl in [412, 443, 490, 510, 555, 670, 765, 865]]
    wavelengths_seawifs, *responses_raw = np.loadtxt("spectral_response/SeaWiFS_RSRs.txt", skiprows=9, unpack=True)
    response_wavelengths = [wavelengths_seawifs for band in band_labels]
    responses_raw = np.array(responses_raw)
    responses = responses_raw / responses_raw.max(axis=1)[:, np.newaxis]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:cyan", "xkcd:lime green", "xkcd:forest green", "xkcd:dark red", "xkcd:dark brown", "k"]

    SeaWiFS = Sensor("SeaWiFS", band_labels, colours, response_wavelengths, responses)

    return SeaWiFS


def load_Sentinel2(AB="A"):
    band_labels = [f"B{nr}" for nr in range(1,9)] + ["B8a", "B9"]
    wavelengths_msi, *responses = np.loadtxt(f"spectral_response/MSI/MSI_S2{AB}.txt", skiprows=1, unpack=True)
    response_wavelengths = [wavelengths_msi for band in band_labels]
    responses = np.array([responses])[0]
    colours = ["xkcd:blue", "xkcd:cyan", "xkcd:green", "xkcd:red", "xkcd:dark red", "xkcd:dark red", "xkcd:dark red", "xkcd:dark brown", "xkcd:dark brown", "k"]

    Sentinel2 = Sensor(f"Sentinel-2{AB}", band_labels, colours, response_wavelengths, responses)

    return Sentinel2


def load_Sentinel2A():
    return load_Sentinel2("A")


def load_Sentinel2B():
    return load_Sentinel2("B")


def load_MODIS(AT="A"):
    band_labels = [f"{wvl} nm" for wvl in [412, 443, 469, 488, 531, 551, 555, 645, 667, 678, 748, 859, 869, 1240]]
    wavelengths_modis, *responses = np.loadtxt(f"spectral_response/HMODIS{AT}_RSRs.txt", skiprows=8, unpack=True, usecols=np.arange(15))
    responses = np.array([responses])[0]
    response_wavelengths = [wavelengths_modis for band in band_labels]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:blue", "xkcd:cyan", "xkcd:bright green", "xkcd:green", "xkcd:forest green", "xkcd:red", "xkcd:dark red", "xkcd:dark brown", "xkcd:dark brown", "k", "k", "k"]

    name = "MODIS Aqua" if AT == "A" else "MODIS Terra"
    MODIS = Sensor(name, band_labels, colours, response_wavelengths, responses)

    return MODIS


def load_MODISA():
    return load_MODIS("A")


def load_MODIST():
    return load_MODIS("T")


def load_MERIS():
    band_labels = [f"b{j}" for j in range(1,16)]
    wavelengths_meris, *responses = np.loadtxt("spectral_response/MERIS_RSRs.txt", skiprows=5, unpack=True)
    responses = np.array(responses)
    responses[responses < 0] = 0
    response_wavelengths = [wavelengths_meris for band in band_labels]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:cyan", "xkcd:bright green", "xkcd:green", "xkcd:red", "xkcd:dark red", "xkcd:dark red", "xkcd:dark brown", "xkcd:dark brown", "xkcd:dark brown", "k", "k", "k", "k"]

    MERIS = Sensor("MERIS", band_labels, colours, response_wavelengths, responses)

    return MERIS


def load_CZCS():
    band_labels = [f"{wvl} nm" for wvl in [443, 520, 550, 670]]
    wavelengths_czcs, *responses = np.loadtxt("spectral_response/CZCS_RSRs.txt", skiprows=56, unpack=True)
    responses = np.array(responses)
    responses[responses < 0] = 0
    response_wavelengths = [wavelengths_czcs for band in band_labels]
    colours = ["xkcd:dark blue", "xkcd:lime green", "xkcd:forest green", "xkcd:dark red"]

    CZCS = Sensor("CZCS", band_labels, colours, response_wavelengths, responses)

    return CZCS


def load_OLCI(AB="A"):
    filename = "S3A_OL_SRF_20160713_mean_rsr.nc4" if AB == "A" else "S3B_OL_SRF_0_20180109_mean_rsr.nc4"
    band_labels = [f"Oa{nr:02d}" for nr in range(1,22)]
    sr_olci = xr.open_dataset(f"spectral_response/OLCI/{filename}")
    response_wavelengths = sr_olci.mean_spectral_response_function_wavelength.data
    responses = sr_olci.mean_spectral_response_function.data
    colours = ["xkcd:dark violet", "xkcd:violet", "xkcd:purple", "xkcd:light blue", "xkcd:lime green", "xkcd:bright green", "xkcd:bright red", "xkcd:dark red", "xkcd:reddish brown", "xkcd:brown", "xkcd:dark brown", "k", "k", "k", "k", "k", "k", "k", "k", "k", "k"]

    OLCI = Sensor(f"OLCI S3{AB}", band_labels, colours, response_wavelengths, responses)

    return OLCI


def load_OLCIA():
    return load_OLCI("A")


def load_OLCIB():
    return load_OLCI("B")


def load_SPECTACLE():
    files = list(Path("spectral_response/SPECTACLE/").glob("*.npy"))
    devices = [file.stem for file in files]
    bands = [[f"{device} {c}" for c in "RGB"] for device in devices]
    band_labels = [band for sublist in bands for band in sublist]
    colours = p.RGB_OkabeIto * len(devices)

    responses = []
    wavelengths_spectacle = np.arange(390, 701, 1)
    for file in files:
        wavelengths_raw = np.load(file)[0]
        responses_raw = np.load(file)[1:4]
        responses_proc = np.vstack([np.interp(wavelengths_spectacle, wavelengths_raw, response_raw) for response_raw in responses_raw])
        responses.extend(responses_proc)
    responses = np.vstack(responses)
    responses = responses / responses.max(axis=1)[:,np.newaxis]
    response_wavelengths = [wavelengths_spectacle for band in band_labels]

    SPECTACLE = Sensor("SPECTACLE", band_labels, colours, response_wavelengths, responses)

    return SPECTACLE


functions = [load_ETM_plus, load_OLI, load_CZCS, load_SeaWiFS, load_MODISA, load_MODIST, load_MERIS, load_VIIRS, load_Sentinel2A, load_Sentinel2B, load_OLCIA, load_OLCIB, load_SPECTACLE]

from_name = {"etm+": load_ETM_plus, "etm plus": load_ETM_plus, "etm": load_ETM_plus, "landsat7": load_ETM_plus, "l7": load_ETM_plus,
             "oli": load_OLI, "landsat8": load_OLI, "l8": load_OLI,
             "czcs": load_CZCS,
             "seawifs": load_SeaWiFS,
             "modisa": load_MODISA, "modis": load_MODISA,
             "modist": load_MODIST,
             "meris": load_MERIS, "envisat": load_MERIS,
             "viirs": load_VIIRS,
             "sentinel2a": load_Sentinel2A, "sentinel2": load_Sentinel2A, "msi": load_Sentinel2A, "msia": load_Sentinel2A,
             "sentinel2b": load_Sentinel2B, "msib": load_Sentinel2B,
             "olcia": load_OLCIA, "olci": load_OLCIA, "sentinel3a": load_OLCIA, "sentinel3": load_OLCIA,
             "olcib": load_OLCIB, "sentinel3b": load_OLCIB,
             "spectacle": load_SPECTACLE}


def load_selected_sensors(*sensor_names):
    sensor_names = [name.lower() for name in sensor_names]
    functions = [from_name[name] for name in sensor_names]
    sensors = [function() for function in functions]
    return sensors


def load_from_name():
    sensor_names = sys.argv[2:]
    return load_selected_sensors(*sensor_names)


def load_all_sensors():
    sensors = [function() for function in functions]
    return sensors
