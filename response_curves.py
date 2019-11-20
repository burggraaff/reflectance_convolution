import xarray as xr
import numpy as np
from matplotlib import pyplot as plt
import bandaveraging as ba
from pathlib import Path

class Sensor(object):
    def __init__(self, name, band_labels, colours, response_wavelengths, responses):
        self.name = name

        assert len(band_labels) == len(colours) == len(response_wavelengths) == len(responses)
        self.band_labels = band_labels
        self.sensor_band_labels = [self.name + " " + label for label in self.band_labels]
        self.colours = colours
        self.responses = responses
        self.wavelengths = response_wavelengths

    def __repr__(self):
        return self.name + " (bands: " + ", ".join(self.band_labels) + ")"

    def plot(self, ax=None, saveto=None):
        if ax is None:
            fig, ax = plt.subplots()
            if saveto is None:
                saveto = f"results/{self.name}_bands.pdf"

        for wavelengths, response, label, colour in zip(self.wavelengths, self.responses, self.band_labels, self.colours):
            ax.plot(wavelengths, response, label=label, c=colour)
        ax.set_xlim(380, 925)
        ax.set_ylim(0, 1.01)
        ax.set_xlabel("Wavelength [nm]")
        ax.set_ylabel("Relative response")
        ax.set_title(self.name)
        ax.legend(loc="best")
        ax.grid(ls="--", color="0.5")

        if saveto is not None:
            plt.savefig(saveto)
            plt.show()
            plt.close()

    def band_average(self, *args, **kwargs):
        result = np.array([ba.bandaverage_multi(wvl, response, *args, **kwargs) for wvl, response in zip(self.wavelengths, self.responses)])
        return result

    def boxplot_relative(self, *args, **kwargs):
        ba.boxplot_relative(*args, band_labels=self.band_labels, sensor_label=self.name, colours=self.colours, **kwargs)

    def boxplot_absolute(self, *args, **kwargs):
        ba.boxplot_absolute(*args, band_labels=self.band_labels, sensor_label=self.name, colours=self.colours, **kwargs)


def load_OLI():
    band_labels = ["CA", "Blue", "Green", "Red", "NIR"]
    responses = [np.loadtxt(f"spectral_response/OLI/{band}.txt", skiprows=1, unpack=True, usecols=[0,1]) for band in band_labels]
    response_wavelengths = [r[0] for r in responses]
    responses = [r[1] for r in responses]
    colours = ["xkcd:dark blue", "xkcd:blue", "xkcd:green", "xkcd:red", "xkcd:dark brown"]

    OLI = Sensor("OLI", band_labels, colours, response_wavelengths, responses)

    return OLI

def load_ETM_plus():
    bands = [1,2,3,4]
    band_labels = [f"ETM+\nband {band}" for band in bands]
    responses = [np.loadtxt(f"spectral_response/ETM+/spectral_b{x}.dat", skiprows=3, unpack=True) for x in bands]
    response_wavelengths = [r[0] for r in responses]
    responses = [r[1] for r in responses]
    colours = ["b", "g", "r", "xkcd:dark red"]

    ETM_plus = Sensor("ETM+", band_labels, colours, response_wavelengths, responses)

    return ETM_plus

def load_VIIRS():
    band_labels = [f"M{j}" for j in np.arange(1,7)]
    wavelengths_viirs, *responses_raw = np.loadtxt("spectral_response/VIIRSN_IDPSv3_RSRs.txt", skiprows=5, unpack=True, usecols=np.arange(7))
    response_wavelengths = [wavelengths_viirs for band in band_labels]
    responses_raw = np.array(responses_raw)
    responses = responses_raw / responses_raw.max(axis=1)[:, np.newaxis]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:cyan", "xkcd:forest green", "xkcd:dark red", "k"]

    VIIRS = Sensor("VIIRS", band_labels, colours, response_wavelengths, responses)

    return VIIRS

def load_SeaWiFS():
    band_labels = [f"{wvl} nm" for wvl in [412, 443, 490, 510, 555, 670, 765]]
    wavelengths_seawifs, *responses_raw = np.loadtxt("spectral_response/SeaWiFS_RSRs.txt", skiprows=9, unpack=True, usecols=np.arange(8))
    response_wavelengths = [wavelengths_seawifs for band in band_labels]
    responses_raw = np.array(responses_raw)
    responses = responses_raw / responses_raw.max(axis=1)[:, np.newaxis]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:cyan", "xkcd:lime green", "xkcd:forest green", "xkcd:dark red", "k"]

    SeaWiFS = Sensor("SeaWiFS", band_labels, colours, response_wavelengths, responses)

    return SeaWiFS

def load_Sentinel2A():
    band_labels = [f"B{nr}" for nr in range(1,8)]
    wavelengths_msi, *responses = np.loadtxt("spectral_response/MSI/MSI_S2A.txt", skiprows=1, unpack=True)
    response_wavelengths = [wavelengths_msi for band in band_labels]
    responses = np.array([responses])[0]
    colours = ["xkcd:blue", "xkcd:cyan", "xkcd:green", "xkcd:red", "xkcd:dark red", "xkcd:dark brown", "k"]

    Sentinel2A = Sensor("Sentinel-2A", band_labels, colours, response_wavelengths, responses)

    return Sentinel2A

def load_Sentinel2B():
    band_labels = [f"B{nr}" for nr in range(1,8)]
    wavelengths_msi, *responses = np.loadtxt("spectral_response/MSI/MSI_S2A.txt", skiprows=1, unpack=True)
    response_wavelengths = [wavelengths_msi for band in band_labels]
    responses = np.array([responses])[0]
    colours = ["xkcd:blue", "xkcd:cyan", "xkcd:green", "xkcd:red", "xkcd:dark red", "xkcd:dark brown", "k"]

    Sentinel2B = Sensor("Sentinel-2B", band_labels, colours, response_wavelengths, responses)

    return Sentinel2B

def load_MODISA():
    band_labels = [f"{wvl} nm" for wvl in [412, 443, 469, 488, 531, 551, 555, 645, 667, 678, 748]]
    wavelengths_modis, *responses = np.loadtxt("spectral_response/HMODISA_RSRs.txt", skiprows=8, unpack=True, usecols=np.arange(12))
    responses = np.array([responses])[0]
    ind = np.where(wavelengths_modis <= 900)[0]
    wavelengths_modis = wavelengths_modis[ind]
    response_wavelengths = [wavelengths_modis for band in band_labels]
    responses = responses[:, ind]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:blue", "xkcd:cyan", "xkcd:bright green", "xkcd:green", "xkcd:forest green", "xkcd:red", "xkcd:dark red", "xkcd:dark brown", "k"]

    MODISA = Sensor("MODIS Aqua", band_labels, colours, response_wavelengths, responses)

    return MODISA

def load_MODIST():
    band_labels = [f"{wvl} nm" for wvl in [412, 443, 469, 488, 531, 551, 555, 645, 667, 678, 748]]
    wavelengths_modis, *responses = np.loadtxt("spectral_response/HMODIST_RSRs.txt", skiprows=8, unpack=True, usecols=np.arange(12))
    responses = np.array([responses])[0]
    ind = np.where(wavelengths_modis <= 900)[0]
    wavelengths_modis = wavelengths_modis[ind]
    response_wavelengths = [wavelengths_modis for band in band_labels]
    responses = responses[:, ind]
    colours = ["xkcd:dark purple", "xkcd:dark blue", "xkcd:blue", "xkcd:cyan", "xkcd:bright green", "xkcd:green", "xkcd:forest green", "xkcd:red", "xkcd:dark red", "xkcd:dark brown", "k"]

    MODISA = Sensor("MODIS Terra", band_labels, colours, response_wavelengths, responses)

    return MODISA

def load_CZCS():
    band_labels = [f"{wvl} nm" for wvl in [443, 520, 550, 670]]
    wavelengths_czcs, *responses = np.loadtxt("spectral_response/CZCS_RSRs.txt", skiprows=56, unpack=True)
    responses = np.array(responses)
    responses[responses < 0] = 0
    response_wavelengths = [wavelengths_czcs for band in band_labels]
    colours = ["xkcd:dark blue", "xkcd:lime green", "xkcd:forest green", "xkcd:dark red"]

    CZCS = Sensor("CZCS", band_labels, colours, response_wavelengths, responses)

    return CZCS

def load_OLCIA():
    band_labels = [f"Oa{nr:02d}" for nr in range(1,17)]
    sr_olci = xr.open_dataset("spectral_response/OLCI/S3A_OL_SRF_20160713_mean_rsr.nc4")
    response_wavelengths = sr_olci.mean_spectral_response_function_wavelength.data[:16]
    responses = sr_olci.mean_spectral_response_function.data[:16]
    colours = ["xkcd:dark violet", "xkcd:violet", "xkcd:purple", "xkcd:light blue", "xkcd:lime green", "xkcd:bright green", "xkcd:bright red", "xkcd:dark red", "xkcd:reddish brown", "xkcd:brown", "xkcd:dark brown", "k", "k", "k", "k", "k"]

    OLCIA = Sensor("OLCI S3A", band_labels, colours, response_wavelengths, responses)

    return OLCIA

def load_OLCIB():
    band_labels = [f"Oa{nr:02d}" for nr in range(1,17)]
    sr_olci = xr.open_dataset("spectral_response/OLCI/S3B_OL_SRF_0_20180109_mean_rsr.nc4")
    response_wavelengths = sr_olci.mean_spectral_response_function_wavelength.data[:16]
    responses = sr_olci.mean_spectral_response_function.data[:16]
    colours = ["xkcd:dark violet", "xkcd:violet", "xkcd:purple", "xkcd:light blue", "xkcd:lime green", "xkcd:bright green", "xkcd:bright red", "xkcd:dark red", "xkcd:reddish brown", "xkcd:brown", "xkcd:dark brown", "k", "k", "k", "k", "k"]

    OLCIB = Sensor("OLCI S3B", band_labels, colours, response_wavelengths, responses)

    return OLCIB

def load_SPECTACLE():
    files = list(Path("spectral_response/SPECTACLE/").glob("*.npy"))
    devices = [file.stem for file in files]
    bands = [[f"{device} {c}" for c in "RGB"] for device in devices]
    band_labels = [band for sublist in bands for band in sublist]
    colours = [*"rgb"] * len(devices)

    responses = []
    wavelengths_spectacle = np.arange(390, 701, 1)
    for file in files:
        wavelengths_raw = np.load(file)[0]
        responses_raw = np.load(file)[1:4]
        responses_proc = np.vstack([np.interp(wavelengths_spectacle, wavelengths_raw, response_raw) for response_raw in responses_raw])
        responses.extend(responses_proc)
    responses = np.vstack(responses)
    response_wavelengths = [wavelengths_spectacle for band in band_labels]

    SPECTACLE = Sensor("SPECTACLE", band_labels, colours, response_wavelengths, responses)

    return SPECTACLE
