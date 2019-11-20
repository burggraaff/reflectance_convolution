import numpy as np
from matplotlib import pyplot as plt
import bandaveraging as ba

class Sensor(object):
    def __init__(self, name, band_labels, colours, response_wavelengths, responses):
        self.name = name

        assert len(band_labels) == len(colours) == len(response_wavelengths) == len(responses)
        self.band_labels = band_labels
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
        ax.set_xlim(380, 900)
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
