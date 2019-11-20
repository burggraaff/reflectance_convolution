import numpy as np

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

def load_OLI():
    band_labels = ["CA", "Blue", "Green", "Red", "NIR"]
    responses = [np.loadtxt(f"spectral_response/OLI/{band}.txt", skiprows=1, unpack=True, usecols=[0,1]) for band in band_labels]
    response_wavelengths = [r[0] for r in responses]
    responses = [r[1] for r in responses]
    colours = ["xkcd:dark blue", "xkcd:blue", "xkcd:green", "xkcd:red", "xkcd:dark brown"]

    OLI = Sensor("OLI", band_labels, colours, response_wavelengths, responses)

    return OLI
