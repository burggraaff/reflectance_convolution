"""
Plot all spectral response functions (SRFs) separately, and unique ones together
"""

import numpy as np
from matplotlib import pyplot as plt
from sba.response_curves import load_selected_sensors

# Make a combined plot of selected sensors
sensors_selected = load_selected_sensors("etm+", "oli", "spectacle", "seawifs", "viirs", "modis", "msi", "meris", "olci")

fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True, sharey="row", figsize=(7,2), gridspec_kw={"hspace": 0.08, "wspace": 0.03})
for sensor, ax in zip(sensors_selected, axs.T.ravel()):
    sensor.plot(ax)
    ax.tick_params(axis="y", left=False, labelleft=False)

for ax in axs[:-1].ravel():
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
for ax in axs[-1].ravel():
    ax.set_xlabel("Wavelength [nm]")
    ax.set_xticks(np.arange(250, 1500, 250))

for ax in axs[:,0].ravel():
    ax.set_ylabel("RSR")
    ax.tick_params(axis="y", left=True, labelleft=True)

axs[0,0].set_yticks([1/3, 2/3, 1])
axs[0,0].set_yticklabels(["0.33", "0.67", "1.00"])
axs[1,0].set_yticks([1/3, 2/3, 1])
axs[1,0].set_yticklabels(["0.33", "0.67", "1.00"])
axs[2,0].set_yticks([0, 1/3, 2/3, 1])
axs[2,0].set_yticklabels(["0.00", "0.33", "0.67", "1.00"])

for ax in axs.ravel():
    ax.set_ylim(ymin=0)

axs[0,0].set_xlim(320, 1320)

axs[0,0].set_title("ETM+, OLI,\nSPECTACLE")
axs[0,1].set_title("SeaWiFS, Suomi NPP/\nVIIRS, Aqua/MODIS")
axs[0,2].set_title("S2A/MSI, MERIS,\nS3A/OLCI")
plt.savefig("spectral_response/SRF_all.pdf", bbox_inches="tight")
plt.show()
plt.close()
