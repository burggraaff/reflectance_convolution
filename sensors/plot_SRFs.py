"""
Plot all spectral response functions (SRFs) separately, and unique ones together
"""

import numpy as np
from matplotlib import pyplot as plt
from sba.response_curves import load_all_sensors, load_selected_sensors

# Make a separate plot for all sensors
sensors_all = load_all_sensors()
for sensor in sensors_all:
    print(sensor)
    sensor.plot()

# Make a combined plot of selected sensors
sensors_selected = load_selected_sensors("etm+", "oli", "spectacle", "seawifs", "viirs", "modis", "msi", "meris", "olci")
fig, axs = plt.subplots(nrows=3, ncols=3, sharex=True, sharey=False, tight_layout=True, figsize=(7,3), gridspec_kw={"hspace": 0, "wspace": 0})
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

axs[0,0].set_yticks([0.25, 0.5, 0.75, 1])
axs[1,0].set_yticks([0.25, 0.5, 0.75, 1])
axs[2,0].set_yticks([0, 0.25, 0.5, 0.75, 1])

axs[0,0].set_xlim(320, 1320)

axs[0,0].set_title("ETM+, OLI,\nSPECTACLE")
axs[0,1].set_title("SeaWiFS, VIIRS,\nMODIS Aqua")
axs[0,2].set_title("Sentinel-2A, MERIS,\nSentinel-3A")
plt.savefig("spectral_response/SRF_all.pdf")
plt.show()
plt.close()
