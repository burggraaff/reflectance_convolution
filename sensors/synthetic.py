"""
Apply synthetic (boxcar and gaussian) spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from sba.bandaveraging import calculate_differences
from sba.io import load_data
from sba.response_curves import read_synthetic_sensor_type, load_synthetic_sensor
from sba.plotting import synthetic_sensor_contourf, synthetic_sensor_contourf_combined

label, wavelengths_data, Ed, Lw, R_rs = load_data()
sensor_type = read_synthetic_sensor_type()

wavelengths_central = np.arange(330, 810, 1)
# FWHMs = np.arange(6, 66, 1)
FWHMs = np.array([30,40,50])

# Arrays for storing results - per FWHM, per central wavelength, absolute/relative
medians = np.tile(np.nan, [2, len(FWHMs), len(wavelengths_central)])  # Median
perc5s = np.copy(medians)  # 5th percentile
perc95s = np.copy(medians)  # 95th percentile

for i,center in enumerate(wavelengths_central):
    print(f"Central wavelength: {center} nm")
    for j,fwhm in enumerate(FWHMs):
        boxcar = load_synthetic_sensor(sensor_type, center, fwhm)
        reflectance_space = boxcar.band_average(wavelengths_data, R_rs)
        radiance_space = boxcar.band_average(wavelengths_data, Lw) / boxcar.band_average(wavelengths_data, Ed)

        difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)
        difference_combined = np.array([difference_absolute, difference_relative])[:,0]

        medians[:,j,i] = np.median(difference_combined, axis=1)
        perc5s[:,j,i] = np.percentile(difference_combined, 5, axis=1)
        perc95s[:,j,i] = np.percentile(difference_combined, 95, axis=1)

results_stacked = np.stack([perc5s, medians, perc95s])
results_absrel = np.moveaxis(results_stacked, 1, 0)
results_absrel[0] *= 1e6  # Convert to 10^-6 sr^-1
quantities = ["P5", "Median", "P95"]

results_thesis = np.stack([wavelengths_central, medians[1,2]])
np.save("thesis.npy", results_thesis)

for results_combined, absrel in zip(results_absrel, ["abs", "rel"]):
    synthetic_sensor_contourf_combined(wavelengths_central, FWHMs, results_combined, sensor_type=sensor_type, absrel=absrel, label=label, quantities=quantities)

inds = np.searchsorted(FWHMs, [6, 10, 15, 20, 30, 40, 50])
line_labels = [f"{FWHM} nm" for FWHM in FWHMs[inds]]
results_inds = results_absrel[1,:,inds]
results_inds = np.moveaxis(results_inds, 0, 1)

fig = plt.figure(figsize=(7,2))
for p5, med, p95, line_label in zip(*results_inds, line_labels):
    plt.plot(wavelengths_central, med.T, label=line_label)
    plt.fill_between(wavelengths_central, p5, p95, alpha=0.35)
plt.xlim(330, 800)
plt.xlabel("Central wavelength [nm]")
plt.ylabel(r"$\Delta \bar R_{rs}$ [%]")
plt.title(f"Convolution error in {label} data with synthetic {sensor_type} bands")
plt.grid(ls="--")
legend = plt.legend(title="FWHM", ncol=1, loc="center right", bbox_to_anchor=(1.2, 0.5))
renderer = fig.canvas.get_renderer()
shift = max([t.get_window_extent(renderer).width for t in legend.get_texts()])
for t in legend.get_texts():
    t.set_ha("right")
    t.set_position((shift,0))
plt.setp(legend.get_lines(), linewidth=2)
plt.savefig(f"results/{label}/{label}_{sensor_type}_line.pdf", bbox_inches="tight")
plt.show()
