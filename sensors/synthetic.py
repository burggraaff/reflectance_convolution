"""
Apply synthetic (boxcar and gaussian) spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sba.bandaveraging import calculate_differences
from sba.io import load_data
from sba.response_curves import read_synthetic_sensor_type, load_synthetic_sensor

label, wavelengths_data, Ed, Lw, R_rs = load_data()
sensor_type = read_synthetic_sensor_type()

wavelengths_central = np.arange(330, 800, 5)
FWHMs = np.concatenate([np.arange(1, 30, 1), np.arange(30, 82, 2)])

result_absolute = np.tile(np.nan, [len(FWHMs), len(wavelengths_central)])
result_relative = result_absolute.copy()

for i,center in enumerate(wavelengths_central):
    print(f"Central wavelength: {center} nm")
    for j,fwhm in enumerate(FWHMs):
        boxcar = load_synthetic_sensor(sensor_type, center, fwhm)
        reflectance_space = boxcar.band_average(wavelengths_data, R_rs)
        radiance_space = boxcar.band_average(wavelengths_data, Lw) / boxcar.band_average(wavelengths_data, Ed)

        difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)

        result_absolute[j,i] = np.median(difference_absolute)
        result_relative[j,i] = np.median(difference_relative)

for result, absrel, unit in zip([result_absolute, result_relative], ["abs", "rel"], ["sr$^{-1}$", "%"]):
    low, high = np.nanmin(result), np.nanmax(result)
    vmin = np.min([low, -high])
    vmax = np.max([-low, high])

    # contourf plot
    im = plt.contourf(result, vmin=vmin, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], levels=np.linspace(vmin, vmax, 25), cmap=plt.cm.seismic)
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title(f"Difference for {sensor_type} responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Difference (Rad. space - Refl. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"results/{label}/{label}_{sensor_type}_{absrel}.pdf")
    plt.show()
