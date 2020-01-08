"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sba.bandaveraging import calculate_differences
from sba.io import load_data
from sba.response_curves import Sensor

label, wavelengths_data, Ed, Lw, R_rs = load_data()

wavelengths_central = np.arange(330, 800, 5)
FWHMs = np.concatenate([np.arange(1, 30, 1), np.arange(30, 82, 2)])

result_absolute = np.tile(np.nan, [len(FWHMs), len(wavelengths_central)])
result_relative = result_absolute.copy()

def generate_boxcar(center, fwhm, boxcar_wavelength_step = 0.1):
    half_width = fwhm / 2.
    wavelengths_in_boxcar = np.arange(center-half_width, center+half_width+boxcar_wavelength_step, boxcar_wavelength_step)
    response = np.ones_like(wavelengths_in_boxcar)
    boxcar_sensor = Sensor("Boxcar", [f"{center:.1f} +- {half_width:.1f} nm"], [""], [wavelengths_in_boxcar], [response])
    return boxcar_sensor

for i,center in enumerate(wavelengths_central):
    print(f"Central wavelength: {center} nm")
    for j,fwhm in enumerate(FWHMs):
        half_width = fwhm / 2.
        if center-half_width < wavelengths_data[0] or center+half_width > wavelengths_data[-1]:
            # Skip if the boxcar response does not fall entirely within the data wavelength range
            continue
        boxcar = generate_boxcar(center, fwhm)
        reflectance_space = boxcar.band_average(wavelengths_data, R_rs)
        radiance_space = boxcar.band_average(wavelengths_data, Lw) / boxcar.band_average(wavelengths_data, Ed)

        difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)

        result_absolute[j,i] = np.median(difference_absolute)
        result_relative[j,i] = np.median(difference_relative)

for result, absrel, unit in zip([result_absolute, result_relative], ["abs", "rel"], ["sr$^{-1}$", "%"]):
    low, high = np.nanmin(result), np.nanmax(result)
    vmin = np.min([low, -high])
    vmax = np.max([-low, high])

    # imshow plot
    im = plt.imshow(result, vmin=vmin, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], aspect="auto", cmap=plt.cm.seismic)
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title("Difference for boxcar responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Difference (Rad. space - Refl. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"results/boxcar/map_{absrel}.pdf")
    plt.show()

    # contourf plot
    im = plt.contourf(result, vmin=vmin, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], levels=np.linspace(vmin, vmax, 25), cmap=plt.cm.seismic)
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title("Difference for boxcar responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Difference (Rad. space - Refl. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"results/boxcar/contours_{absrel}.pdf")
    plt.show()

    # contourf plot of absolute differences
    im = plt.contourf(np.abs(result), vmin=0, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], levels=np.linspace(0, vmax, 25))
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title("Absolute difference for boxcar responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Abs. difference (Rad. space - Refl. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"results/boxcar/contours_absolute_{absrel}.pdf")
    plt.show()
