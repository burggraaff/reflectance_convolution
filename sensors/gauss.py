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

wavelengths_band = np.arange(320, 800, 0.1)

wavelengths_central = np.arange(360, 780, 5)
FWHMs = np.concatenate([np.arange(1, 10, 1), np.arange(10, 36, 2), np.arange(36, 104, 4)])

result_absolute = np.tile(np.nan, [len(FWHMs), len(wavelengths_central)])
result_relative = result_absolute.copy()

def generate_gaussian(center, fwhm):
    half_width = fwhm / 2.
    response = np.exp(-(wavelengths_band-center)**2 / (2 * fwhm**2))
    gaussian_sensor = Sensor("Gaussian", [f"{center:.1f} +- {half_width:.1f} nm"], [""], [wavelengths_band], [response])
    return gaussian_sensor

for i,center in enumerate(wavelengths_central):
    print(f"Central wavelength: {center} nm")
    for j,fwhm in enumerate(FWHMs):
        if center-1.5*fwhm < wavelengths_band[0] or center+1.5*fwhm > wavelengths_band[-1]:
            # Skip combination if the gaussian response does not fall
            # within the data wavelength range up to 3 stds out
            continue
        gaussian = generate_gaussian(center, fwhm)
        reflectance_space = gaussian.band_average(wavelengths_data, R_rs)
        radiance_space = gaussian.band_average(wavelengths_data, Lw) / gaussian.band_average(wavelengths_data, Ed)

        difference_absolute, difference_relative = calculate_differences(reflectance_space, radiance_space)

        result_absolute[j,i] = np.median(difference_absolute)
        result_relative[j,i] = np.median(difference_relative)

for gauss_result, absrel, unit in zip([result_absolute, result_relative], ["abs", "rel"], ["sr$^{-1}$", "%"]):
    low, high = np.nanmin(gauss_result), np.nanmax(gauss_result)
    vmin = np.min([low, -high])
    vmax = np.max([-low, high])

    # imshow plots
    im = plt.imshow(gauss_result, vmin=vmin, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], aspect="auto", cmap=plt.cm.seismic)
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title("Difference for Gaussian responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Difference (Refl. space - Rad. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"results/gaussian/map_{absrel}.pdf")
    plt.show()

    # contourf plots
    im = plt.contourf(gauss_result, vmin=vmin, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], levels=np.linspace(vmin, vmax, 11), cmap=plt.cm.seismic)
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title("Difference for Gaussian responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Difference (Refl. space - Rad. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"results/gaussian/contours_{absrel}.pdf")
    plt.show()

    # contourf plot of absolute differences
    im = plt.contourf(np.abs(gauss_result), vmin=0, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], levels=np.linspace(0, vmax, 11))
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title("Absolute difference for Gaussian responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Abs. difference (Refl. space - Rad. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"results/gaussian/contours_absolute_{absrel}.pdf")
    plt.show()
