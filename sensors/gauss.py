"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from sba.bandaveraging import split_spectrum, bandaverage_multi
from sba.bandaveraging import load_data

label, wavelengths_data, Ed, Lw, R_rs = load_data()

raise Exception

wavelengths_band = np.arange(320, 800, 0.1)

wavelengths_central = np.arange(360, 780, 5)
FWHMs = np.arange(1, 75, 3)

result_absolute = np.tile(np.nan, [len(FWHMs), len(wavelengths_central)])
result_relative = result_absolute.copy()

wavelengths_interp = np.arange(380, 800.5, 0.5)

def generate_gaussian(x, center, fwhm):
    response = np.exp(-(x-center)**2 / (2 * fwhm**2))
    return response

for i,center in enumerate(wavelengths_central):
    print(f"Central wavelength: {center} nm")
    for j,fwhm in enumerate(FWHMs):
        if center-1.5*fwhm < wavelengths_band[0] or center+1.5*fwhm > wavelengths_band[-1]:
            # Skip combination if the gaussian response does not fall
            # within the data wavelength range up to 3 stds out
            continue
        gaussian_response = generate_gaussian(wavelengths_band, center, fwhm)
        gauss_reflectance_space = bandaverage_multi(wavelengths_band, gaussian_response, wavelengths, R_rs)
        gauss_radiance_space = bandaverage_multi(wavelengths_band, gaussian_response, wavelengths, Lw) / bandaverage_multi(wavelengths_band, gaussian_response, wavelengths, Ed)
        result_absolute[j,i] = np.median((gauss_reflectance_space - gauss_radiance_space))
        result_relative[j,i] = 100*np.median((gauss_reflectance_space - gauss_radiance_space) / gauss_radiance_space)

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
    plt.savefig(f"results/gauss_map_{absrel}.pdf")
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
    plt.savefig(f"results/gauss_contours_{absrel}.pdf")
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
    plt.savefig(f"results/gauss_contours_absolute_{absrel}.pdf")
    plt.show()
