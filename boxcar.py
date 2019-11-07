"""
Generate boxcar and gaussian spectral response functions
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io.ascii import read
from astropy import table
from mpl_toolkits.axes_grid1 import make_axes_locatable
from bandaveraging import split_spectrum, bandaverage_multi

wavelengths_central = np.arange(330, 800, 5)
FWHMs = np.arange(1, 75, 2)

result_absolute = np.tile(np.nan, [len(FWHMs), len(wavelengths_central)])
result_relative = result_absolute.copy()

def generate_boxcar(center, fwhm, boxcar_wavelength_step = 0.1):
    wavelengths_in_boxcar = np.arange(center-fwhm/2., center+fwhm/2.+boxcar_wavelength_step, boxcar_wavelength_step)
    response = np.ones_like(wavelengths_in_boxcar)
    return wavelengths_in_boxcar, response

data_norcohab = read("data/norcohab_processed.tab")
data_archemhab = read("data/archemhab_processed.tab")

data_all = table.vstack([data_norcohab, data_archemhab])

wavelengths, Ed = split_spectrum(data_all, "Ed")
wavelengths, Lw = split_spectrum(data_all, "Lw")
wavelengths, R_rs = split_spectrum(data_all, "R_rs")

for i,center in enumerate(wavelengths_central):
    print(f"Central wavelength: {center} nm")
    for j,fwhm in enumerate(FWHMs):
        if center-fwhm/2 < wavelengths[0] or center+fwhm/2 > wavelengths[-1]:
            # Skip combination if the boxcar response does not fall entirely
            # within the data wavelength range
            continue
        wavelengths_boxcar, response_boxcar = generate_boxcar(center, fwhm)
        boxcar_reflectance_space = bandaverage_multi(wavelengths_boxcar, response_boxcar, wavelengths, R_rs)
        boxcar_radiance_space = bandaverage_multi(wavelengths_boxcar, response_boxcar, wavelengths, Lw) / bandaverage_multi(wavelengths_boxcar, response_boxcar, wavelengths, Ed)
        result_absolute[j,i] = np.median((boxcar_reflectance_space - boxcar_radiance_space))
        result_relative[j,i] = 100*np.median((boxcar_reflectance_space - boxcar_radiance_space) / boxcar_radiance_space)

for boxcar_result, absrel, unit in zip([result_absolute, result_relative], ["abs", "rel"], ["sr$^{-1}$", "%"]):
    low, high = np.nanmin(boxcar_result), np.nanmax(boxcar_result)
    vmin = np.min([low, -high])
    vmax = np.max([-low, high])

    # imshow plot
    im = plt.imshow(boxcar_result, vmin=vmin, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], aspect="auto", cmap=plt.cm.seismic)
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title("Difference for boxcar responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Difference (Rad. space - Refl. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"boxcar_map_{absrel}.pdf")
    plt.show()

    # contourf plot
    im = plt.contourf(boxcar_result, vmin=vmin, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], levels=np.linspace(vmin, vmax, 25), cmap=plt.cm.seismic)
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title("Difference for boxcar responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Difference (Rad. space - Refl. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"boxcar_contours_{absrel}.pdf")
    plt.show()

    # contourf plot of absolute differences
    im = plt.contourf(np.abs(boxcar_result), vmin=0, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], levels=np.linspace(0, vmax, 25))
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title("Absolute difference for boxcar responses")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Abs. difference (Rad. space - Refl. space, {unit})")
    plt.tight_layout()
    plt.savefig(f"boxcar_contours_absolute_{absrel}.pdf")
    plt.show()
