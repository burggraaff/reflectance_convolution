"""
Module with functions for plotting
"""

from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.axes_grid1 import make_axes_locatable, ImageGrid
from .bandaveraging import calculate_median_and_errors
from .data_processing import split_spectrum
from pathlib import Path
import numpy as np
import warnings


RrsR = r"$\bar{R}_{rs}^R$"
RrsL = r"$\bar{R}_{rs}^L$"


def double_boxplot(data, label="", unit="", sensor_label="", data_label="", band_labels=None, colours=None):
    save_folder = Path(f"results/{data_label}")
    save_folder.mkdir(exist_ok=True)

    if band_labels is None:
        band_labels = [""] * len(data)
    if colours is None:
        colours = ["k"] * len(data)

    fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 1+0.3*data.shape[0]), tight_layout=True, sharey=True, gridspec_kw={"hspace": 0, "wspace": 0})
    bplot_main = axs[0].boxplot(data.T, vert=False, showfliers=False, whis=[5,95], patch_artist=True, labels=band_labels)
    for patch, colour in zip(bplot_main["boxes"], colours):
        patch.set_facecolor(colour)
    bplot_zoom = axs[1].boxplot(data.T, vert=False, showfliers=False, whis=[5,95], patch_artist=True, labels=band_labels)
    for patch, colour in zip(bplot_zoom["boxes"], colours):
        patch.set_facecolor(colour)

    for ax in axs:
        ax.set_xlabel(f"Difference [{unit}]")
        ax.grid(ls="--", color="0.5")
        ax.axvline(0, ls="--", color="k")

    axs[0].set_title(f"Data: {data_label} ($N$ = {data.shape[1]}) ; Sensor: {sensor_label}")

    xlim = axs[0].get_xlim()
    if xlim[0] > 0:
        axs[0].set_xlim(xmin=0)
    elif xlim[-1] < 0:
        axs[0].set_xlim(xmax=0)

    axs[1].tick_params(axis="y", left=False, labelleft=False)

    if label == "rel": # relative plot
        axs[1].set_xlim(-5.1, 5.1)
        axs[1].set_xticks(np.arange(-5, 6, 1))
    else: # absolute plot
        axs[1].set_xlim(-30.1, 30.1)

    plt.savefig(f"{save_folder}/{data_label}_{sensor_label}_{label}_double.pdf")
    plt.show()
    plt.close()


def make_boxplot(data, label="", unit="", sensor_label="", data_label="", band_labels=None, colours=None):
    save_folder = Path(f"results/{data_label}")
    save_folder.mkdir(exist_ok=True)

    if band_labels is None:
        band_labels = [""] * len(data)
    if colours is None:
        colours = ["k"] * len(data)

    plt.figure(figsize=(5, 1+0.3*data.shape[0]), tight_layout=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bplot = plt.boxplot(data.T, vert=False, showfliers=False, whis=[5,95], patch_artist=True, labels=band_labels)
    for patch, colour in zip(bplot["boxes"], colours):
        patch.set_facecolor(colour)
    plt.xlabel(f"Difference [{unit}]")
    plt.title(f"Data: {data_label} ($N$ = {data.shape[1]})  ; Sensor: {sensor_label}")
    plt.grid(ls="--", color="0.5")
    plt.axvline(0, ls="--", color="k")

    xlim = plt.xlim()
    if xlim[0] > 0:
        plt.xlim(xmin=0)
    elif xlim[-1] < 0:
        plt.xlim(xmax=0)

    plt.tight_layout()
    plt.savefig(f"{save_folder}/{data_label}_{sensor_label}_{label}.pdf")
    plt.show()
    plt.close()


def boxplot_relative(differences, band_labels=None, sensor_label="", **kwargs):
    if band_labels is None:
        band_labels = [""] * len(differences)

    make_boxplot(differences, label="rel", unit="%", sensor_label=sensor_label, band_labels=band_labels, **kwargs)


def boxplot_absolute(differences, band_labels=None, sensor_label="", scaling_exponent=6, **kwargs):
    if band_labels is None:
        band_labels = [""] * len(differences)

    differences_scaled = differences * 10**scaling_exponent
    unit = "$10^{-" + f"{scaling_exponent}" + "}$ sr$^{-1}$"

    make_boxplot(differences_scaled, label="abs", unit=unit, sensor_label=sensor_label, band_labels=band_labels, **kwargs)


def plot_spectra(data, data_label="", alpha=0.1):
    # Plot all Es, Lw, R_rs spectra
    wavelengths, Ed = split_spectrum(data, "Ed")
    wavelengths, Lw = split_spectrum(data, "Lw")
    wavelengths, R_rs = split_spectrum(data, "R_rs")

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,7))

    for ax, spectrum in zip(axs.ravel(), [Ed, Lw, R_rs]):
        ax.plot(wavelengths, spectrum.T, c="k", alpha=alpha, zorder=1)

    for ax, label in zip(axs.ravel(), ["$E_d$", "$L_w$", "$R_{rs}$"]):
        ax.set_ylabel(label)
        ax.grid(ls="--", zorder=0)

    axs[-1].set_xlabel("Wavelength [nm]")
    axs[-1].set_xlim(wavelengths[0], wavelengths[-1])

    axs[0].set_title(f"{data_label} spectra ({len(data)})")
    plt.savefig(f"data/plots/spectra_{data_label}.pdf")
    plt.show()
    plt.close()


def map_data(data, data_label="", projection="moll", figsize=(10, 6), parallels=np.arange(-90, 95, 15), meridians=np.arange(-180, 180, 30), **kwargs):
    # Plot map of observations
    plt.figure(figsize=figsize, tight_layout=True)

    m = Basemap(projection=projection, **kwargs)
    m.fillcontinents(color="#FFDDCC", lake_color='#DDEEFF')
    m.drawmapboundary(fill_color="#DDEEFF")
    m.drawcoastlines()

    m.drawparallels(parallels, labels=[1,1,0,0])
    m.drawmeridians(meridians, labels=[0,0,0,1])

    m.scatter(data["Longitude"], data["Latitude"], latlon=True, c="r", edgecolors="k", s=60, zorder=10)

    plt.title(f"Locations of {data_label} data ($N = {len(data)}$)")

    plt.savefig(f"data/plots/map_{data_label}.pdf")
    plt.show()
    plt.close()


def synthetic_sensor_contourf(wavelengths_central, FWHMs, result, sensor_type="", absrel="", label="", quantity=""):
    if absrel == "rel":
        unit = "%"
    elif absrel == "abs":
        unit = "sr$^{-1}$"
    else:
        unit = ""

    low, high = np.nanmin(result), np.nanmax(result)
    vmin = np.min([low, -high])
    vmax = np.max([-low, high])

    # contourf plot
    im = plt.contourf(result, vmin=vmin, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], levels=np.linspace(vmin, vmax, 25), cmap=plt.cm.seismic)
    plt.xlabel("Central wavelength [nm]")
    plt.ylabel("FWHM [nm]")
    plt.title(f"Difference for {sensor_type}s: {quantity}")
    divider = make_axes_locatable(plt.gca())
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    cax.set_ylabel(f"Difference ({RrsL} - {RrsR}, {unit})")
    plt.tight_layout()
    plt.savefig(f"results/{label}/{label}_{sensor_type}_{quantity}_{absrel}.pdf")
    plt.show()


def synthetic_sensor_contourf_combined(wavelengths_central, FWHMs, results, sensor_type="", absrel="", label="", quantities=None):
    if absrel == "rel":
        unit = "%"
    elif absrel == "abs":
        unit = "$10^{-6}$ sr$^{-1}$"
    else:
        unit = ""

    if quantities is None:
        quantities = [""] * len(results)

    center = len(results)//2

    low, high = np.nanmin(results), np.nanmax(results)
    vmin = np.floor(np.min([low, -high]))
    vmax = np.ceil(np.max([-low, high]))

    fig = plt.figure(figsize=(8, 3.5))
    grid = ImageGrid(fig, 111, nrows_ncols=(1,3), axes_pad=0.15, share_all=True, aspect=False, cbar_location="right", cbar_mode="single", cbar_size="7%", cbar_pad=0.15)

    for ax, result, quantity in zip(grid, results, quantities):
        im = ax.contourf(result, vmin=vmin, vmax=vmax, origin="lower", extent=[wavelengths_central[0], wavelengths_central[-1], FWHMs[0], FWHMs[-1]], levels=np.linspace(vmin, vmax, 25), cmap=plt.cm.seismic)
        ax.set_title(quantity)
        ax.grid(ls="--")
        ax.set_xticks(np.arange(300, 1500, 100))
        ax.set_xlim(wavelengths_central[0], 2*wavelengths_central[-1]-wavelengths_central[-2])

    grid[0].set_ylabel("FWHM [nm]")
    grid[center].set_xlabel("Central wavelength [nm]")
    for ax in grid[1:]:
        ax.tick_params(axis="y", left=False, labelleft=False)
    cax = ax.cax.colorbar(im)
    cax.set_label_text(f"Difference ({RrsL} - {RrsR}, {unit})")

    plt.savefig(f"results/{label}/{label}_{sensor_type}_{absrel}.pdf")
    plt.show()
