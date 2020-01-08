"""
Module with functions for plotting
"""
from matplotlib import pyplot as plt
from mpl_toolkits.basemap import Basemap
from .bandaveraging import calculate_median_and_errors
from pathlib import Path
import numpy as np

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

    medians, lower_error, upper_error = calculate_median_and_errors(differences)
    for band, med, low, up in zip(band_labels, medians, lower_error, upper_error):
        print(f"{sensor_label} {band} band: {med:+.2f} (+{up:.2f}, -{low:.2f}) %")

    make_boxplot(differences, label="rel", unit="%", sensor_label=sensor_label, band_labels=band_labels, **kwargs)
    double_boxplot(differences, label="rel", unit="%", sensor_label=sensor_label, band_labels=band_labels, **kwargs)

def boxplot_absolute(differences, band_labels=None, sensor_label="", scaling_exponent=6, **kwargs):
    if band_labels is None:
        band_labels = [""] * len(differences)

    differences_scaled = differences * 10**scaling_exponent
    unit = "$10^{-" + f"{scaling_exponent}" + "}$ sr$^{-1}$"

    medians, lower_error, upper_error = calculate_median_and_errors(differences_scaled)
    for band, med, low, up in zip(band_labels, medians, lower_error, upper_error):
        print(f"{sensor_label} {band} band: {med:+.2f} (+{up:.2f}, -{low:.2f}) x 10^-6 sr^-1")

    make_boxplot(differences_scaled, label="abs", unit=unit, sensor_label=sensor_label, band_labels=band_labels, **kwargs)
    double_boxplot(differences_scaled, label="abs", unit=unit, sensor_label=sensor_label, band_labels=band_labels, **kwargs)

def plot_spectra(data, data_label="", alpha=0.1):
    # Plot all Es, Lw, R_rs spectra
    wavelengths = [float(key[3:]) for key in data.keys() if "Lw" in key]

    fig, axs = plt.subplots(nrows=3, ncols=1, sharex=True, tight_layout=True, gridspec_kw={"wspace":0, "hspace":0}, figsize=(5,7))

    for row in data:
        spec_Ed = [row[key] for key in data.keys() if "Ed" in key]
        spec_Lw = [row[key] for key in data.keys() if "Lw" in key]
        spec_R_rs = [row[key] for key in data.keys() if "R_rs" in key]

        for ax, spec in zip(axs.ravel(), [spec_Ed, spec_Lw, spec_R_rs]):
            ax.plot(wavelengths, spec, c="k", alpha=alpha, zorder=1)

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
