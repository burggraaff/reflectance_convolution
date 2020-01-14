#
## Ha+17, Sentinel-2
## Input: Hyperspectral data
## Load S2A, convolve both ways, calculate Chla both ways
#Chla = 0.80 × exp(0.35 × B3/B4)
#
#
## OC2-4
#https://oceancolor.gsfc.nasa.gov/atbd/chlor_a/
#
import numpy as np

from .response_curves import load_SeaWiFS, load_MERIS, load_MODISA, load_VIIRS, load_CZCS


def OCx(R_rs_blue, R_rs_green, a):
    powers = np.arange(1,5)
    log_ratio = np.log10(R_rs_blue / R_rs_green)
    sum_evaluated = np.sum(a[1:] * log_ratio[:,np.newaxis]**powers, axis=1)
    chla = 10**(a[0] + sum_evaluated)
    return chla


def OC4(wavelengths, Ed, Lw, R_rs):
    sensor = load_SeaWiFS()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    a = np.array([0.3272, -2.9940, 2.7218, -1.2259, -0.5683])
    blue_R = np.max(reflectance_space[1:4], axis=0)
    blue_L = np.max(radiance_space[1:4], axis=0)
    green_R = reflectance_space[4]
    green_L = radiance_space[4]

    chl_R = OCx(blue_R, green_R, a)
    chl_L = OCx(blue_L, green_L, a)
    return chl_R, chl_L


def OC4E(wavelengths, Ed, Lw, R_rs):
    sensor = load_MERIS()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    a = np.array([0.3255, -2.7677, 2.4409, -1.1288, -0.4990])
    blue_R = np.max(reflectance_space[1:4], axis=0)
    blue_L = np.max(radiance_space[1:4], axis=0)
    green_R = reflectance_space[4]
    green_L = radiance_space[4]

    chl_R = OCx(blue_R, green_R, a)
    chl_L = OCx(blue_L, green_L, a)
    return chl_R, chl_L


def OC3M(wavelengths, Ed, Lw, R_rs):
    sensor = load_MODISA()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    a = np.array([0.2424, -2.7423, 1.8017, 0.0015, -1.2280])
    blue_R = np.max(reflectance_space[1:4:2], axis=0)
    blue_L = np.max(radiance_space[1:4:2], axis=0)
    green_R = reflectance_space[5]
    green_L = radiance_space[5]

    chl_R = OCx(blue_R, green_R, a)
    chl_L = OCx(blue_L, green_L, a)
    return chl_R, chl_L


def OC3V(wavelengths, Ed, Lw, R_rs):
    sensor = load_VIIRS()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    a = np.array([0.2228, -2.4683, 1.5867, -0.4275, -0.7768])
    blue_R = np.max(reflectance_space[1:3], axis=0)
    blue_L = np.max(radiance_space[1:3], axis=0)
    green_R = reflectance_space[3]
    green_L = radiance_space[3]

    chl_R = OCx(blue_R, green_R, a)
    chl_L = OCx(blue_L, green_L, a)
    return chl_R, chl_L


def OC3C(wavelengths, Ed, Lw, R_rs):
    sensor = load_CZCS()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    a = np.array([0.3330, -4.3770, 7.6267, -7.1457, 1.6673])
    blue_R = np.max(reflectance_space[:2], axis=0)
    blue_L = np.max(radiance_space[:2], axis=0)
    green_R = reflectance_space[2]
    green_L = radiance_space[2]

    chl_R = OCx(blue_R, green_R, a)
    chl_L = OCx(blue_L, green_L, a)
    return chl_R, chl_L


all_algorithms = [OC4, OC4E, OC3M, OC3V, OC3C]
all_algorithm_labels = ["SeaWiFS", "MERIS", "MODIS", "VIIRS", "CZCS"]
