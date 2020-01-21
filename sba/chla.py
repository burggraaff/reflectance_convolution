import numpy as np

from .response_curves import load_SeaWiFS, load_MERIS, load_MODISA, load_VIIRS, load_CZCS, load_Sentinel2A, load_SPECTACLE


def KT16_algorithm(B4, B5, B6):
    return 169* (B5 - (B4 + B6)/2) + 19.5


def KT16(wavelengths, Ed, Lw, R_rs):
    sensor = load_Sentinel2A()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    chl_R = Ha17_algorithm(*reflectance_space[3:6])
    chl_L = Ha17_algorithm(*radiance_space[3:6])

    return chl_R, chl_L


def GM09_algorithm(B, G):
    return 0.8 * (B/G)**(-4.3)


def GM09(wavelengths, Ed, Lw, R_rs):
    sensor = load_SPECTACLE()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    chl_R = GM09_algorithm(reflectance_space[5], reflectance_space[4])
    chl_L = GM09_algorithm(radiance_space[5], radiance_space[4])

    return chl_R, chl_L


def Ha17_algorithm(B3, B4):
    return 0.80 * np.exp(0.35 * B3/B4)


def Ha17(wavelengths, Ed, Lw, R_rs):
    sensor = load_Sentinel2A()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    chl_R = Ha17_algorithm(reflectance_space[2], reflectance_space[3])
    chl_L = Ha17_algorithm(radiance_space[2], radiance_space[3])

    return chl_R, chl_L


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

    a = np.array([0.32814, -3.20725, 3.22969, -1.36769, -0.81739])
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

    a = np.array([0.42487, -3.20974, 2.89721, -0.75258, -0.98259])
    blue_R = np.max(reflectance_space[1:4], axis=0)
    blue_L = np.max(radiance_space[1:4], axis=0)
    green_R = reflectance_space[4]
    green_L = radiance_space[4]

    chl_R = OCx(blue_R, green_R, a)
    chl_L = OCx(blue_L, green_L, a)
    return chl_R, chl_L


def OC6M(wavelengths, Ed, Lw, R_rs):
    sensor = load_MODISA()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    a = np.array([1.22914, -4.99423, 5.64706, -3.53426, 0.69266])
    blue_R = np.max(reflectance_space[[0,1,3,4]], axis=0)
    blue_L = np.max(radiance_space[[0,1,3,4]], axis=0)
    green_R = np.mean(reflectance_space[[6,8]], axis=0)
    green_L = np.mean(radiance_space[[6,8]], axis=0)

    chl_R = OCx(blue_R, green_R, a)
    chl_L = OCx(blue_L, green_L, a)
    return chl_R, chl_L


def OC3M(wavelengths, Ed, Lw, R_rs):
    sensor = load_MODISA()
    reflectance_space = sensor.band_average(wavelengths, R_rs)
    radiance_space = sensor.band_average(wavelengths, Lw) / sensor.band_average(wavelengths, Ed)

    a = np.array([0.26294, -2.64669, 1.28364, 1.08209, -1.76828])
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

    a = np.array([0.23548, -2.63001, 1.65498, 0.16117, -1.37247])
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

    a = np.array([0.31841, -4.56386, 8.63979, -8.41411, 1.91532])
    blue_R = np.max(reflectance_space[:2], axis=0)
    blue_L = np.max(radiance_space[:2], axis=0)
    green_R = reflectance_space[2]
    green_L = radiance_space[2]

    chl_R = OCx(blue_R, green_R, a)
    chl_L = OCx(blue_L, green_L, a)
    return chl_R, chl_L


satellite_algorithms = [OC6M, OC3M, OC4, OC4E, OC3V, OC3C, Ha17]
satellite_algorithm_labels = ["MODIS OC6", "MODIS OC3", "SeaWiFS OC4", "MERIS OC4", "VIIRS OC3", "CZCS OC3", "S2A (Ha+17)"]
