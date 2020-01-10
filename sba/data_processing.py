import numpy as np


def get_keys_with_label(data, *labels):
    """
    Get the keys (column names) in an AstroPy table that contain a phrase
    `label`. Any number of labels can be given. The output will contain the
    list for each. If only one label is given, the output is a single list;
    otherwise, it is a list of lists.
    """
    keys = [[key for key in data.keys() if label in key] for label in labels]
    if len(labels) == 1:  # If only one label is given, return one list
        return keys[0]
    else:  #
        return keys


def split_spectrum(data_table, label):
    keys_relevant = get_keys_with_label(data_table, label)
    wavelengths = np.array([float(key.split("_")[-1]) for key in keys_relevant])
    try:
        spectra = np.array([data_table[key]._data for key in keys_relevant]).T
    except AttributeError:
        spectra = np.array([data_table[key].data for key in keys_relevant]).T
    return wavelengths, spectra
