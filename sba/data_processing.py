"""
Module with functions for processing (ir)radiance and reflectance data, e.g.
to clean spectra
"""

import numpy as np
import operator as op


comparators = {">": op.gt, ">=": op.ge, "==": op.eq, "<": op.lt, "<=": op.le}


def get_keys_with_label(data, *labels, exclude="sd"):
    """
    Get the keys (column names) in an AstroPy table that contain a phrase
    `label`. Any number of labels can be given. The output will contain the
    list for each. If only one label is given, the output is a single list;
    otherwise, it is a list of lists.
    """
    keys = [[key for key in data.keys() if (label in key and exclude not in key)] for label in labels]
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


def clip_to_zero(data, threshold=-1e-4):
    Lw_keys, R_rs_keys = get_keys_with_label(data, "Lw", "R_rs")
    for Lw_k, R_rs_k in zip(Lw_keys, R_rs_keys):
        ind = np.where((data[R_rs_k] < 0) & (data[R_rs_k] > threshold))
        data[R_rs_k][ind] = 0
        data[Lw_k][ind] = 0


def remove_rows_based_on_threshold(data, quantity, comparator, threshold, quiet=False):
    operator = comparators[comparator]
    keys = get_keys_with_label(data, quantity)
    remove_indices = [i for i, row in enumerate(data) if any(operator(row[key], threshold) for key in keys)]
    data.remove_rows(remove_indices)
    if not quiet:
        print(f"Removed {len(remove_indices)} rows with {quantity} {comparator} {threshold}")


def remove_negative_R_rs(data, clip=-1e-4, **kwargs):
    clip_to_zero(data, threshold=clip)
    remove_rows_based_on_threshold(data, "R_rs", "<", 0, **kwargs)
