"""
Module with functions for processing (ir)radiance and reflectance data, e.g.
to clean spectra
"""

import numpy as np
import operator as op
from astropy import table


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


def convert_to_unit_single(data, key, unit_old="", unit_new=None):
    data[key].unit = unit_old
    if unit_new is not None:
        data[key] = data[key].to(unit_new)


def convert_to_unit(data, label, unit_old="", unit_new=None):
    keys = get_keys_with_label(data, label)
    for key in keys:
        convert_to_unit_single(data, key, unit_old=unit_old, unit_new=unit_new)


def rename_columns(data, label_old, label_new, strip=False, **kwargs):
    keys_old = get_keys_with_label(data, label_old, **kwargs)
    for key in keys_old:
        key_new = key.replace(label_old, label_new)
        if strip:  # Strip units from key
            key_new = key_new.split(" ")[0]
        data.rename_column(key, key_new)


def add_Lw_from_Ed_Rrs(data, Ed_label="Ed", R_rs_label="R_rs"):
    wavelengths, Ed = split_spectrum(data, Ed_label)
    wavelengths, R_rs = split_spectrum(data, R_rs_label)

    Lw = Ed * R_rs

    Ed_keys, R_rs_keys = get_keys_with_label(data, Ed_label, R_rs_label)
    Lw_keys = [key.replace("Ed", "Lw") for key in Ed_keys]

    Lw_unit = data[Ed_keys[0]].unit * data[R_rs_keys[0]].unit

    Lw_table = table.Table(data=Lw, names=Lw_keys)
    convert_to_unit(Lw_table, "Lw", Lw_unit)

    data = table.hstack([data, Lw_table])
    return data

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
