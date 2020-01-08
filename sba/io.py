from astropy.io.ascii import read
from numpy import loadtxt, genfromtxt

def write_data(data, label, **kwargs):
    label_lowercase = label.lower()
    data.write(f"data/{label_lowercase}_processed.tab", format="ascii.fast_tab", overwrite=True, **kwargs)
