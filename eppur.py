import numpy as np
import sys
# Plotting
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcol
import matplotlib.cm as cm
# Importing from astropy
# To load the data downloaded from Gaia
from astropy.table import Table
# To convert the units
from astropy import units as u
# SkyCoord --> convert coordinates between different references systems
from astropy.coordinates import SkyCoord, Distance
from astropy.coordinates import Angle
# movie
from astropy.time import Time
# Aitoff projection ​ in !!galactic coordinates!!.

def aitoff(coordinates,title,s,alpha):
    # https://matplotlib.org/gallery/subplots_axes_and_figures/geo_demo.html
    plt.figure()
    plt.subplot(111, projection="aitoff")
    plt.title(title)
    plt.grid(True)
    # s:The marker size in points**2
    # alpha: The alpha blending value, between 0 (transparent) and 1 (opaque).
    plt.scatter(coordinates.data.lon.radian, coordinates.data.lat.radian,marker='*', s=s, alpha=alpha)
    plt.savefig(title + ".png")
