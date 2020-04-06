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
    plt.close()

def aitoff_color_bar(coordinates,values,colors,n_bin,title,label,s,alpha, evol=False):
    # https://matplotlib.org/3.1.1/gallery/color/custom_cmap.html
    # https://stackoverflow.com/questions/25748183/python-making-color-bar-that-runs-from-red-to-blue
    # Defining the color map
    #n_bin = 7  # Discretizes the interpolation into bins
    # Fewer bins will result in "coarser" colomap interpolation
    cmap_name = 'my_list'
    #for n_bin in range(2,11):
    my_cm = mcol.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
    vmin,vmax = np.min(values),np.max(values)
    normalize = mcol.Normalize(vmin=vmin,vmax=vmax)
    mappable = cm.ScalarMappable(norm=normalize, cmap = my_cm)
    mappable.set_array([])
    # Now plotting
    plt.figure()
    plt.subplot(111, projection="aitoff")
    plt.title(title)
    plt.scatter(coordinates.data.lon.radian,coordinates.data.lat.radian ,\
    marker='*', s=s, alpha=alpha,c=mappable.to_rgba(values))
    plt.colorbar(mappable, label=label)
    if evol:
        plt.savefig( title + ".png")
    else:
        plt.savefig( title + "n_bin_" + str(n_bin)+ ".png")
    plt.close()
