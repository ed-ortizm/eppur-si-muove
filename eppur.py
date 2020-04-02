#!/usr/bin/env python3
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
from astropy.coordinates import SkyCoord

# Aitoff projection ​ in !!galactic coordinates!!.
data_file = "gaia_data.vot"
data = Table.read(data_file)
# Summary information about the table
data.info
# http://learn.astropy.org/rst-tutorials/Coordinates-Transform.html
# create a SkyCoord() object from the right ascension and the declination.
# coordinates: ICRS(Ra,Dec), the default frame
ra = data['ra']
dec = data['dec']
coord = SkyCoord(ra,dec)
gal_coor = coord.galactic
# wrapping angle
gal_coor.data.lon.wrap_angle = 180. * u.degree
# I need to save this info, I'll need it when plotting the velocity map
data['Gal longitude'] = gal_coor.l
data['Gal latitude'] = gal_coor.b
# https://matplotlib.org/gallery/subplots_axes_and_figures/geo_demo.html
#plt.figure()
#plt.subplot(111, projection="aitoff")
#plt.title("Aitoff")
#plt.grid(True)
# s:The marker size in points**2
# alpha: The alpha blending value, between 0 (transparent) and 1 (opaque).
#plt.scatter(gal_coor.data.lon.radian, gal_coor.data.lat.radian,marker='*', s=1., alpha=0.05)
#plt.show()
#plt.savefig("250K_brigtest_gaia_sources.pdf")

## Radial velocity: it is already computed respect to the sun
# https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html
# radial_velocity : Radial velocity (double, Velocity[km/s] )
# Spectroscopic radial velocity in the solar barycentric reference frame.
# sort method updates the table --> .sort??: Sort the table according to one or
# more keys. This operates on the existing table and does not return a new table.
data.sort('phot_g_mean_mag')# careful, check mag lower brighter??
# The data is sorted but there are missing values
rad_vel = data['radial_velocity']
# https://docs.astropy.org/en/stable/table/masking.html
# I'm going to create a mask to exclude missing values.
# The .mask attribute outputs a True for missing values, to get a True for
# non missing values I turn to np.logicalnot()
mask = np.logical_not(rad_vel.mask)
data_masked = data[mask][:5000]
# I'll need the min and max radial velocity to create the color bar
min_vel = np.min(data_masked['radial_velocity'])
max_vel = np.max(data_masked['radial_velocity'])

## Now let's plot the radial velocities it with a color bar.
# Actually I found something really cool
# https://matplotlib.org/3.1.1/gallery/color/custom_cmap.html
# https://stackoverflow.com/questions/25748183/python-making-color-bar-that-runs-from-red-to-blue

# Here I define my color map
colors = [(0, 0, 1), (1, 0, 0)]  # B -> R
n_bin = 7  # Discretizes the interpolation into bins
# Fewer bins will result in "coarser" colomap interpolation
cmap_name = 'my_list'
my_cm = mcol.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
normalize = mcol.Normalize(vmin=min_vel,vmax=max_vel)
mappable = cm.ScalarMappable(norm=normalize, cmap = my_cm)
mappable.set_array([])
# Now plotting
plt.figure()
plt.subplot(111, projection="aitoff")
plt.title("Aitoff")
plt.scatter(data_masked['Gal longitude'],data_masked['Gal latitude'] ,\
marker='*', s=6., alpha=1,c=mappable.to_rgba(data_masked['radial_velocity']))
plt.colorbar(mappable, label="radial velocity")
plt.show()
plt.savefig("radial_velocity_" + str(n_bin) + ".pdf")
