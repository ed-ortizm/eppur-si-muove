#!/usr/bin/env python3


import numpy as np
import sys
# Plotting
import matplotlib
import matplotlib.pyplot as plt
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
# masks for blue and red shifted stars (it is binary but well)
blue_shifted_mask = data_masked['radial_velocity'] < 0
blue_shifted = data_masked[blue_shifted_mask]
red_shifted_mask  = data_masked['radial_velocity'] > 0
red_shifted = data_masked[red_shifted_mask]
# https://matplotlib.org/gallery/subplots_axes_and_figures/geo_demo.html
plt.figure()
plt.subplot(111, projection="aitoff")
plt.title("Aitoff: radial velocity")
plt.grid(True)
# s:The marker size in points**2
# alpha: The alpha blending value, between 0 (transparent) and 1 (opaque).
plt.scatter(blue_shifted['Gal longitude'],blue_shifted['Gal latitude'] ,\
marker='*', s=1., alpha=1,c='blue')
plt.scatter(red_shifted['Gal longitude'],red_shifted['Gal latitude'] ,\
marker='*', s=1., alpha=1,c='red')
plt.show()

# The Sun has an apparent magnitude of −27
# https://matplotlib.org/3.1.0/tutorials/colors/colorbar_only.html
