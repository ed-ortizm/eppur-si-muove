#!/usr/bin/env python3
from eppur import *
## Loading the data into a Table
data_file = "../gaia_data.vot"
data = Table.read(data_file)

## ​ Image for the radial velocity of the 5000 brightest stars relative to the Sun
# Radial velocity: it is already computed relative to the sun
# https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html
# radial_velocity : Radial velocity (double, Velocity[km/s] )
# Spectroscopic radial velocity in the solar barycentric reference frame.
# Table().sort('key') method updates the table
data.sort('phot_g_mean_mag')
# Theere are missing values for the velocity, therefore I'll work with the 5K
# most luminous stars for which threre is a radial velocity available.
# https://docs.astropy.org/en/stable/table/masking.html
rad_vel = data['radial_velocity']
# I'm going to create a mask to exclude missing values. The .mask attribute
# outputs a True for missing values, to get a True for non missing values I turn
# to np.logicalnot()
mask = np.logical_not(rad_vel.mask)
data_masked = data[mask][:5000]
ra = data_masked['ra']
dec = data_masked['dec']
ra = [ra[i] for i in range(len(ra))]*u.degree
dec = [dec[i] for i in range(len(dec))]*u.degree
## Now let's plot the radial velocities it with a color bar.
coord = SkyCoord(ra,dec)
gal_coor = coord.galactic
# wrapping angle
gal_coor.data.lon.wrap_angle = 180. * u.degree
values = data_masked['radial_velocity']
colors = [(0, 0, 1), (1, 0, 0)]  # B -> R
n_bin=7
title = 'Radial velocity map'
label= 'Rad vel Km/s'
aitoff_color_bar(gal_coor, values, colors, n_bin, title, label, s=3., alpha=1.)
