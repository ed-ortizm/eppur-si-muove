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
from astropy.coordinates import Angle
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
#colors = [(0, 0, 1), (1, 0, 0)]  # B -> R
#n_bin = 7  # Discretizes the interpolation into bins
# Fewer bins will result in "coarser" colomap interpolation
#cmap_name = 'my_list'
#for n_bin in range(2,11):
#    my_cm = mcol.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
#    normalize = mcol.Normalize(vmin=min_vel,vmax=max_vel)
#    mappable = cm.ScalarMappable(norm=normalize, cmap = my_cm)
#    mappable.set_array([])
#    # Now plotting
#    plt.figure()
#    plt.subplot(111, projection="aitoff")
#    plt.title("Aitoff")
#    plt.scatter(data_masked['Gal longitude'],data_masked['Gal latitude'] ,\
#    marker='*', s=6., alpha=1,c=mappable.to_rgba(data_masked['radial_velocity']))
#    plt.colorbar(mappable, label="radial velocity")
#    plt.savefig("radial_velocity_" + str(n_bin) + ".png")
    #plt.show()

## Images showing how the sky will look around some of the famous stars
# creating the targets
data_1K = data[:1000]
min_mag = np.min(data_1K['phot_g_mean_mag'])
max_mag = np.max(data_1K['phot_g_mean_mag'])

# Here I define my color map
colors = [(1, 0, 0), (0, 1, 0),(0, 0, 1)]
n_bin = 1000  # Discretizes the interpolation into bins
# Fewer bins will result in "coarser" colomap interpolation
cmap_name = 'my_list'
my_cm = mcol.LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bin)
normalize = mcol.Normalize(vmin=min_mag,vmax=max_mag)
mappable = cm.ScalarMappable(norm=normalize, cmap = my_cm)
mappable.set_array([])

# to do the offsets, I'll follow this:
# https://docs.astropy.org/en/stable/coordinates/matchsep.html
# For this, I'll chose the three brightest stars seen from earth
# http://www.astro.wisc.edu/~dolan/constellations/extra/brightest.html
# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sirius
#(1) Sirius: Gal coord. (ep=J2000) :	227.23028548 -08.89028243
## Creating Sirius frame of reference
sirius_lon = Angle(227.23028548*u.degree).wrap_at(180*u.degree)
sirius_lon = sirius_lon.value *u.degree
sirius_lat = -08.89028243 *u.degree
center = SkyCoord(sirius_lon, sirius_lat)
sirius_frame = center.skyoffset_frame()
target = SkyCoord(data_1K['Gal longitude'],data_1K['Gal latitude'],frame='galactic')
target = target.transform_to(sirius_frame)
data_1K['Sirius l'] = target.lon
data_1K['Sirius b'] = target.lat
plt.figure()
plt.subplot(111, projection="aitoff")
plt.title("Aitoff")
plt.grid(True)
# s:The marker size in points**2
# alpha: The alpha blending value, between 0 (transparent) and 1 (opaque).
plt.scatter(data_1K['Sirius l'], data_1K['Sirius b'],marker='*',\
c=mappable.to_rgba(data_1K['phot_g_mean_mag']), s=10., alpha=1)
plt.colorbar(mappable, label="magnitude")
plt.savefig("1K_brigtest_sources_Sirius.png")

# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=canopus&submit=submit+id
# (2) Canopus: Gal coord. (ep=J2000) :	261.21209728 -25.29220489
canopus_lon = Angle(261.21209728*u.degree).wrap_at(180*u.degree)
canopus_lon = sirius_lon.value *u.degree
canopus_lat = -25.29220489 *u.degree
center = SkyCoord(canopus_lon, canopus_lat)
canopus_frame = center.skyoffset_frame()
target = SkyCoord(data_1K['Gal longitude'],data_1K['Gal latitude'],frame='galactic')
target = target.transform_to(canopus_frame)
data_1K['Canopus l'] = target.lon
data_1K['Canopus b'] = target.lat
plt.figure()
plt.subplot(111, projection="aitoff")
plt.title("Aitoff")
plt.grid(True)
# s:The marker size in points**2
# alpha: The alpha blending value, between 0 (transparent) and 1 (opaque).
plt.scatter(data_1K['Canopus l'], data_1K['Canopus b'],marker='*',\
c=mappable.to_rgba(data_1K['phot_g_mean_mag']) ,s=10., alpha=1)
plt.colorbar(mappable, label="magnitude")
plt.savefig("1K_brigtest_sources_Canopus.png")

# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Rigil+Kentaurus&submit=submit+id
#(3) Rigil Kentaurus: Gal coord. (ep=J2000) :	315.73416514 -00.67965348
rigil_lon = Angle(315.73416514*u.degree).wrap_at(180*u.degree)
rigil_lon = sirius_lon.value *u.degree
rigil_lat = -0.67965348 *u.degree
center = SkyCoord(rigil_lon, rigil_lat)
rigil_frame = center.skyoffset_frame()
target = SkyCoord(data_1K['Gal longitude'],data_1K['Gal latitude'],frame='galactic')
target = target.transform_to(rigil_frame)
data_1K['rigil l'] = target.lon
data_1K['rigil b'] = target.lat
plt.figure()
plt.subplot(111, projection="aitoff")
plt.title("Aitoff")
plt.grid(True)
# s:The marker size in points**2
# alpha: The alpha blending value, between 0 (transparent) and 1 (opaque).
plt.scatter(data_1K['rigil l'], data_1K['rigil b'],marker='*',\
c=mappable.to_rgba(data_1K['phot_g_mean_mag']), s=10., alpha=1)
plt.colorbar(mappable, label="magnitude")
plt.savefig("1K_brigtest_sources_Rigil_Kentaurus.png")
# phot_g_mean_mag : G-band mean magnitude (float, Magnitude[mag])
# Mean magnitude in the G band. This is computed from the G-band mean
# flux applying the magnitude zero-point in the Vega scale.
