#!/usr/bin/env python3
from eppur import *
## Loading the data into a Table
data_file = "gaia_data.vot"
data = Table.read(data_file)

## Create an image of the sky with the 250000 brightest stars in ​ aitoff
## projection in ​ galactic coordinates
# http://learn.astropy.org/rst-tutorials/Coordinates-Transform.html

# Retrieving ra and dec (IRCS) as provided by Gaia
# ra = data['ra']
# dec = data['dec']
# # Transforming them into a Quantity data type
# ra = [ra[i] for i in range(len(ra))]*u.degree
# dec = [dec[i] for i in range(len(dec))]*u.degree
# # Converting them into galactic coordinates
# coord = SkyCoord(ra,dec)
# gal_coor = coord.galactic
# # wrapping angle
# gal_coor.data.lon.wrap_angle = 180. * u.degree
# # I need to save this info, I'll need it when plotting the velocity map
# data['l'] = gal_coor.l
# data['b'] = gal_coor.b
#
# #aitoff(coordinates=gal_coor,title='250K brightest soucrces from Gaia',s=1,alpha=0.05)
#
# ## ​ Image for the radial velocity of the 5000 brightest stars relative to the Sun
# # Radial velocity: it is already computed relative to the sun
# # https://gea.esac.esa.int/archive/documentation/GDR2/Gaia_archive/chap_datamodel/sec_dm_main_tables/ssec_dm_gaia_source.html
# # radial_velocity : Radial velocity (double, Velocity[km/s] )
# # Spectroscopic radial velocity in the solar barycentric reference frame.
# # Table().sort('key') method updates the table
# data.sort('phot_g_mean_mag')
# # Theere are missing values for the velocity, therefore I'll work with the 5K
# # most luminous stars for which threre is a radial velocity available.
# # https://docs.astropy.org/en/stable/table/masking.html
# rad_vel = data['radial_velocity']
# # I'm going to create a mask to exclude missing values. The .mask attribute
# # outputs a True for missing values, to get a True for non missing values I turn
# # to np.logicalnot()
# mask = np.logical_not(rad_vel.mask)
# data_masked = data[mask][:5000]
# ra = data_masked['ra']
# dec = data_masked['dec']
# ra = [ra[i] for i in range(len(ra))]*u.degree
# dec = [dec[i] for i in range(len(dec))]*u.degree
# ## Now let's plot the radial velocities it with a color bar.
# coord = SkyCoord(ra,dec)
# gal_coor = coord.galactic
# # wrapping angle
# gal_coor.data.lon.wrap_angle = 180. * u.degree
# values = data_masked['radial_velocity']
# colors = [(0, 0, 1), (1, 0, 0)]  # B -> R
# n_bin=7
# title = 'Radial velocity map'
# label= 'Rad vel Km/s'
# aitoff_color_bar(gal_coor, values, colors, n_bin, title, label, s=3., alpha=1.)

## Images showing how the sky (1K brightest stars) will look around some of the
#famous stars
data.sort('phot_g_mean_mag')
data_1K = data[:1000]
ra = data_1K['ra']
dec = data_1K['dec']
ra = [ra[i] for i in range(len(ra))]*u.degree
dec = [dec[i] for i in range(len(dec))]*u.degree
coord = SkyCoord(ra,dec)
gal_coor = coord.galactic
l = gal_coor.l
b = gal_coor.b
l = [l[i] for i in range(len(l))]*u.degree
b = [b[i] for i in range(len(b))]*u.degree
# wrapping angle
gal_coor.data.lon.wrap_angle = 180. * u.degree
values = data_1K['phot_g_mean_mag']
colors = [(1, 0, 0), (0, 1, 0),(0, 0, 1)]
n_bin = 1000  # Discretizes the interpolation into bins
label= 'g magnitude'

# To do the offsets, I'll follow this:
# https://docs.astropy.org/en/stable/coordinates/matchsep.html
# For this, I'll chose the three brightest stars seen from earth
# http://www.astro.wisc.edu/~dolan/constellations/extra/brightest.html
# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Sirius
#(1) Sirius: Gal coord. (ep=J2000) :	227.23028548 -08.89028243

## Creating Sirius frame of reference
sirius_lon = Angle(227.23028548*u.degree).wrap_at(180*u.degree)
sirius_lon = sirius_lon.value *u.degree
sirius_lat = -08.89028243 *u.degree
center = SkyCoord(l = sirius_lon, b =sirius_lat, frame = 'galactic')
sirius_frame = center.skyoffset_frame()
target = SkyCoord(l, b, frame='galactic')
target = target.transform_to(sirius_frame)
title = '1K brigtest stars from Sirius'
aitoff_color_bar(target, values, colors, n_bin, title, label, s=10., alpha=1.)


# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=canopus&submit=submit+id
# (2) Canopus: Gal coord. (ep=J2000) :	261.21209728 -25.29220489
canopus_lon = Angle(261.21209728*u.degree).wrap_at(180*u.degree)
canopus_lon = sirius_lon.value *u.degree
canopus_lat = -25.29220489 *u.degree
center = SkyCoord(l= canopus_lon, b= canopus_lat,frame = 'galactic')
canopus_frame = center.skyoffset_frame()
target = SkyCoord(l,b,frame='galactic')
target = target.transform_to(canopus_frame)
title = '1K brigtest stars from Canopus'
aitoff_color_bar(target, values, colors, n_bin, title, label, s=10., alpha=1.)

# http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Rigil+Kentaurus&submit=submit+id
#(3) Rigil Kentaurus: Gal coord. (ep=J2000) :	315.73416514 -00.67965348
rigil_lon = Angle(315.73416514*u.degree).wrap_at(180*u.degree)
rigil_lon = sirius_lon.value *u.degree
rigil_lat = -0.67965348 *u.degree
center = SkyCoord(l= rigil_lon, b= rigil_lat, frame = 'galactic')
rigil_frame = center.skyoffset_frame()
target = SkyCoord(l,b,frame='galactic')
target = target.transform_to(rigil_frame)
title = '1K brigtest stars from Rigil Kentaurus'
aitoff_color_bar(target, values, colors, n_bin, title, label, s=10., alpha=1.)

### 10 second movie (at 60 frames per second!) showing how the apparent position
# of the 5000 brightest stars will change over the next 10000 years
# https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html
### Got this error:
# ValueError: Some parallaxes are negative, which are notinterpretable as
# distances. See the discussion in this paper: https://arxiv.org/abs/1507.02105 .
# If you want parallaxes to pass through, with negative parallaxes instead
# becoming NaN, use the `allow_negative=True` argument.
## then I'll mask the negative values
# data.sort('phot_g_mean_mag')
# parallax_mask = data['parallax'] > 0
# data_5K = data[parallax_mask][:5000]
# parallax = [data_5K['parallax'][i] for i in range(5000)]*u.mas
# distance = Distance(parallax=parallax)
# Gaia data is in equatorial system so we pass all to galactic system to use
# the tutorial
# ra = data_5K['ra']
# ra = [ra[i] for i in range(5000)]*u.degree
# dec = data_5K['dec']
# dec = [dec[i] for i in range(5000)]*u.degree
# pmra = data_5K['pmra']
# pmra = [pmra[i] for i in range(5000)]* u.mas/u.year
# pmdec = data_5K['pmdec']
# pmdec = [pmdec[i] for i in range(5000)]* u.mas/u.year
# coord = SkyCoord(ra=ra,dec=dec,pm_ra_cosdec=pmra, pm_dec=pmdec)
# gal_coor = coord.galactic
# wrapping angle
# gal_coor.data.lon.wrap_angle = 180. * u.degree
# I need to save this info, I'll need it when plotting the velocity map
# data_5K['l'] = gal_coor.l
# data_5K['b'] = gal_coor.b
# data_5K['pm_l_cosb'] = gal_coor.pm_l_cosb
# data_5K['pm_b'] = gal_coor.pm_b
# c = SkyCoord(data_5K['l'],data_5K['b'], distance = distance,\
# pm_l_cosb=data_5K['pm_l_cosb'],pm_b=data_5K['pm_b'],frame='galactic',\
# obstime=Time('2015-06-24 11:00:00'))
# #pmra : Proper motion in right ascension direction (double, Angular Velocity[mas/year])
# # pmra = Proper motion in right ascension μα*≡μαcosδ of the source in ICRS at the
# # reference epoch ref_epoch. This is the local tangent plane projection of the
# # proper motion vector in the direction of increasing right ascension.
# #pmdec : Proper motion in declination direction (double, Angular Velocity[mas/year] )
# #Proper motion in declination μδ of the source at the reference epoch ref_epoch.
# #This is the projection of the proper motion vector in the direction of increasing declination.
#
# ## I need 600 frames, 60 frames per second (10 seconds). Therefore a single
# # frame contains 100/6 years --> dt = 100/6 years
#
# #frames = np.linspace(0.,10_000.,600)
# for time in range(6):
#     dt = time*100./6.
#     c_new = c.apply_space_motion(dt= dt * u.year)
#     c_new.data.lon.wrap_angle = 180. * u.degree
#     b = [c_new.b.value[i] for i in range(5000)]
#     l = [c_new.l.value[i] for i in range(5000)]
#     plt.figure()
#     plt.subplot(111, projection="aitoff")
#     plt.title("5K motion 10K years")
#     plt.grid(True)
#     # s:The marker size in points**2
#     # alpha: The alpha blending value, between 0 (transparent) and 1 (opaque).
#     plt.scatter(l,b,marker='*', s=2., alpha=1)
#     plt.savefig(f"{time:    03}.png")
