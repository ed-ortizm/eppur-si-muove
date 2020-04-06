#!/usr/bin/env python3
from eppur import *
## Loading the data into a Table
data_file = "../gaia_data.vot"
data = Table.read(data_file)

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
