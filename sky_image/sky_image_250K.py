#!/usr/bin/env python3
from eppur import *
## Loading the data into a Table
data_file = "../gaia_data.vot"
data = Table.read(data_file)

## Create an image of the sky with the 250000 brightest stars in ​ aitoff
## projection in ​ galactic coordinates
# http://learn.astropy.org/rst-tutorials/Coordinates-Transform.html

# Retrieving ra and dec (IRCS) as provided by Gaia
ra = data['ra']
dec = data['dec']
# Transforming them into a Quantity data type
ra = [ra[i] for i in range(len(ra))]*u.degree
dec = [dec[i] for i in range(len(dec))]*u.degree
# Converting them into galactic coordinates
coord = SkyCoord(ra,dec)
gal_coor = coord.galactic
# wrapping angle
gal_coor.data.lon.wrap_angle = 180. * u.degree
# I need to save this info, I'll need it when plotting the velocity map
data['l'] = gal_coor.l
data['b'] = gal_coor.b

aitoff(coordinates=gal_coor,title='250K brightest soucrces from Gaia',s=1,alpha=0.05)
