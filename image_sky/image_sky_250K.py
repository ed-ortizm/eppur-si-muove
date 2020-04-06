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

aitoff(coordinates=gal_coor,title='250K brightest soucrces from Gaia',s=0.05,alpha=0.05)

# ss = np.linspace(0.,0.5.,10)
# alphas = np.linspace(0.,1,10)
# for s in ss:
#     for alpha in alphas:
#         title = '250K brightest soucrces from Gaia' + '_s_' + str(s)[:4] + '_alpha_' + str(alpha)[:4]
#         aitoff(coordinates=gal_coor,title=title,s=s,alpha=alpha)
