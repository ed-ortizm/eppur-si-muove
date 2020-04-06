#!/usr/bin/env python3
from eppur import *
## Loading the data into a Table
data_file = "../gaia_data.vot"
data = Table.read(data_file)


### 10 second movie (at 60 frames per second!) showing how the apparent position
# of the 5000 brightest stars will change over the next 10000 years
# https://docs.astropy.org/en/stable/coordinates/apply_space_motion.html
### Got this error:
# ValueError: Some parallaxes are negative, which are notinterpretable as
# distances. See the discussion in this paper: https://arxiv.org/abs/1507.02105 .
# If you want parallaxes to pass through, with negative parallaxes instead
# becoming NaN, use the `allow_negative=True` argument.
## then I'll mask the negative values
data.sort('phot_g_mean_mag')
parallax_mask = data['parallax'] > 0
data_5K = data[parallax_mask][:5000]
parallax = [data_5K['parallax'][i] for i in range(5000)]*u.mas
distance = Distance(parallax=parallax)
# Gaia data is in equatorial system so we pass all to galactic system to use
# the tutorial
ra = data_5K['ra']
ra = [ra[i] for i in range(5000)]*u.degree
dec = data_5K['dec']
dec = [dec[i] for i in range(5000)]*u.degree
pmra = data_5K['pmra']
pmra = [pmra[i] for i in range(5000)]* u.mas/u.year
pmdec = data_5K['pmdec']
pmdec = [pmdec[i] for i in range(5000)]* u.mas/u.year
coord = SkyCoord(ra=ra,dec=dec,pm_ra_cosdec=pmra, pm_dec=pmdec)
gal_coor = coord.galactic
# wrapping angle
gal_coor.data.lon.wrap_angle = 180. * u.degree
data_5K['l'] = gal_coor.l
data_5K['b'] = gal_coor.b
data_5K['pm_l_cosb'] = gal_coor.pm_l_cosb
data_5K['pm_b'] = gal_coor.pm_b
c = SkyCoord(data_5K['l'],data_5K['b'], distance = distance,\
pm_l_cosb=data_5K['pm_l_cosb'],pm_b=data_5K['pm_b'],frame='galactic',\
obstime=Time(2015.15, format= 'decimalyear'))
# pmra : Proper motion in right ascension direction (double, Angular Velocity[mas/year])
# # pmra = Proper motion in right ascension μα*≡μαcosδ of the source in ICRS at the
# # reference epoch ref_epoch. This is the local tangent plane projection of the
# # proper motion vector in the direction of increasing right ascension.
# #pmdec : Proper motion in declination direction (double, Angular Velocity[mas/year] )
# #Proper motion in declination μδ of the source at the reference epoch ref_epoch.
# #This is the projection of the proper motion vector in the direction of increasing declination.
#
## I need 600 frames, 60 frames per second (10 seconds). Therefore a single
# frame contains 100/6 years --> dt = 100/6 years
# for Plotting
colors = [(1, 0, 0), (0, 1, 0),(0, 0, 1)]
n_bin = 1000  # Discretizes the interpolation into bins
label= 'g magnitude'
values = data_5K['phot_g_mean_mag']
for time in range(601):
    dt = time*100./6.
    c_new = c.apply_space_motion(dt= dt * u.year)
    c_new.data.lon.wrap_angle = 180. * u.degree
    title= f"{time:03}"
    aitoff_color_bar(c_new,values,colors,n_bin,title,label,s=2,alpha=1, evol=True)
    # b = [c_new.b.value[i] for i in range(5000)]
    # l = [c_new.l.value[i] for i in range(5000)]
    # plt.figure()
    # plt.subplot(111, projection="aitoff")
    # plt.title("5K motion 10K years")
    # plt.grid(True)
    # plt.scatter(l,b,marker='*', s=2., alpha=1)
    # plt.savefig(f"{time:03}.png")
