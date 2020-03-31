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
# create a SkyCoord() object from the right ascension and the declination.
# coordinates: ICRS(RA,Dec), the default frame
ra = data['ra']
dec = data['dec']
coordinates = SkyCoord(ra,dec)



# https://matplotlib.org/gallery/subplots_axes_and_figures/geo_demo.html
#plt.figure()
#plt.subplot(111, projection="aitoff")
#plt.title("Aitoff")
#plt.grid(True)
