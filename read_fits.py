"""Import and display FITs image"""

import os
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits
from astropy.wcs import WCS


# Finds FITs files and created image_path
cwd = os.getcwd()
data_dir = '/Users/amanchokshi/Desktop/Huntsman/Data'
os.chdir(data_dir)
for file in glob.glob("*.fits"):
    print(file)
os.chdir(cwd)
image_path = os.path.join(data_dir, file)

# Read in FITs file
hdulist = fits.open(image_path)

# Extract image and header
image = hdulist[0].data
header = hdulist[0].header

# reads header and creates World Coordinate System (WCS) object
w = WCS(header)

# Figure 1 - Linear
fig_1 = plt.figure(1)
fig_1.add_subplot(111, projection=w)
plt.imshow(image, origin='lower', cmap='gray')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.colorbar()

# Figure 2 - Logarithmic
fig_2 = plt.figure(2)
fig_2.add_subplot(111, projection=w)
plt.imshow(image, origin='lower', cmap='gray', norm=LogNorm())
plt.xlabel('RA')
plt.ylabel('Dec')
plt.colorbar()
plt.show()
