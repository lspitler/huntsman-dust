"""Import and display FITs image"""

import os
import glob
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.io import fits


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
header = hdulist[0].data

# print(hdu.data.shape)
# print(hdu.header)

plt.figure(1)
plt.imshow(image, origin='lower', cmap='gray')
plt.colorbar()

plt.figure(2)
plt.imshow(image, origin='lower', cmap='gray', norm=LogNorm())
plt.colorbar()
plt.show()
