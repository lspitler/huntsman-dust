"""Import and display FITs image"""

import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from im_path import im_path
from image_load import image_load


source_dir = '/Users/amanchokshi/Desktop/Huntsman/huntsman_dust'
os.chdir(source_dir)
image_path = im_path('/Users/amanchokshi/Desktop/Huntsman/Data')
image, header, wcs = image_load(image_path)


# Figure 1 - Linear
fig_1 = plt.figure(1)
fig_1.add_subplot(111, projection=wcs)
plt.imshow(image, origin='lower', cmap='gray')
plt.xlabel('RA')
plt.ylabel('Dec')
plt.colorbar()

# Figure 2 - Logarithmic
fig_2 = plt.figure(2)
fig_2.add_subplot(111, projection=wcs)
plt.imshow(image, origin='lower', cmap='gray', norm=LogNorm())
plt.xlabel('RA')
plt.ylabel('Dec')
plt.colorbar()
plt.show()
