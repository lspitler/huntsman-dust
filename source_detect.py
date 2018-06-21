"""Detecting and masking sources"""

import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from im_path import im_path
from image_load import image_load
from background_2D import background_2D
from find_objects import find_objects

# Finds fits file, reads it to import image data, header, WCS
source_dir = '/Users/amanchokshi/Desktop/Huntsman/huntsman_dust'
os.chdir(source_dir)
image_path = im_path('/Users/amanchokshi/Desktop/Huntsman/Data')
image, header, wcs = image_load(image_path)

# Selects sub region of image 
image = image[1600:2300, 2000:3000].astype(float)

# 2D background estimation
bkg, bkgrms = background_2D(image, sigma=3., iters=10, box_size=20, 
                            filter_size=5)
# detects a threshold of 3 sigma above background
threshold = bkg.background + (3.0 * bkgrms)

# finds objects using image segmentation
segm = find_objects(image, threshold, FWHM=2.0)


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7))
ax1.imshow(image, origin='lower', cmap='Greys_r', norm=LogNorm())
ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345))

plt.show()