"""Import and display FITs image"""

import os
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


from im_path import im_path
from image_load import image_load
from plt_fits import plt_fits

source_dir = '/Users/amanchokshi/Desktop/Huntsman/huntsman_dust'
os.chdir(source_dir)
image_path, file = im_path('/Users/amanchokshi/Desktop/Huntsman/Data')
image, header, wcs = image_load(image_path)

image = image[1600:2300, 2000:3000]
wcs = wcs[1600:2300, 2000:3000]

plt_fits(image, wcs, title="Sample Title", cmap='Greys_r', norm=LogNorm())
