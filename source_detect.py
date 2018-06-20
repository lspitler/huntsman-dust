"""Detecting and masking sources"""

import os
import glob

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import Gaussian2DKernel
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from astropy.stats import SigmaClip

from photutils import detect_threshold
from photutils import detect_sources
from photutils import Background2D, MedianBackground
from photutils import deblend_sources


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

# Selects sub image
image = image[1600:2300, 2000:3000].astype(float)

# reads header and creates World Coordinate System (WCS) object
w = WCS(header)

# 2D background estimation
sigma_clip = SigmaClip(sigma=3., iters=10)
mask = (image == 0)
bkg_estimator = MedianBackground()
bkg = Background2D(image, (20, 20), filter_size=(10, 10),
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                   mask=mask, edge_method=u'pad')
print(bkg.background_median)
print(bkg.background_rms_median)
plt.imshow(bkg.background, origin='lower', cmap='Greys')
bkg.plot_meshes(outlines=True, color='#1f77b4')
bkgrms = bkg.background_rms

# detects a threshold of 3 sigma above background
# threshold = detect_threshold(image, snr=3.)
threshold = bkg.background + (3.0 * bkgrms)

sigma = 2.0 * gaussian_fwhm_to_sigma    # FWHM = 2.
kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
kernel.normalize()
segm = detect_sources(image, threshold, npixels=6, filter_kernel=kernel)


norm = ImageNormalize(stretch=SqrtStretch())
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 7))
ax1.imshow(image, origin='lower', cmap='Greys_r', norm=LogNorm())
ax2.imshow(segm, origin='lower', cmap=segm.cmap(random_state=12345))


plt.show()
