"""Function to load Fits images"""

from astropy.io import fits
from astropy.wcs import WCS

def image_load(image_path):
    hdulist = fits.open(image_path)
    image = hdulist[0].data
    header = hdulist[0].header
    wcs = WCS(header)
    return image, header, wcs
    