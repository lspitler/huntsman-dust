"""Find objects using image segmentation"""

from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
from photutils import detect_sources


def find_objects(image, threshold, FWHM=2.0, npixels=6):
    sigma = FWHM * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(image, threshold, npixels=npixels, filter_kernel=kernel)
    return segm
