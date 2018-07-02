import os
import glob
import numpy as np

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import SigmaClip
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import proj_plane_pixel_scales
from astropy import units as u
from photutils import detect_sources
from photutils import Background2D, MedianBackground
from util_plot import util_plot


def find_fits_files(data_dir, file_name=None):
    """Creates image path in data directory.

        Args:
            data_dir(str, required):    Directory which contains FITs. images

        Returns:
            image_path(str):    Path to fits. image
            file(srt):    file name
    """
    cwd = os.getcwd()
    os.chdir(data_dir)
    for file in glob.glob("*"+file_name+"*"):
        print('File name: ' + file)
    os.chdir(cwd)
    image_path = os.path.join(data_dir, file)
    return image_path, file


def image_load(image_path):
    """Returns image, header and wcs objects.

        Args:
            image_path(str, required):    Image path to particular FITs. File

        Returns:
            image(array):    This is the image data
            header(table):    This is the header object
            wcs:    World Coordinte System object
    """
    hdulist = fits.open(image_path)
    image = hdulist[0].data
    header = hdulist[0].header
    wcs = WCS(header)
    return image, header, wcs


def background_2D(image, sigma=None, iters=None, box_size=None, filter_size=None,
                  plt_grid=None):
    """2D background estimation.

    This function creates a 2D background estimate by dividing the image into a grid,
    defined by box_size.

        Args:
            image(array, required):   This is the image data
            sigma(float, required):   Sigma level
            iters(int, required):   Number of iterations
            box_size(int, required):    Defines the box dimesions, in pixels
            filter_size(int, required):    Defines the filter reach in pixels
            plt_grid(boolean):    Overplot grid on image

        Returns:
            bkg(array):    2D background level
            bkgrms(array):    RMS background
    """
    sigma_clip = SigmaClip(sigma=sigma, iters=iters)
    mask = (image == 0)
    bkg_estimator = MedianBackground()
    bkg = Background2D(image, box_size=box_size, filter_size=filter_size,
                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                       mask=mask, edge_method=u'pad')
    print('Background Median: ' + str(bkg.background_median))
    print('Background RMS median: ' + str(bkg.background_rms_median))
    if plt_grid is True:
        plt.imshow(bkg.background, origin='lower', cmap='Greys')
        bkg.plot_meshes(outlines=True, color='#1f77b4')
    bkgrms = bkg.background_rms
    return bkg, bkgrms


def find_objects(image, threshold, FWHM=None, npixels=None):
    """Find sources in image by a segmentation process.

    This function detects sources a given sigma above a threshold,
    only if it has more that npixels that are interconnected.

        Args:
            image(array, required):        This is the image data
            threshold(array, required):    This is the threshold above which detection occurs
            FWHM(int, required):           Full Width Half Maximum of 2D circular gaussian
                                           kernel used to filter the image prior to
                                           thresholding. Input is in terms of pixels.
            npixels(int, required):        The minimum number of pixels to define a sources

        Returns:
            Segm:    The segmentation image
    """
    sigma = FWHM * gaussian_fwhm_to_sigma
    kernel = Gaussian2DKernel(sigma, x_size=3, y_size=3)
    kernel.normalize()
    segm = detect_sources(image, threshold, npixels=npixels, filter_kernel=kernel)
    return segm


def mask_galaxy(image, wcs, name=None, Ra=None, Dec=None, radius=None):
    """Masks galaxy at Ra, Dec within a radius given in arcminutes

    Creates a circular mask centered at a given Ra, Dec. The radius
    is given in arcmins. The wcs object is used to convert these inputs
    to pixel locations. A pixel scale is also determined. If the object
    name is suppled, SESAME is used to find object center. If no active
    internet connection is available, center location must be manually
    entered, in degrees. If no center coordinates are supplied, (0, 0)
    is the default center.

        Args:
            image(array, required):    Image data
            wcs:    World Coordinte System object
            name(str, optional):    Name of galaxy or object
            Ra(str):    Right Ascention
            Dec(str):    Declination
            Radius(float, required):    Radius to be masked, in arcminutes

        Returns:
            masked_img(array):    Image which has been masked
            mask(boolean array):    Mask of the given object"""
    # Radius must be given in arcminutes
    # Dimentions of the image
    dim = (image.shape)
    y, x = dim[0], dim[1]

    # Finds the center of an object by inputting its name into SESAME
    # This requires an active internet connection
    # a, b are the coordinates of the center given in pixels

    try:
        center = SkyCoord.from_name(name)
    except Exception:
        print("No active internet connection. Manually enter Ra, Dec.")
        center = SkyCoord(Ra, Dec, unit="deg")

    c_pix = skycoord_to_pixel(center, wcs)
    a, b = c_pix[0], c_pix[1]
    print(center)

    radius = radius*u.arcmin

    # Finds pixel scale using WSC object. The default units can be found by
    # unit = header['CUNIT1'], they are degrees by convention
    # degrees are converted to arcmins, after which radius in computed in pixels
    scale = proj_plane_pixel_scales(wcs)
    pix_scale = scale[0]*u.deg.to(u.arcmin)
    print('Image Scale: ' + str(pix_scale)+' arcmin/pix')

    rad_pix = (radius/pix_scale).value

    # Indexes each pixel and checks if its is >= radius from center
    Y, X = np.ogrid[:y, :x]
    dist_from_center = np.sqrt((X - a)**2 + (Y - b)**2)
    mask = dist_from_center <= rad_pix
    masked_img = image.copy()
    masked_img[(mask > 0)] = 0

    return masked_img, mask


def plt_fits(image, wcs, figure=None, title=None, cmap=None, norm=None):
    """Plots FITs images with axis given in Ra, Dec.

        Args:
            image(array):   Image data
            wcs:    World Coordinte System object
            figure(optional):    Figure Number
            title(str, optional):    Title of the figure
            cmap(str, optiona):    Color map
            norm:    Image normalizatuion

    """
    util_plot()
    fig = plt.figure(num=figure)
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    ax.imshow(image, origin='lower', cmap=cmap, norm=norm)
    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('DEC')
    ax.set_title(title)
