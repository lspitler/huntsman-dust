import os
import utils
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm


def detect_mask_sources(script_dir=None,
                        data_dir=None,
                        slice_in=None, slice_out=None,
                        sigma=None, iters=None, box_size=None,
                        filter_size=None, plt_grid=None,
                        FWHM=None, npixels=None,
                        obj_name=None,
                        Ra=None, Dec=None, obj_radius=None):
    """Detecting and masking sources.

    This program detects and masks sources on two levels.
        (1.) A 2D background is determined by creating a grid of desired dimensions,
             sigma clipping sources within each box and iteratively determining
             background levels. By interpolating between these grids, a 2D background
             array is created. Discrete sources are identified based on two criteria:
                (i)  Sources must be a fixed sigma above the background. By convetions,
                     a source is identified if it is 3.0 sigma above the threshold,
                     but any other value of sigma is also acceptable.
                (ii) A source must have a minimum number of interconnected pixels above
                     the threshold described in section (i) for it to be considered a
                     source.
        (2.) If there is a galaxy or large source present in the image, this is
             seperately masked. The center of the galaxy is either determined using
             SESAME, for which an internet connection is required, or by supplying
             the Ra, Dec coordinates of the center. The radius to be masked is supplied
             in arcmins. A circular mask is created.

                 Agrs:
                     script_dir(str, required):  Directory where scripts are saved
                     data_dir(srt, required):    Directory with fits Data
                     sliceing image: This runs analysis on a small subset of the image
                                slice_in(int, optional):    Start value to slice image
                                slice_out(int, optional):   End value to slice image
                     sigma(float, required):     Sigma threshold
                     iters(int, required):       Iteration to be used to find background
                     box_size(int, required):    Grid size in 2D background estimation
                     filter_size(int, required): Filter size in pixels
                     plt_grid(boolean, required):If true show grid used for 2D background
                     FWHM(int, required):        Full Width Half Maximum of 2D circular gaussian
                                                 kernel used to filter the image prior
                                                 to thresholding. Input is in terms of pixels.
                     npixels(int, required):     Number of pixels to classify a source
                     obj_name(str, required):    Galaxy name. Ex: 'NGC6822'
                     Ra(str, optional):          Ra of galaxy in degrees
                     Dec(str, optional):         Dec of galaxy in degrees
                     obj_radius(float, required):Radius of galaxy in arcmins
    """

    # Finds fits file, reads it to import image data, header, WCS
    script_dir = script_dir
    os.chdir(script_dir)
    image_path, file = utils.find_fits_files(data_dir)
    image, header, wcs = utils.image_load(image_path)

    # Slices image and corresponding WCS object
    image = image[slice_in:slice_out, slice_in:slice_out]
    wcs = wcs[slice_in:slice_out, slice_in:slice_out]

    # 2D background estimation
    bkg, bkgrms = utils.background_2D(image, sigma=sigma, iters=iters, box_size=box_size,
                                      filter_size=filter_size, plt_grid=plt_grid)

    # Detects a threshold of sigma above background
    threshold = bkg.background + (sigma * bkgrms)

    # Finds objects using image segmentation
    segm = utils.find_objects(image, threshold, FWHM=FWHM, npixels=npixels)

    # Mask galaxy at given Ra, Dec, within raduis in arcmins
    masked_img, mask = utils.mask_galaxy(image, wcs, name=obj_name,
                                         Ra=0., Dec=0., radius=obj_radius)

    # Dispays Log stretched image, segmented image, Galaxy masked image
    utils.plt_fits(image, wcs, figure=1, title="Logarithmic scaled FITs image",
                   cmap='Greys_r', norm=LogNorm())
    utils.plt_fits(masked_img, wcs, figure=2, title="Masked Galaxy",
                   cmap='Greys_r', norm=LogNorm())
    utils.plt_fits(segm, wcs, figure=3, title="Sources detected in the image",
                   cmap=segm.cmap(random_state=12345), norm=None)


# This is activated iff this module is run, not if it is imported.
if __name__ == '__main__':

    detect_mask_sources(script_dir='/Users/amanchokshi/Desktop/Huntsman/Scripts/huntsman_dust',
                        data_dir='/Users/amanchokshi/Desktop/Huntsman/Data',
                        slice_in=2000, slice_out=2500,
                        sigma=3.0, iters=10, box_size=20,
                        filter_size=5, plt_grid=False,
                        FWHM=2.0, npixels=6,
                        obj_name='NGC6822',
                        Ra=296.234, Dec=-14.797, obj_radius=10)

    plt.show()
