import numpy as np
import matplotlib.pyplot as plt
import huntsman_dust.utils as utils

from matplotlib.colors import LogNorm


def detect_mask_sources(img_path,
                        obj_rad,
                        obj_name,
                        s_in,
                        s_out,
                        mask_fits_name=None,
                        img_mask_fits_name=None,
                        sigma=3.0,
                        iters=10,
                        box_size=20,
                        filter_s=5,
                        plt_grid=False,
                        FWHM=2.0,
                        npixels=6,
                        ds9_reg=False,
                        Ra=None,
                        Dec=None):
    """Detecting and masking sources.

    This program detects and masks sources on two levels.
        (1.) A 2D background is determined by creating a grid of desired
             dimensions, sigma clipping sources within each box and iteratively
             determining background levels. By interpolating between these
             grids, a 2D background array is created. Discrete sources are
             identified based on two criteria:
                (i)  Sources must be a fixed sigma above the background.
                     By convention, a source is identified if it is 3.0 sigma
                     above the threshold, but any other value of sigma is also
                     acceptable.
                (ii) A source must have a minimum number of interconnected
                     pixels above the threshold described in section (i) for
                     it to be considered a source.
        (2.) If there is a galaxy or large source present in the image, this is
             seperately masked. The center of the galaxy is either determined
             using SESAME, for which an internet connection is required, or by
             supplying the Ra, Dec coordinates of the center. The radius to be
             masked is supplied in arcmins. A circular mask is created.

                 Agrs:
                     image_path(srt, required):   Directory with fits Data
                     s_in(int, optional):         Slices image
                     s_out(int, optional):        Slices image
                     sigma(float, required):      Sigma threshold
                     iters(int, required):        Iterations to find background
                     box_size(int, required):     Grid size for 2D background
                     filter_size(int, required):  Filter size in pixels
                     plt_grid(boolean, required): Show 2D background grid
                     FWHM(int, required):         Full Width Half Maximum of 2D
                                                  circular gaussian kernel used
                                                  to filter the image prior to
                                                  thresholding. Input is in
                                                  terms of pixels.
                     npixels(int, required):      Min pixels to classify source
                     ds9_region(boolean, opt):    Creates ds9 region file
                     obj_name(str, required):     Galaxy name. Ex: 'NGC6822'
                     Ra(str, optional):           Ra of galaxy in degrees
                     Dec(str, optional):          Dec of galaxy in degrees
                     obj_radius(float, required): Radius of galaxy in arcmins
                     mask_fits_name(str, opt):    Name of output mask .fits
                     img_mask_fits_name(str, opt):Name of output msk img .fits
    """

    # Finds fits file, reads it to import image data, header, WCS
    image, header, wcs = utils.image_load(img_path)

    # Slices image and corresponding WCS object
    image = image[s_in:s_out, s_in:s_out]
    wcs = wcs[s_in:s_out, s_in:s_out]

    # 2D background estimation
    bkg, bkgrms = utils.background_2D(image,
                                      sigma=sigma,
                                      iters=iters,
                                      box_size=box_size,
                                      filter_size=filter_s,
                                      plt_grid=plt_grid)

    # Detects a threshold of sigma above background
    threshold = bkg.background + (sigma * bkgrms)

    # Finds objects using image segmentation
    segm = utils.find_objects(image,
                              threshold,
                              FWHM=FWHM,
                              npixels=npixels)

    # Creates ds9 region file if ds9_region is True
    utils.ds9_region(img_path,
                     image,
                     segm,
                     wcs,
                     ds9_region=ds9_reg)

    seg_mask = np.ones((image.shape[0], image.shape[1]), dtype=bool)
    seg_mask[segm.data_masked < 1] = 0

# Mask galaxy at given Ra, Dec, within raduis in arcmins
    if obj_name is not None:
        mask = utils.mask_galaxy(image,
                                 wcs,
                                 name=obj_name,
                                 radius=obj_rad,
                                 Ra=Ra,
                                 Dec=Dec)

        mask = ~np.ma.mask_or(mask, seg_mask)
    else:
        mask = seg_mask

    img_masked = image.copy()
    img_masked[mask < 1] = 0

# writes mask to binary fits file
    if mask_fits_name is not None:
        utils.fits_write(mask, header, img_path, name=mask_fits_name)

# writes masked image to binary fits file
    if img_mask_fits_name is not None:
        utils.fits_write(img_masked, header, img_path, name=img_mask_fits_name)

    return image, img_masked, wcs, segm


# This is activated iff this module is run, not if it is imported.
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description="Detects & masks source in 2D image")

    parser.add_argument('--img_path', required=True, metavar='\b',
                        help='Path to .fits file')
    parser.add_argument('--s_in', type=int, default=None, metavar='\b',
                        help='Slices image, in value')
    parser.add_argument('--s_out', type=int, default=None, metavar='\b',
                        help='Slices image, out value')
    parser.add_argument('--sigma', default=3.0, metavar='\b',
                        help='Sigma value to determine threshold')
    parser.add_argument('--iters', default=10, metavar='\b',
                        help='Number of iterations used to determine 2D' +
                        'background')
    parser.add_argument('--box_size', default=20, metavar='\b',
                        help='Grid size for 2D background')
    parser.add_argument('--filter_s', default=5, metavar='\b',
                        help='size of filter used on background')
    parser.add_argument('--plt_grid', default=False, metavar='\b',
                        help='Plots the 2D background grid')
    parser.add_argument('--FWHM', default=2.0, metavar='\b',
                        help='FWHM of 2D gaussian kernel used for filtering')
    parser.add_argument('--npixels', default=6, metavar='\b',
                        help='Min no. pix used to classify a source')
    parser.add_argument('--ds9_reg', default=False, metavar='\b',
                        help='Export ds9 regions file')
    parser.add_argument('--obj_name', default=None, metavar='\b',
                        help='Name of object to be masked')
    parser.add_argument('--obj_rad', type=float, default=10, metavar='\b',
                        help='Radius of object to be masked')
    parser.add_argument('--Ra', default=None, metavar='\b',
                        help='Inactive internet. Provide Ra of object center')
    parser.add_argument('--Dec', default=None, metavar='\b',
                        help='Inactive internet. Provide Dec of object center')
    parser.add_argument('-mask', '--mask_fits_name', type=str, default=None,
                        metavar='\b', help='Writes mask to .fits file with'
                        + 'name=mask_fits_name')
    parser.add_argument('-img_mask', '--img_mask_fits_name', type=str,
                        default=None, help='Writes img_masked to .fits file' +
                        'with name = img_mask_fits_name', metavar='\b')

    args = parser.parse_args()

    image, img_masked, wcs, segm = detect_mask_sources(**vars(args))

    utils.plt_fits(img_masked,
                   wcs,
                   figure=1,
                   title="Logarithmic scaled masked FITs image",
                   cmap='Greys_r',
                   norm=LogNorm())

    utils.plt_fits(image,
                   wcs,
                   figure=2,
                   title="Logarithmic scaled FITs image",
                   cmap='Greys_r',
                   norm=LogNorm())

    utils.plt_fits(segm,
                   wcs,
                   figure=3,
                   title="Sources detected in the image",
                   cmap=segm.cmap(random_state=12345),
                   norm=None)
    plt.show()
# '/Users/amanchokshi/Desktop/Huntsman/Data/ngc6822_r.fits'
