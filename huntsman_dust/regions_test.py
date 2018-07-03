import utils
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from photutils import source_properties, EllipticalAperture


def detect_mask_sources(image_path=None,
                        slice_in=None,
                        slice_out=None,
                        sigma=None,
                        iters=None,
                        box_size=None,
                        filter_size=None,
                        plt_grid=None,
                        FWHM=None,
                        npixels=None):
    """Tests whether all sources are detected.

        Agrs:
            data_dir(srt, required):    Directory with fits Data
            sliceing image: Runs analysis on a small subset of image
                    slice_in(int, optional):    Start value
                    slice_out(int, optional):   End value
            sigma(float, required):     Sigma threshold
            iters(int, required):       Iterations to find background
            box_size(int, required):    Grid size for 2D background
            filter_size(int, required): Filter size in pixels
            plt_grid(boolean, required):Show 2D background grid
            FWHM(int, required):        Full Width Half Maximum of 2D
                                        circular gaussian kernel used
                                        to filter the image prior to
                                        thresholding. Input is in
                                        terms of pixels.
            npixels(int, required):     Min pixels to classify source
    """

    # Finds fits file, reads it to import image data, header, WCS
    image, header, wcs = utils.image_load(image_path)

    # Slices image and corresponding WCS object
    image = image[slice_in:slice_out, slice_in:slice_out]
    wcs = wcs[slice_in:slice_out, slice_in:slice_out]

    # 2D background estimation
    bkg, bkgrms = utils.background_2D(image,
                                      sigma=sigma,
                                      iters=iters,
                                      box_size=box_size,
                                      filter_size=filter_size,
                                      plt_grid=plt_grid)

    # Detects a threshold of sigma above background
    threshold = bkg.background + (sigma * bkgrms)

    # Finds objects using image segmentation
    segm = utils.find_objects(image,
                              threshold,
                              FWHM=FWHM,
                              npixels=npixels)

    cat = source_properties(image, segm)
    r = 3
    apertures = []
    for obj in cat:
        position = (obj.xcentroid.value, obj.ycentroid.value)
        a = obj.semimajor_axis_sigma.value * r
        b = obj.semiminor_axis_sigma.value * r
        theta = obj.orientation.value
        apertures.append(EllipticalAperture(position, a, b, theta=theta))

    utils.plt_fits(image,
                   wcs,
                   figure=1,
                   title="Image with sources",
                   cmap='Greys',
                   norm=LogNorm())
    for aperture in apertures:
        aperture.plot(color='#ff7733', lw=2.0)

    utils.plt_fits(segm,
                   wcs,
                   figure=3,
                   title="Segmented image with sources",
                   cmap=segm.cmap(random_state=12345),
                   norm=None)
    for aperture in apertures:
        aperture.plot(color='white', lw=1.5)


# This is activated iff this module is run, not if it is imported.
if __name__ == '__main__':

    detect_mask_sources(image_path='/Users/amanchokshi/Desktop/Huntsman/Data'
                        + '/ngc6822_r.fits',
                        slice_in=1500,
                        slice_out=1600,
                        sigma=3.0,
                        iters=10,
                        box_size=20,
                        filter_size=5,
                        plt_grid=False,
                        FWHM=2.0,
                        npixels=6)
    plt.show()
