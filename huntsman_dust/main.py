"""Detecting and masking sources"""


def main():
    import os
    from matplotlib.colors import LogNorm

    # imports huntsman_dust functions
    from im_path import im_path
    from image_load import image_load
    from background_2D import background_2D
    from find_objects import find_objects
    from plt_fits import plt_fits

    # Finds fits file, reads it to import image data, header, WCS
    source_dir = '/Users/amanchokshi/Desktop/Huntsman/Scripts/huntsman_dust'
    os.chdir(source_dir)
    image_path, file = im_path('/Users/amanchokshi/Desktop/Huntsman/Data')
    image, header, wcs = image_load(image_path)

    # Selects sub region of image
    image = image[1600:2300, 2000:3000]
    wcs = wcs[1600:2300, 2000:3000]

    # 2D background estimation
    bkg, bkgrms = background_2D(image, sigma=3., iters=10, box_size=20,
                                filter_size=5, plt_grid=False)

    # detects a threshold of 3 sigma above background
    threshold = bkg.background + (3.0 * bkgrms)

    # finds objects using image segmentation
    segm = find_objects(image, threshold, FWHM=2.0, npixels=6)

    # dispays Log stretched image & segmented image
    plt_fits(image, wcs, figure=1, title="Logarithmic scaled FITs image",
             cmap='Greys_r', norm=LogNorm())
    plt_fits(segm, wcs, figure=2, title="Sources detected in the image",
             cmap=segm.cmap(random_state=12345), norm=None)


if __name__ == '__main__':
    main()
    import matplotlib.pyplot as plt
    plt.show()
