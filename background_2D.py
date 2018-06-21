"""2D background estimation"""

import matplotlib.pyplot as plt
from astropy.stats import SigmaClip
from photutils import Background2D, MedianBackground

def background_2D(image, sigma=3, iters=10, boxsize=20, filter_size=10):
    sigma_clip = SigmaClip(sigma=3., iters=10)
    mask = (image == 0)
    bkg_estimator = MedianBackground()
    bkg = Background2D(image,box_size=20, filter_size=10,
                   sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
                   mask=mask, edge_method=u'pad')
    print('Background Median: ' + str(bkg.background_median))
    print('Background RMS median: ' + str(bkg.background_rms_median))
    plt.imshow(bkg.background, origin='lower', cmap='Greys')
    bkg.plot_meshes(outlines=True, color='#1f77b4')
    bkgrms = bkg.background_rms
    return bkg, bkgrms