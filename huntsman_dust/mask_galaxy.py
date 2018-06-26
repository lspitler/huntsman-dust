"""This function masks a specific location (Ra, Dec)
within a certain radius given is arcminutes"""

from astropy.coordinates import SkyCoord
from astropy.coordinates import SkyCoord
from astropy.wcs.utils import skycoord_to_pixel
from astropy.wcs.utils import proj_plane_pixel_scales

import numpy as np
from astropy import units as u


def mask_galaxy(image, wcs, name=None, Ra=None, Dec=None, radius=None):
    # radius must be given in arcminutes
    # Dimentions of the image
    dim = (image.shape)
    y, x = dim[0], dim[1]

    # Finds the center of an object by inputting its name into SESAME
    # This requires an active internet connection
    # a, b are the coordinates of the center given in pixels
    center = SkyCoord.from_name(name)
    if name is None:
        center = SkyCoord(Ra, Dec)

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

    # indexes each pixel and checks if its is >= radius from center
    Y, X = np.ogrid[:y, :x]
    dist_from_center = np.sqrt((X - a)**2 + (Y - b)**2)
    mask = dist_from_center <= rad_pix
    masked_img = image.copy()
    masked_img[(mask > 0)] = 0

    return masked_img, mask
