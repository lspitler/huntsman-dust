import matplotlib.pyplot as plt
from latexify import latexify


def plt_fits(image, wcs, figure=None, title=None, cmap=None, norm=None):
    latexify()
    fig = plt.figure(num=figure)
    ax = fig.add_subplot(1, 1, 1, projection=wcs)
    ax.imshow(image, origin='lower', cmap=cmap, norm=norm)
    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('DEC')
    ax.set_title(title)
