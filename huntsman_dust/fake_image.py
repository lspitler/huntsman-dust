import numpy as np
import matplotlib.pyplot as plt
import huntsman_dust.utils as utils

from collections import OrderedDict
from scipy import fftpack
from photutils.datasets import make_noise_image
from photutils.datasets import make_random_models_table
from photutils.datasets import make_gaussian_sources_image


def fake_image(n_sources=100,
               shape=[512, 512],
               amplitude_r=[0, 20000],
               std_dev=[0, 7],
               random_state=666,
               noise={'type': None, 'mean': None, 'stddev': None}):
    """Creates fake image with Gaussian sources.

    Creates a fake image with gaussian sources, whose parameters can be
    adjusted based on the following arguments.

    A background with spatial fluctuations at various scales is
    created seperately. This is achived by first creating the desired 2D power
    spectrum, which is a radial power law in Fourier space. According to
    litrature, the index of ISM power law distribution is -2.9. Taking the
    inverse FFT of this p_law array gives a background with the desired levels
    of spatial fluctuations in real space.

        Args:
            n_sources(int):       Number of sources
            shape(2-tuple):       Dimensions of the image
            amplitude_r(list):    Range of amplitudes of sources
            std_dev(list):        Range of standard deviations of sources
            random_state(int):    Seed for random number generator
            noise(dictionary):    Parameters for noise
                (i)   type:       Gaussian or Poisson
                (ii)  mean:       Mean value of noise
                (iii) stddev:     Standard deviation of gaussian noise
    """

    param_ranges = [('amplitude', [amplitude_r[0], amplitude_r[1]]),
                    ('x_mean', [0, shape[1]]),
                    ('y_mean', [0, shape[0]]),
                    ('x_stddev', [std_dev[0], std_dev[1]]),
                    ('y_stddev', [std_dev[0], std_dev[1]]),
                    ('theta', [0, np.pi])]
    param_ranges = OrderedDict(param_ranges)
    sources = make_random_models_table(n_sources,
                                       param_ranges,
                                       random_state=random_state)

    if noise['type'] is None:
        sources = make_gaussian_sources_image(shape, sources)
    else:
        sources = sources + make_noise_image(shape,
                                             type=noise['type'],
                                             mean=noise['mean'],
                                             stddev=noise['stddev'])

    # CREATING BACKGROUNG (ISM)
    # The objective is to create a background with different levels of
    # spatial fluctuations built in. This is achived by first creating the
    # desired 2D power spectrum, which is a radial power law in Fourier space.
    # According to litrature, the index of ISM power law distribution is -2.9
    # Taking the inverse FFT of this p_law array gives a background with the
    # desired levels of spatial fluctuations in real space.

    p_law = np.zeros(shape, dtype=float)
    y, x = np.indices(p_law.shape)
    center = np.array([(y.max()-y.min())/2.0, (x.max()-x.min())/2.0])
    r = np.hypot(x - center[1], y - center[0])

    r_ind = r.astype(int)
    r_max = r.max().astype(int)

    a = np.arange(0.1, r_max+1.1, 1)      # These values control size of clouds
    b = 10**11*a**(-2.9)              # This controls magnitude of background

    for i in range(0, r_max+1):
        p_law[r_ind > i] = b[i]

    magnitude = np.sqrt(p_law)
    phase = 2*np.pi*np.random.randn(shape[0], shape[1])
    FFT = magnitude * np.exp(1j*phase)
    background = np.abs((fftpack.ifft2(fftpack.fftshift(FFT))))

    sim_sky = sources + background

    return sources, background, sim_sky


if __name__ == '__main__':
    sources, background, sim_sky = fake_image()

    utils.plt_image(sim_sky,
                    figure=1,
                    title="Simulated Sources With Background",
                    xlabel="Pixels",
                    ylabel="Pixels",
                    cmap='Greys_r')

    plt.show()
