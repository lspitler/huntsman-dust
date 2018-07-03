# Huntsman Dust

This project is being created to study Galactic cirrus. Galactic cirrus are
clouds of Interstellar Matter (ISM) extending above and below the plane of
the Milky Way. This project is being developed to utilise data from Macquarie
Universities Huntsman Telescope, but can be adapted to other data.

## Detecting and Masking Sources

This program first aims to efficiently disentangle foreground ISM from
background discrete sources. This program detects and masks sources at two
levels.

1.  A 2D background is determined by creating a grid of desired
    dimensions, sigma clipping sources within each box and iteratively
    determining background levels. By interpolating between these
    grids, a 2D background array is created. Discrete sources are
    identified based on two criteria:
     * Sources must be a fixed sigma above the background.
       By convention, a source is identified if it is 3.0 sigma
       above the threshold, but any other value of sigma is also
       acceptable.
     * A source must have a minimum number of interconnected
       pixels above the threshold for it to be considered a source.

2. If there is a galaxy or large source present in the image, this is
   seperately masked. The center of the galaxy is either determined
   using SESAME, for which an internet connection is required, or by
   supplying the Ra, Dec coordinates of the center. The radius to be
   masked is supplied in arcmins. A circular mask is created.       



## Power Spectrum Analysis

Now that sources have been masked, we can begin a power spectrum analysis
of the dusty data.

**Under Development**


## Authors

* **Aman Chokshi**
* **Lee Spitler**
