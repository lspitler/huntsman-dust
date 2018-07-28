[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1323168.svg)](https://doi.org/10.5281/zenodo.1323168)

# Huntsman Dust

This project is being created to study Galactic cirrus. Galactic cirrus are
clouds of Interstellar Matter (ISM) extending above and below the plane of
the Milky Way. This project is being developed to utilise data from Macquarie
Universities Huntsman Telescope, but can be adapted to other data.

# Getting Started

To start using the Huntsman-Dust package, there are two easy steps.

1. **Setup** the huntsman-dust package
2. **Running Code** from your huntsman-dust package

## Setup
### Installing from Source

The project source is in a GitHub repository at https://github.com/lspitler/huntsman-dust. To install using git on the command line:  

```
$ cd ~/Build  
$ git clone https://github.com/lspitler/huntsman-dust.git  
$ cd huntsman-dust  
$ pip install -r requirements.txt  
$ python setup.py install   
```

## Running Code

* The functions in the Huntsman-Dust package are designed to be run from the terminal.  
* Create a Symbolic Link to the `~/Build` folder using the command  
`ln -s ~/Build`  
* Scripts can now be run as follows   
`python huntsman_dust/power_spectrum.py`  
* Use the help flag to view the arguments for each function   
 `python huntsman_dust/power_spectrum.py -h`


# Functions

This package aims to help you detect and mask sources in your image, and
perform a power spectrum analysis on the masked data.

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
   separately masked. The centre of the galaxy is either determined
   using SESAME, for which an internet connection is required, or by
   supplying the Ra, Dec coordinates of the centre. The radius to be
   masked is supplied in arc-minutes. A circular mask is created.       



## Power Spectrum Analysis

Now that sources have been masked, we can begin a power spectrum analysis
of the dusty data. A 2D FFT is performed of the data and the masked data.
The 2D psd obtained is now azimuthally averaged to create a 1D psd, which
represents the radial power spectrum of the masked data.

## Fake Images

A set of simulated is generated. The fake image has gaussian sources, whose
parameters can be adjusted based on the following arguments. The amplitude(flux)
of the sources are distributed according to a power-law, with exponent gamma. The
default value of gamma is -1.25, in accordance with the Schechter luminosity
function. The standard deviations of the sources follow a 1/r^2 distribution.

A background with spatial fluctuations at various scales is created separately.
This is achieved by first creating the desired 2D power spectrum, which is a
radial power law in Fourier space. According to literature, the index of ISM
power law distribution is -2.9. Taking the inverse FFT of this p_law array gives
a background with the desired levels of spatial fluctuations in real space.

This background is added to the set of galaxies, to create an accurate simulated
data set


## Authors

* **Aman Chokshi**
* **Lee Spitler**
