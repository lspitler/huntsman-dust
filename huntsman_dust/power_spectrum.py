import os
import matplotlib.pyplot as plt
import huntsman_dust.utils as utils
import huntsman_dust.util_plot as util_plot
from huntsman_dust.fake_image import fake_image
from huntsman_dust.detect_mask_sources import detect_mask_sources
import numpy as np

def power_spectrum(s_in,
                   s_out,
                   obj_name,
                   obj_rad,
                   Ra, Dec,
                   img_path=None,
                   plt_name=None):

    if img_path is not None:
        image, img_masked, _, _ = detect_mask_sources(img_path=img_path,
                                                      s_in=s_in,
                                                      s_out=s_out,
                                                      obj_name=obj_name,
                                                      Ra=Ra,
                                                      Dec=Dec,
                                                      obj_rad=obj_rad)

        psd1D_im = utils.p_spec(image)
        psd1D_msk = utils.p_spec(img_masked)

        util_plot.util_plot()
        plt.figure(1)

        plt.plot(psd1D_im, color='#ff8822', label='Unmasked Spectrum')
        plt.plot(psd1D_msk, color='#338899', label='Masked Spectrum')

        plt.title('Angular Power Spectrum of Input Image')
        plt.xlabel('Spatial Frequency')
        plt.ylabel('Power Spectrum')
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()

    else:
        sources, background, sim_sky = fake_image()
        psd_sou = utils.p_spec(sources)
        psd_bac = utils.p_spec(background)
        psd_sky = utils.p_spec(sim_sky)

        util_plot.util_plot()
        plt.figure(1)

        plt.plot(psd_sou, color='#ff8822', label='Sources')
        plt.plot(psd_bac, color='#338899', label='Background')
        plt.plot(psd_sky, color='#990022', label='Combined')

        plt.title('Angular Power Spectrum of Fake Image')
        plt.xlabel('Spatial Frequency')
        plt.ylabel('Power Spectrum')
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()

    if plt_name is not None:
        if img_path is not None:
            data_path, file = os.path.split(img_path)
        else:
            data_path = './'
        file_path = os.path.join(data_path, plt_name + "."+'png')
        plt.savefig(file_path)
    plt.show()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description="Creates radially average power spectrum")

    parser.add_argument('--img_path', default=None, metavar='\b',
                        help='Path to .fits file')
    parser.add_argument('--s_in', type=int, default=None, metavar='\b',
                        help='Slices image, in value')
    parser.add_argument('--s_out', type=int, default=None, metavar='\b',
                        help='Slices image, out value')
    parser.add_argument('--obj_name', default=None, metavar='\b',
                        help='Name of object to be masked')
    parser.add_argument('--obj_rad', type=float, default=10, metavar='\b',
                        help='Radius of object to be masked')
    parser.add_argument('--Ra', default=None, metavar='\b',
                        help='Inactive internet. Provide Ra of object center')
    parser.add_argument('--Dec', default=None, metavar='\b',
                        help='Inactive internet. Provide Dec of object center')
    parser.add_argument('--plt_name', '--output_plot_filename',
                        default=None, metavar='\b',
                        help='Name of plot file to be saved')

    args = parser.parse_args()

    power_spectrum(**vars(args))
