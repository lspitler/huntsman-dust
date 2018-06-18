"""Import and display FITs image"""

import os
import glob
import matplotlib.pyplot as plt
from astropy.io import fits

cwd = os.getcwd()
data_dir = '/Users/amanchokshi/Desktop/Huntsman/Data'
os.chdir(data_dir)
for file in glob.glob("*.fits"):
    print(file)
os.chdir(cwd)
    
image = os.path.join(data_dir, file)
hdulist = fits.open(image)
hdu = hdulist[0]

print(hdu.data.shape)

plt.imshow(hdu.data, origin='lower')
plt.show()