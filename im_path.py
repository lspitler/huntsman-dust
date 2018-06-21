"""Function to find image path"""

import os
import glob

def im_path(data_dir):
    cwd = os.getcwd()
    os.chdir(data_dir)
    for file in glob.glob("*.fits"):
        print(file)
    os.chdir(cwd)
    image_path = os.path.join(data_dir, file)        
    return image_path
