"""
This file reads and displays SPIDER 
'single' image using PIL (pillow) module 

"""

# import PIL module
from PIL import Image
import glob as glob
import numpy as np
import matplotlib.pyplot as plt
#import mrcfile
from tqdm import tqdm
import torch
from pwem.emlib.image import ImageHandler

"""
def spi2array(f_name) -> object:
    spi_image = Image.open(f_name, 'r')
    spi_array = np.array(spi_image, dtype='float32')
    spi_array = normalize(spi_array)
    return spi_array
"""
def spi2array(f_name) -> object:
    spi_array = ImageHandler().read(f_name).getData()
    return spi_array


# read SPIDER dataset from directory
def read_from_directory(dataset_path) -> object:
    spi_dataset = []
    for img_path in glob.glob(dataset_path + '/*.spi'):
        spi_dataset.append(spi2array(img_path))
    return spi_dataset


# read SPIDER dataset from list
def read_from_list(xmd_path, img_name_list):
    spi_dataset = []
    for img_path in img_name_list:
        spi_dataset.append(spi2array(str(xmd_path) + '/' + img_path))
    return spi_dataset


# show SPIDER image format
def imshow(spi_image):
    plt.imshow(spi_image)
    plt.show()
    return

def normalize(spi_array):
    _sdv = np.std(spi_array)
    _mean = np.mean(spi_array)
    _min = np.min(spi_array)
    _max = np.max(spi_array)
    spi_array = (spi_array)/max(_max, np.abs(_min))
    #spi_array = (spi_array - _mean) / _sdv
    #spi_array = (spi_array - _min) / (_max - _min)
    return spi_array
    
    
def torch_normalize(spi_array):
    _sdv = torch.std(spi_array)
    _mean = torch.mean(spi_array)
    _min = torch.min(spi_array)
    _max = torch.max(spi_array)
    spi_array = (spi_array - _min) / (_max - _min)
    return spi_array

"""
def mrc_stack_reader(path):
    mrc = mrcfile.mmap(path, mode='r+')
    return mrc
"""
