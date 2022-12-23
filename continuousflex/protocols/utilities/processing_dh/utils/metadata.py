"""
This file reads, cleans an 'xmd' file and 
create a dataframe and a numpy array that
contains only the normal modes amplitudes
"""

import re
import numpy as np
from math import cos, sin, radians
from .euler2quaternion import eul2quat
import torch

def header(path):
    num_chars = 20
    head = 0
    with open(path, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if len(line) < num_chars:
                head+=1
            else:
                pass
    return head

def read_file(path):
    head = header(path)
    f = open(path, 'r')
    file_list = f.readlines()
    column_names = [col.replace('\n', '') for col in file_list[4:head]]
    for i in range(head):
        file_list.pop(0)
    return file_list, column_names


def create_array(path, flag='nma'):
    file_list, column_names = read_file(path)
    columns = len(list(filter(None,re.split("\s|'", file_list[0]))))
    num_modes = columns-(len(column_names)-1)
    nma_index = column_names.index(' _nmaDisplacements')
    for i in range(num_modes):
        column_names.insert(nma_index, 'mode '+str(num_modes - i))
    column_names.remove(' _nmaDisplacements')
    img_index = column_names.index(' _image')
    rot_index = column_names.index(' _angleRot')
    tilt_index = column_names.index(' _angleTilt')
    psi_index = column_names.index(' _anglePsi')
    shiftx_index = column_names.index(' _shiftX')
    shifty_index = column_names.index(' _shiftY')
    
    for i in range(len(file_list)):
        file_list[i] = file_list[i].replace('\n', ' ')
        file_list[i] = list(filter(None, re.split("\s|'",file_list[i])))
    
    data_array=np.reshape(file_list,(len(file_list),columns))
    img_names=data_array[:,img_index]
    nm_amplitudes = data_array[:, nma_index: nma_index+num_modes].astype('float32') 
    nm_amplitudes, nma_min, nma_max = min_max(nm_amplitudes)
    angles = data_array[:, [rot_index, tilt_index, psi_index]].astype('float32')
    shifts = data_array[:, [shiftx_index, shifty_index]].astype('float32')
    shifts, shf_min, shf_max = min_max(shifts)
    quaternions = np.zeros((angles.shape[0], 4), dtype='float32')
    for i in range(len(angles)):
        quaternions[i,:] = eul2quat(angles, i)
    if flag=='nma':
        print("Number of Normal Modes detected is: ",num_modes)
        return nm_amplitudes, nma_min, nma_max, img_names
    elif flag=='ang':
        return quaternions, img_names
    elif flag=='shf':
        return shifts, shf_min, shf_max, img_names
    else:
        raise ValueError('Unknown flag, you must select nma for Normal mode amplitudes, ang for euler angles, shf for shifts (X and Y)')


def min_max(arr):
    _min = torch.min(arr, dim=0)
    _max = torch.max(arr, dim=0)
    arr = (arr - _min[0]) / (_max[0] - _min[0])
    return arr, _min[0], _max[0]


def standardization(arr, params: bool = False, num_modes: int = 3):
    _mean = []
    _mu = []
    num_params = 0
    if params:
        num_params = num_modes + 5
    else:
        num_params = num_modes
    for i in range(num_params):
        _mean.append(np.mean(arr[:, i]))
        _mu.append(np.std(arr[:, i]))
    for i in range(num_params):
        for j in range(len(arr)):
            tmp = arr[j, i]
            arr[j, i] = (tmp-_mean[i])/_mu[i]
    return arr, _mean, _mu

def reverse_min_max(arr, _min, _max):
    """
    This function rescale back the target values to its original range
    it rescaled it back and put it in a list then reshape it to an array
    of the same shape as the input
    Parameters
    ----------
    arr : numpy array float32
        a numpy array for example (100,3).

    Returns
    -------
    rescaled_output : numpy array float32
        rescaled_output: a numpy array of the same shape as input for
        example (100, 3).

    """
    return (arr * (_max - _min)) + _min

def reverse_standardization(arr, _mean, _mu, params: bool = False, num_modes:int =3):
    num_params = 0
    if params:
        num_params = num_modes + 5
    else:
        num_params = num_modes
    rescaled_list = []
    for i in range(len(arr)):
        for j in range(params):
            rescaled_list.append((arr[i,j]*_mu[j])+_mean[j])
    rescaled_output = np.array(rescaled_list).reshape((len(arr), num_params))
    return rescaled_output



def rotation_matrix(euler_angles):
    rot = radians(euler_angles[0])
    tilt = radians(euler_angles[1])
    psi = radians(euler_angles[2])
    rot_mat = np.array([[(cos(rot)*cos(tilt)*cos(psi))-(sin(rot)*sin(psi)),
                         (-cos(psi)*sin(rot))-(cos(rot)*cos(tilt)*sin(psi)), 
                         (cos(rot)*sin(tilt))],
                        [(cos(rot)*sin(psi))+(cos(tilt)*cos(psi)*sin(rot)),
                         (cos(rot)*cos(psi))-(cos(tilt)*sin(rot)*sin(psi)),
                          (sin(rot)*sin(tilt))],
                        [-cos(psi)*sin(tilt), sin(psi)*sin(tilt),cos(tilt)]])

    return rot_mat
