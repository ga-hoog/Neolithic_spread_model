import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import matplotlib.animation as animation
import sys
import json
import numpy as np
import sys, os


def calcDiffusionArray(self):
    """Takes self.base_array and returns the diffusion_array.

    Each value in the diffusion_array is the sum of the four surrounding
    values in the self.base_array. This is used in the diffusion part of
    the animal reaction/diffusion algorithm.
    """
    array_with_border = \
        np.zeros((self.size_y+2, self.size_x+2))
    array_with_border[1:1+self.size_y,1:1+self.size_x]= \
        self.base_array
    diffusion_array = \
        np.zeros((self.size_y, self.size_x))

    for i in range (self.size_y):
        for j in range (self.size_x):
            diffusion_array[i,j]= \
                array_with_border[i,j+1] + array_with_border[i+2, j+1] +\
                array_with_border[i+1,j] + array_with_border[i+1, j+2]
    return diffusion_array

def calculateAdjacency(array):
    """ Calculates the adjacancy values of an array
    """
    shift0 = np.roll(array, 1, 0)
    shift0[0,:] = 0
    shift1 = np.roll(array, 1, 1)
    shift1[:,0] = 0
    shift2 = np.roll(array, -1, 0)
    shift2[array.shape[0]-1, :] = 0
    shift3 = np.roll(array, -1, 1)
    shift3[:, array.shape[1]-1] = 0
    adjacency = shift0+shift1+shift2+shift3
    return adjacency

def calculateFlight(array, cap_array):
    """ Calculates the flight
    """

    land_array = np.where(cap_array>0, 1, 0)

    shift0 = np.roll(array, 1, 0)
    shift0[0,:] = 0
    shift1 = np.roll(array, 1, 1)
    shift1[:,0] = 0
    shift2 = np.roll(array, -1, 0)
    shift2[array.shape[0]-1, :] = 0
    shift3 = np.roll(array, -1, 1)
    shift3[:, array.shape[1]-1] = 0

    shift_land_0 = np.roll(land_array, 1, 0)
    shift_land_0[0,:] = 0
    shift_land_1 = np.roll(land_array, 1, 1)
    shift_land_1[:,0] = 0
    shift_land_2 = np.roll(land_array, -1, 0)
    shift_land_2[land_array.shape[0]-1,:] = 0
    shift_land_3 = np.roll(land_array, -1, 1)
    shift_land_3[:, land_array.shape[1]-1] = 0

    flight = land_array*(shift_land_0*(array-shift0)+\
                        shift_land_1*(array-shift1)+\
                        shift_land_2*(array-shift2)+\
                        shift_land_3*(array-shift3))
    return flight

def calculateDifAdjacency(array, cap_array):
    """ Calculates diffusion
    """
    # Rolled population arrays
    shift0 = np.roll(array, 1, 0)
    shift0[0,:] = 0
    shift1 = np.roll(array, 1, 1)
    shift1[:,0] = 0
    shift2 = np.roll(array, -1, 0)
    shift2[array.shape[0]-1, :] = 0
    shift3 = np.roll(array, -1, 1)
    shift3[:, array.shape[1]-1] = 0


    # Rolled cap arrays
    shift_cap0 = np.roll(cap_array, 1, 0)
    shift_cap0[0,:] = 0
    shift_cap1 = np.roll(cap_array, 1, 1)
    shift_cap1[:,0] = 0
    shift_cap2 = np.roll(cap_array, -1, 0)
    shift_cap2[cap_array.shape[0]-1, :] = 0
    shift_cap3 = np.roll(cap_array, -1, 1)
    shift_cap3[:, cap_array.shape[1]-1] = 0

    sum0 = cap_array+shift_cap0
    ratio0 = np.divide(shift_cap0, sum0,\
        out=np.zeros_like(shift_cap0), where = sum0 != 0)
    sum1 = cap_array+shift_cap1
    ratio1 = np.divide(shift_cap1, sum1,\
        out=np.zeros_like(shift_cap1), where = sum1 != 0)
    sum2 = cap_array+shift_cap2
    ratio2 = np.divide(shift_cap2, sum2,\
        out=np.zeros_like(shift_cap2), where = sum2 != 0)
    sum3 = cap_array+shift_cap3
    ratio3 = np.divide(shift_cap3, sum3,\
        out=np.zeros_like(shift_cap3), where = sum3 != 0)

    adjacency = ratio0*(array - shift0) +\
                ratio1*(array - shift1) +\
                ratio2*(array - shift2) +\
                ratio3*(array - shift3)
    return adjacency
