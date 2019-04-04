import numpy as np
import os
from functions import calculateAdjacency
import matplotlib.pyplot as plt
from Population import Population
import scipy.signal as sig

class Landscape():
    """Landscape is a subclass of ArrayBase, and holds the landscape and
    permanent diffusion_array.

    In addition to the ArrayBase methods it has getLandscapeDiffusion, which
    returns the diffusion_array.
    """
    def __init__(self, landscape_file):
        """
        Class constructor
        Inputs: landscape_file <npy file>
        """
        # Load landscape and lower resolution
        self.gis = self.raster(landscape_file)-2
        self.gis = self.lowerRes()
        self.gis = self.lowerRes()
        self.gis = self.lowerRes()

        # Save shape
        self.size_y = self.gis.shape[0]
        self.size_x = self.gis.shape[1]

        # Convert the gis data to land and water (1s and 0s)
        self.gisToLand()

        # Convert this data to a coastline
        self.landToCoast()
        # Generate fertility arrays
        # Generate a coast of only 1s and 0s
        self.coastal_land =  np.where(self.coastToFert()>0, 1, 0).astype(np.int32)
        self.Fert()

    def landscapeFileToArray(self, filename):
        '''landscapeFileToArray takes in a filename and returns the landscape as as
        numpy array (datatype float64) of landscape heights.
        '''
        # Checks the filename exists
        if not os.path.isfile(filename):
            raise NameError("\"{}\" is not a valid landscape file path" .format(filename))

        with open(filename, 'r') as f:
            lines = f.readlines()
        i = 0
        j = 0
        check = 0
        for line in lines:
            if check == 0:
                size = line.split(',')
                x_min = int(size[0])
                x_max = int(size[1])
                y_max = int(size[2])
                y_min = int(size[3])
                grids_per_int = int(size[4])

                x_grids = (x_max - x_min) * grids_per_int + 1
                y_grids = (y_max - y_min) * grids_per_int + 1
                landscape = np.empty([y_grids, x_grids])
                check = 1
            else:
                landscape[i,j] = float(line.split(',')[2])

                j += 1
                if j == x_grids:
                    j = 0
                    i += 1

        # landscape1 = landscape[174:177,62:66]
        # landscape1 = landscape[10:77,215:300]

        # return(landscape1)
        return(landscape)

    def lowerRes(self):
        """ Lower the resolution of gis data. Simply saves every second array
            index.
        """
        gis1 = self.gis[0:self.gis.shape[0]:2, 0:self.gis.shape[1]:2]
        # gis2 = self.gis[1:self.gis.shape[0]:2, 1:self.gis.shape[1]:2]
        # gis3 = (gis1 + gis2)/2
        return gis1

    def cutLand(self):
        """ Cut the land into relevant section and recalculate coast and fertility.
        """

        self.gis = self.gis[0:80,200:340]
        self.size_y = self.gis.shape[0]
        self.size_x = self.gis.shape[1]
        self.gisToLand()
        self.landToCoast()
        self.Fert()

    def gisToLand(self):
        """ Initialises land array.
            1s for land, 0s for water.
        """
        self.land = np.zeros([self.size_y, self.size_x])
        for i in range(self.size_y):
            for j in range(self.size_x):
                if (self.gis[i,j] >= 0):
                    self.land[i,j] = 1

    def landToCoast(self):
        """Finds the landscape's coastline
        """
        adjacency = calculateAdjacency(self.land)
        adjacency = np.where(adjacency>0, 1, 0)
        land_inverse = np.subtract(\
                        np.ones([self.size_y, self.size_x]),self.land)
        self.coast = np.multiply(adjacency, land_inverse)

    def coastToFert(self):
        """ Turn the coastline into preliminary fishing fertility
        """
        coast_adjacency = calculateAdjacency(self.coast)
        fish_fert = np.multiply(coast_adjacency, self.land)
        return fish_fert

    def coastToFert2(self):
        """ Turn the preliminary fertility into the final fertility
        """

        # conv filter is the area around the considered array point that
        # will be summed to create the final fertility using a convolution

        conv_filter = np.ones(49).reshape(7,7)
        coast_to_fert = self.coastToFert()
        coast_fert_2 = sig.convolve2d(coast_to_fert, conv_filter, mode = 'same',\
                        boundary = 'fill', fillvalue = 0)
        coast_fert_2 = np.multiply(coast_fert_2, self.land)

        coast_fert_2 = coast_fert_2/(np.max(coast_fert_2))
        return coast_fert_2

    def Fert(self):
        """ Calculates the overall fish and farm fertility values
        """
        # Random farm fertility for altitudes under 150m
        self.fert_farm = np.zeros([self.size_y, self.size_x])
        for i in range(self.size_y):
            for j in range(self.size_x):
                if (self.gis[i,j] > 0) and \
                    (self.gis[i,j] < 150):
                    self.fert_farm[i,j] = 0.5 + np.random.rand()*0.5

        # Calculate fish fertility
        self.fert_fish = np.ones_like(self.land) * (self.coastToFert2()>0) *0.3\
        + self.coastToFert2()

    def fertToCap(self, pop_1, pop_2, fish_constant, farm_constant):
        """ Convert the fertilities into carrying capacities.
            fish constant is carrying capacity per fish fertility
            farm constant is carrying capacity per farm fertility
        """
        self.cap_1 = (self.fert_fish * pop_1.fish_exploitation \
                                        * fish_constant) +\
                    (self.fert_farm * pop_1.farm_exploitation \
                                        * farm_constant)
        self.cap_2 = (self.fert_fish * pop_2.fish_exploitation \
                                        * fish_constant) +\
                    (self.fert_farm * pop_2.farm_exploitation \
                                        * farm_constant)

        self.mean_1 = np.mean(self.cap_1[(self.cap_1>0)])
        self.mean_2 = np.mean(self.cap_2[(self.cap_2>0)])
        print(self.mean_1, self.mean_2)
        # print(np.mean(self.cap_1))
        # print(hello)

    def calcDiffusionArray(self):
        """Takes self.base_array and returns the diffusion_array.

        Each value in the diffusion_array is the sum of the four surrounding
        values in the land array.
        """
        array_with_border = \
            np.zeros((self.size_y+2, self.size_x+2))
        array_with_border[1:1+self.size_y,1:1+self.size_x]= \
            self.land
            # np.ones((self.size_y, self.size_x))
        diffusion_array = np.zeros((self.size_y, self.size_x))

        for i in range (self.size_y):
            for j in range (self.size_x):
                diffusion_array[i,j]= \
                    array_with_border[i,j+1] + array_with_border[i+2, j+1] +\
                    array_with_border[i+1,j] + array_with_border[i+1, j+2]
        return diffusion_array

    def raster(self, landscape_file):
        """ Load the array file
        """

        data = np.load(landscape_file)
        # select scandinavia
        scand = data[3400:4550, 2950:4710]
        return(scand)
