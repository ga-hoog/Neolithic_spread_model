import numpy as np

class Population():
    """
    """
    def __init__(self, landscape, init_setting,\
                fish_exploitation, farm_exploitation):
        """
        Class constructor
        Inputs:
                :landscape:      <landscape object>
                :init_setting:   <string>
                :fish_exploitation: <float>
                :farm_exploitation: <float>

        """
        self.farm_exploitation = farm_exploitation
        self.fish_exploitation = fish_exploitation
        self.population = self.initialisePopulation(landscape, init_setting)
        self.trade_array = np.zeros_like(landscape.land)


    def initialisePopulation(self, landscape, init_setting):
        """ Initialise the population for different init settings
        """

        if init_setting == 'zero':
            population = np.zeros((landscape.size_y, landscape.size_x))

        if init_setting == 'basic_fish':
            population = np.zeros((landscape.size_y, landscape.size_x))
            for i in range(landscape.size_y):
                for j in range(landscape.size_x):
                    if (landscape.fert_fish[i,j] > 0 and np.random.rand()<0.12):
                        population[i,j]= np.random.rand()


        if init_setting == 'basic_farm':
            population = np.zeros((landscape.size_y, landscape.size_x))
            population[landscape.size_y-3:landscape.size_y-1, 75:125] = \
            0.1*(0.2 + 0.3*np.random.rand())
        return population
