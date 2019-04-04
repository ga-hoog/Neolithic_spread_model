import numpy as np
import matplotlib.pyplot as plt
from Landscape import Landscape
from Population import Population
from functions import calculateAdjacency, calculateDifAdjacency, calculateFlight
from matplotlib.animation import FuncAnimation
import scipy.signal as sig
import math
import time
import os
# from trade import calctrade
from trade import calctrade
from boat import boatflight, boatflight2
from flight import calcnj

class Simulation:
    """Simulation class, takes in the landscape and two population objects.
        Updates the populations using timeStep1 function.
    """
    def __init__(self, pop_1, pop_2, landscape, tg=50, tk=10,\
                D=0.7, Df=0.05, Dtrade = 0.15, dt=0.05, lflight=2, ltrade=4):

        # Save the relevant values from population and landscape objects
        self.pop_1 = pop_1.population
        self.pop_2 = pop_2.population
        self.trade_1_array = pop_1.trade_array
        self.trade_2_array = pop_2.trade_array
        self.land = landscape.land
        self.coastal_land = landscape.coastal_land
        self.coast = landscape.coast
        self.cap_1 = landscape.cap_1
        self.cap_2 = landscape.cap_2
        self.mean_cap_1 = landscape.mean_1
        self.mean_cap_2 = landscape.mean_2
        self.size_x = landscape.size_x
        self.size_y = landscape.size_y
        self.tg = tg
        self.tk = tk
        self.D = D
        self.Df = Df
        self.Dtrade = Dtrade
        self.dt = dt
        self.l = lflight
        self.ltrade = ltrade

        # Generate filters used in trade and flight calculations
        self.flight_filter = self.generateFlightFilter()
        self.trade_filter = self.generateTradeFilter()

        self.timeStep1()


    def calculateGrowth(self, pop_1_array, pop_2_array,\
                                cap_1_array, cap_2_array):
        """ Calculate the growth
        """
        N = np.add(pop_1_array, pop_2_array)

        growth_1 = (pop_1_array*(1-np.divide(N,cap_1_array,\
                    out=np.ones_like(pop_1_array)*100, where = cap_1_array != 0)))\
                    /self.tg

        growth_2 = (pop_2_array*(1-np.divide(N,cap_2_array,\
                    out=np.ones_like(pop_2_array)*100, where = cap_2_array != 0)))\
                    /self.tg
        return growth_1, growth_2

    def calculateConflictDeathRates(self, pop_1_array, pop_2_array):
        """ Calculate the conflict death rates
        """
        encounter_ratio = -pop_1_array*pop_2_array/\
                            (self.mean_cap_1*self.mean_cap_2)
        death_1 = encounter_ratio*pop_2_array/self.tk
        death_2 = encounter_ratio*pop_1_array/self.tk
        return death_1, death_2

    def calculateDiffusion(self, pop_array, cap_array):
        """ Calculate the diffusion of a single population
        """
        # Convert into densities
        densities = np.divide(pop_array, cap_array,\
                    out=np.zeros_like(pop_array), where = cap_array != 0)

        adjacency = calculateDifAdjacency(densities, cap_array)

        diffusion = -self.D * np.multiply(adjacency, cap_array)
        return diffusion

    def calculateFilterSize(self, l):
        """ Calculate the correct filter_size depending on the
            trade/flight distance parameter
        """
        found = False
        i = 1
        while found == False:
            if np.exp(-i/l) < 0.01:
                found = True
                filter_size = i
            else:
                i += 1
        if filter_size > 60:
            filter_size = 61
        if filter_size % 2 == 0:
            filter_size += 1
        return(filter_size)

    def generateFlightFilter(self):
        """ Generate the filter used in the flight calculation
        """
        filter_size = self.calculateFilterSize(self.l)
        fil = np.zeros((filter_size,filter_size))
        half_size = int(filter_size/2)

        for i in range(filter_size):
            for j in range(filter_size):
                fil[i,j] = np.exp(-np.sqrt(((half_size-i)*(half_size-i) +\
                (half_size-j)*(half_size-j)))/self.l)
        return fil

    def generateTradeFilter(self):
        """ Generate the filter used in the trade calculation
        """
        filter_size = self.calculateFilterSize(self.ltrade)
        fil = np.zeros((filter_size,filter_size))
        half_size = int(filter_size/2)
        for i in range(filter_size):
            for j in range(filter_size):
                fil[i,j] = self.Dtrade * np.exp(-np.sqrt(((half_size-i)*(half_size-i) +(half_size-j)*(half_size-j)))/self.ltrade)
        fil[half_size, half_size] = 0
        return fil

    def calculateNj(self, pop_1_array, pop_2_array):
        """ CalculateNj
        """
        Nj_1 = np.zeros_like(pop_1_array)
        Nj_2 = np.zeros_like(pop_1_array)
        Nj_1 = np.asfortranarray(Nj_1, dtype= np.float32)
        Nj_2 = np.asfortranarray(Nj_2, dtype= np.float32)
        calcnj(pop_1_array.shape[0],pop_1_array.shape[1], pop_1_array, pop_2_array,\
                self.flight_filter, Nj_1, Nj_2, self.land)
        Nj_1 = (Nj_1/(2*math.pi*(self.l)**2))
        Nj_2 = (Nj_2/(2*math.pi*(self.l)**2))
        return Nj_1, Nj_2

    def calculateFlight(self, pop_1_array, pop_2_array, cap_1, cap_2):
        """ Calculate the flight
        """
        Nj_1, Nj_2 = self.calculateNj(pop_1_array, pop_2_array)
        Nj_1 = np.where(Nj_1<(pop_1_array*0.25/(self.Df)), Nj_1, (pop_1_array*0.25/(self.Df)))
        Nj_2 = np.where(Nj_2<(pop_2_array*0.25/(self.Df)), Nj_2, (pop_2_array*0.25/(self.Df)))

        self.Nj_1 = Nj_1
        self.Nj_2 = Nj_2

        flight_1 = calculateFlight(Nj_1 ,cap_1)
        flight_2 = calculateFlight(Nj_2 ,cap_2)

        flight_1 = -self.Df * flight_1
        flight_2 = -self.Df * flight_2
        return flight_1, flight_2

    def timeStep1(self):
        """ Update the populations
        """
        # update capacities
        cap_1 = self.cap_1 + self.trade_1_array
        cap_2 = self.cap_2 + self.trade_2_array
        self.cap_trade_1 = cap_1
        self.cap_trade_2 = cap_2

        diffusion_1 = self.calculateDiffusion(self.pop_1, cap_1)
        diffusion_2 = self.calculateDiffusion(self.pop_2, cap_2)
        self.diffusion_1 = diffusion_1*self.dt
        self.diffusion_2 = diffusion_2*self.dt

        flight_1, flight_2 = self.calculateFlight(self.pop_1, self.pop_2,\
                                                    cap_1, cap_2)
        self.flight_1 = flight_1 *self.dt
        self.flight_2 = flight_2 *self.dt

        growth_1, growth_2 = self.calculateGrowth(self.pop_1, self.pop_2,\
                                                cap_1, cap_2)

        self.growth_1 = growth_1 * self.dt
        self.growth_2 = growth_2 * self.dt

        deathrate_1, deathrate_2 = self.calculateConflictDeathRates(\
                                                        self.pop_1, self.pop_2)


        self.deathrate_1 = np.where(abs(deathrate_1*self.dt)<self.pop_1,\
                            deathrate_1*self.dt, -self.pop_1)
        self.deathrate_2 = np.where(abs(deathrate_2*self.dt)<self.pop_2, \
                            deathrate_2*self.dt, -self.pop_2)

        pop_1_delta = self.dt * (growth_1 + deathrate_1 + diffusion_1 + flight_1)
        pop_2_delta = self.dt * (growth_2 + deathrate_2 + diffusion_2 + flight_2)
        self.pop_1 += pop_1_delta
        self.pop_2 += pop_2_delta
        self.pop_1 = np.where(self.pop_1 >= 0, self.pop_1, 0)
        self.pop_2 = np.where(self.pop_2 >= 0, self.pop_2, 0)


        # Calculate the new trade arrays
        trade_1_del_test = np.zeros_like(self.trade_1_array)
        trade_1_del_test = np.asfortranarray(trade_1_del_test, dtype = np.float32)
        trade_1_args = -np.ones((2, self.size_y, self.size_x))
        self.trade_1_args = np.asfortranarray(trade_1_args, dtype = np.float32)
        trade_2_args = -np.ones((2, self.size_y, self.size_x))
        self.trade_2_args = np.asfortranarray(trade_2_args, dtype = np.float32)
        calctrade(trade_1_del_test.shape[0], trade_1_del_test.shape[1], self.pop_1,\
            self.cap_1, self.trade_filter, trade_1_del_test, self.trade_1_args)
        trade_2_del_test = np.zeros_like(self.trade_2_array)
        trade_2_del_test = np.asfortranarray(trade_2_del_test, dtype = np.float32)
        calctrade(trade_2_del_test.shape[0], trade_2_del_test.shape[1], self.pop_2,\
            self.cap_2, self.trade_filter, trade_2_del_test, self.trade_2_args)
        trade_1_del_test = np.where(self.cap_1 ==0, 0, trade_1_del_test)
        trade_2_del_test = np.where(self.cap_2 ==0, 0, trade_2_del_test)

        self.trade_1_array = trade_1_del_test
        self.trade_2_array = trade_2_del_test


        self.pop_1 = np.asfortranarray(self.pop_1, dtype = np.float32)
        self.pop_2 = np.asfortranarray(self.pop_2, dtype = np.float32)
        boatflight2(trade_1_del_test.shape[0], trade_1_del_test.shape[1],\
                    30, 0.1, 0.05,\
                    np.asfortranarray(self.coastal_land, dtype = np.int32),\
                    np.asfortranarray(cap_1, dtype = np.float32),\
                    self.pop_1, self.pop_2)
        boatflight(trade_1_del_test.shape[0], trade_1_del_test.shape[1],\

                    # 20, 0.1, 0.005,\
                    30, 0.5, 0.0003,\
                    np.asfortranarray(self.coastal_land, dtype = np.int32),\
                    np.asfortranarray(cap_2, dtype = np.float32),\
                    self.pop_2)

    def animate1(self, frame):
        self.timeStep1()
        print(frame)
        self.matrice1.set_array(self.pop_1+self.coast*45)
        self.matrice2.set_array(self.pop_2+self.coast*90)
        self.matrice3.set_array(self.trade_1_array)
        self.matrice4.set_array(self.trade_2_array)
        return [self.matrice1, self.matrice2, self.matrice3, self.matrice4]

    def display1(self, frames):
        fig, axes = plt.subplots(2,2)
        ax1=  axes[0,0]
        ax2=  axes[0,1]
        ax3=  axes[1,0]
        ax4=  axes[1,1]

        self.matrice1 = ax1.matshow(self.pop_1, vmax = 5)
        self.matrice2 = ax2.matshow(self.pop_2, vmax = 5)
        self.matrice3 = ax3.matshow(-self.Nj_1, vmax = 10, vmin = -10)
        self.matrice4 = ax4.matshow(-self.Nj_2, vmax = 10, vmin = -10)

        fig.colorbar(self.matrice1, ax=ax1)
        fig.colorbar(self.matrice2, ax=ax2)
        fig.colorbar(self.matrice3, ax=ax3)
        fig.colorbar(self.matrice4, ax=ax4)
        ax1.axis('scaled')
        ax2.axis('scaled')
        ax3.axis('scaled')
        ax4.axis('scaled')
        anim = FuncAnimation(fig, self.animate1, frames\
                , repeat = False, interval = 0, blit = True)
        plt.show()

def initialiseSimulation(fish_set, farm_set, fish_exploit, farm_exploit,\
                        lflight, ltrade, tg=30.0,\
                        tk=10.0, D=0.7, Df=0.05, Dtrade = 0.015, dt=0.05):

    land = Landscape("europe.npy")
    pop1 = Population(land, fish_set, fish_exploit, 0)
    # pop1 = Population(land, fish_set, 0, fish_exploit)
    pop2 = Population(land, farm_set, 0, farm_exploit)
    land.fertToCap(pop1, pop2, 1, 1)
    pop1.population = pop1.population * (land.cap_1>0)
    pop2.population = pop2.population * (land.cap_2>0)
    sim = Simulation(pop1, pop2, land, tg, tk, D, Df, \
                    Dtrade, dt, lflight, ltrade)
    parameters = [fish_exploit, farm_exploit, tg, tk, D, \
                    Df, Dtrade, dt, lflight, ltrade]

    return sim, parameters

def presentSimulation():
    sim, parameters = initialiseSimulation('basic_fish','basic_farm', 25, 60,2,4)
    sim.display1(100000)
#
if __name__ == "__main__":
    presentSimulation()
