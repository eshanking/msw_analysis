#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 18:04:47 2022

@author: beckettpierce
"""

import numpy as np
import matplotlib.pyplot as plt
import math as m
import random
import seascapes_figures
import os
import imageio
np.random.seed(22)
random.seed(22)
class Sim():
    
    def __init__(self,
                 time_step = None,
                 diffusion_constant = None,
                 end_time = None,
                 x_step = None,
                 vessel_count = None,
                 lattice_length = None,
                 vessel_distribution = None,
                 **kwargs
                 ):
        self.time_step = time_step
        self.diffusion_constant = diffusion_constant
        self.end_time = end_time
        self.x_step = x_step
        self.vessel_count = vessel_count
        self.lattice_length = lattice_length
        self.vessel_distribution = vessel_distribution
        
    def gen_vessel_matrix(self):
        vessel_matrix = np.zeros((self.lattice_length,self.lattice_length))
        if self.vessel_distribution == 'random':
            xcords = random.sample(range(0,self.lattice_length),self.vessel_count)
            ycords = random.sample(range(0,self.lattice_length),self.vessel_count)
            for z in range(self.vessel_count):
                vessel_matrix[xcords[z],ycords[z]] = 1
        if self.vessel_distribution == 'fixed':
            spacing = 12
            for j in range(self.lattice_length):
                for i in range(self.lattice_length):
                    if ((((i==0) or (i%(2*spacing)==0)) and ((j+spacing))%(2*spacing)==0)):
                        vessel_matrix[i,j] =1
                    if ((((j==0) or (j%(2*spacing)==0)) and ((i+spacing))%(2*spacing)==0)):
                        vessel_matrix[i,j] =1
        return vessel_matrix
     #   if self.vessel_distribution == 'fixed':
            
        
        return vessel_matrix

    def run_sim(self):
        vessel_matrix = self.gen_vessel_matrix()
        self.vessel_matrix = vessel_matrix
        #initialize conc matrix as vessel matrix
        before_matrix = np.copy(vessel_matrix)
        after_matrix = np.copy(vessel_matrix)
        for iterator in range(int(self.end_time/self.time_step)):
            for xcord in range(self.lattice_length):
                for ycord in  range(self.lattice_length):
                    #interior points
                    # if xcord != 0 & ycord!= 0 & ycord!=self.lattice_length & xcord!=self.lattice_length & vessel_matrix[xcord,ycord]==0:
                    #     first_term = before_matrix[xcord,ycord]*(1-(4*self.diffusion_constant*self.time_step/(self.x_step**2)))
                    #     surrounding_sum = before_matrix[xcord+1,ycord] + before_matrix[xcord-1,ycord] + before_matrix[xcord,ycord+1] + before_matrix[xcord,ycord-1]
                    #     second_term = (self.diffusion_constant*self.time_step/(self.x_step**2))*surrounding_sum
                    #     after_matrix[xcord,ycord] = first_term + second_term
                    # #what if, x at max or min, y at max or min, both, vessel
                    # if xcord == 0 & ycord!=0 & vessel_matrix[xcord,ycord]==0:
                    #     first_term = before_matrix[xcord,ycord]*(1-(4*self.diffusion_constant*self.time_step/(self.x_step**2)))
                    #     surrounding_sum = before_matrix[xcord+1,ycord] + before_matrix[xcord-1,ycord] + before_matrix[xcord,ycord+1] +before_matrix[xcord,ycord+1]
                    #     second_term = (self.diffusion_constant*self.time_step/(self.x_step**2))*surrounding_sum
                    #     after_matrix[xcord,ycord] = first_term + second_term
                    # if xcord != 0 & ycord==0 & vessel_matrix[xcord,ycord]==0:
                    #     first_term = before_matrix[xcord,ycord]*(1-(4*self.diffusion_constant*self.time_step/(self.x_step**2)))
                    #     surrounding_sum = before_matrix[xcord+1,ycord] + before_matrix[xcord+1,ycord] + before_matrix[xcord,ycord+1] +before_matrix[xcord,ycord+1]
                    #     second_term = (self.diffusion_constant*self.time_step/(self.x_step**2))*surrounding_sum
                    #     after_matrix[xcord,ycord] = first_term + second_term
                    # if xcord == self.lattice_length & ycord != self.lattice_length:
                    # if xcord =! self.lattice_length & ycord == self.lattice_length:
                    if vessel_matrix[xcord,ycord]!=1:
                        cumulative_sum = 0
                        if xcord == 0:
                            cumulative_sum = cumulative_sum + 2*before_matrix[(xcord+1),ycord]
                        if ycord == 0:
                            cumulative_sum = cumulative_sum + 2*before_matrix[xcord,(ycord+1)]
                        if xcord == self.lattice_length-1:
                            cumulative_sum = cumulative_sum + 2*before_matrix[(xcord-1),ycord]
                        if ycord == self.lattice_length-1:
                            cumulative_sum = cumulative_sum + 2*before_matrix[xcord,(ycord-1)]
                        if self.lattice_length - 1 > xcord & xcord > 0:
                             cumulative_sum = cumulative_sum + before_matrix[xcord+1,ycord] + before_matrix[(xcord-1),ycord]
                        if self.lattice_length - 1 > ycord & ycord > 0:
                             cumulative_sum = cumulative_sum + before_matrix[xcord,(ycord+1)] + before_matrix[xcord,(ycord-1)]
                        first_term = before_matrix[xcord,ycord]*(1-(4*self.diffusion_constant*self.time_step/(self.x_step**2)))
                        surrounding_sum = cumulative_sum
                        second_term = (self.diffusion_constant*self.time_step/(self.x_step**2))*surrounding_sum
                        # if first_term!=0:
                        #     print(first_term)
                            
                        #print(second_term)
                        after_matrix[xcord,ycord] = first_term + second_term                        
                    else:
                        after_matrix[xcord,ycord]=1
                        cumulative_sum=0
            before_matrix = np.copy(after_matrix)        
                    
        return after_matrix




#i = Sim.gen_vessel_matrix()
time_list = np.linspace(0,200,30)
filenames = []
for x in time_list:
    
    np.random.seed(22)
    random.seed(22)
    test = Sim(      time_step = .25,
                     diffusion_constant = .01,
                     end_time = x,
                     x_step = .1,
                     vessel_count = 22,
                     lattice_length = 100,
                     vessel_distribution = 'fixed',
                     )
    result = test.run_sim()
    
    fig, ax = plt.subplots()
    im = ax.imshow(result, cmap=plt.get_cmap('hot'), interpolation='nearest',
                   vmin=0, vmax=1)
    ax.set(xlabel='x length', ylabel='y length', ylim = (0,100), xlim =  (0,100),
                   title='Steady State Drug Diffusion w/Constant Source & Reaction')


    fig.colorbar(im)
    
    
    filename = f'{x}.png'
    filenames.append(filename)

    # save frame
    plt.savefig(filename)
    plt.close()

with imageio.get_writer('mygif.gif', mode='I') as writer:
    for filename in filenames:
        image = imageio.imread(filename)
        writer.append_data(image)
for filename in set(filenames):
    os.remove(filename)



      