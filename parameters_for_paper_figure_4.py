#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 28 18:04:47 2022
@author: beckettpierce
"""
import numpy as np
import matplotlib.pyplot as plt
import random
import os
import imageio
#import seascapes_figures
from seascapes_figures.classes.population_class import Population

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
                 max_dose = None,
                 **kwargs
                 ):
        
        self.time_step = time_step
        self.diffusion_constant = diffusion_constant
        self.end_time = end_time
        self.x_step = x_step
        self.vessel_count = vessel_count
        self.lattice_length = lattice_length
        self.vessel_distribution = vessel_distribution
        self.max_dose = max_dose
        
    def gen_vessel_matrix(self):
        vessel_matrix = np.zeros((self.lattice_length,self.lattice_length))
        if self.vessel_distribution == 'random':
            xcords = random.sample(range(0,self.lattice_length),self.vessel_count)
            ycords = random.sample(range(0,self.lattice_length),self.vessel_count)
            for z in range(self.vessel_count):
                vessel_matrix[xcords[z],ycords[z]] = self.max_dose
        if self.vessel_distribution == 'fixed':
            spacing = 12
            for j in range(self.lattice_length):
                for i in range(self.lattice_length):
                    if ((((i==0) or (i%(2*spacing)==0)) and ((j+spacing))%(2*spacing)==0)):
                        vessel_matrix[i,j] = self.max_dose
                    if ((((j==0) or (j%(2*spacing)==0)) and ((i+spacing))%(2*spacing)==0)):
                        vessel_matrix[i,j] = self.max_dose
        return vessel_matrix
     #   if self.vessel_distribution == 'fixed':
            
        
        return vessel_matrix

    def run_sim(self):
        vessel_matrix = self.gen_vessel_matrix()
        self.vessel_matrix = vessel_matrix
        #initialize conc matrix as vessel matrix
        before_matrix = np.copy(vessel_matrix)
        after_matrix = np.copy(vessel_matrix)
        for iterator in range(int(self.end_time / self.time_step)):
            print(iterator)
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
                        
                    if vessel_matrix[xcord,ycord] != self.max_dose:
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
                        first_term = before_matrix[xcord,ycord] * (1- (4*self.diffusion_constant* self.time_step / (self.x_step**2)))
                        surrounding_sum = cumulative_sum
                        second_term = (self.diffusion_constant * self.time_step / (self.x_step ** 2)) *surrounding_sum
                        # if first_term!=0:
                        #     print(first_term)
                            
                        #print(second_term)
                        after_matrix[xcord,ycord] = first_term + second_term                        
                    else:
                        after_matrix[xcord,ycord] = self.max_dose
                        cumulative_sum = 0
            before_matrix = np.copy(after_matrix)        
                    
        return after_matrix

#i = Sim.gen_vessel_matrix()
#time_list = np.linspace(0,40,20)

time_list = np.linspace(0,1000,20)
filenames = []
options = {    'k_abs':.95,
    'k_elim':.00839,
    'max_dose':5,
    'n_timestep':4000,
    'timestep_scale':.25,
    'fitness_data':'estimate',
    'curve_type':'pulsed',
    'prob_drop':.4
    }

p = Population(**options)
p.plot_fitness_curves()
cc = p.gen_color_cycler()
cc_dict = cc.by_key()
colors = cc_dict['color']
counter = np.zeros((np.size(time_list),16))
for x in time_list:
    
    np.random.seed(22)
    random.seed(22)
    test = Sim(      time_step = 1,
                     diffusion_constant = 3.64 * 10**-6 ,
                     #antibiotic diffusion
                     end_time = x,
                     x_step = .005,
                     #50 microns
                     vessel_count = 22,
                     lattice_length = 100,
                     vessel_distribution = 'fixed',
                     max_dose = 5
                     )
    
    result = test.run_sim()
    #color gen
    most_fit_at_conc = np.zeros([test.lattice_length,test.lattice_length])
    color_matrix = np.zeros([test.lattice_length,test.lattice_length,3])
    most_fit_at_conc = np.zeros([test.lattice_length,test.lattice_length])
    color_matrix = np.zeros([test.lattice_length,test.lattice_length,3])

    
    fig, ax = plt.subplots()

    for result_x in range(test.lattice_length):
        for result_y in range(test.lattice_length):
            
            p_fit_list = p.gen_fit_land(result[result_x,result_y])
            most_fit_at_conc[result_x,result_y] = (np.argmax(p_fit_list))
            
            optimal_at_x = most_fit_at_conc[result_x,result_y]
            optimal_at_x = int(optimal_at_x)
            color_of_optimal = colors[optimal_at_x]

           # counter[time_list[int(x)]][optimal_at_x] =+1 
            location_on_time_list = np.where(time_list ==x)[0][0]
            counter[location_on_time_list,optimal_at_x] +=1
            for z in range(3):
                color_matrix[result_x,result_y,z]= color_of_optimal[z]
    #         im = ax.scatter(result_x,result_y, color=color_matrix[result_x,result_y])
    #         ax.set(xlabel='x length', ylabel='y length', ylim = (0,100), xlim =  (0,100),
    #                title='Steady State Drug Diffusion w/Constant Source & Reaction')
    
    # filename = f'{x}.png'
    # filenames.append(filename)

    # # save frame
    # plt.savefig(filename)
    # plt.close()

# with imageio.get_writer('tester.gif', mode='I') as writer:
#     for filename in filenames:
#         image = imageio.imread(filename)
#         writer.append_data(image)
# for filename in set(filenames):
#     os.remove(filename)

#fig2, ax2 = plt.subplots()

counter = counter/10000

for jj in range(p.n_allele**2):
    if jj>8:
        ax.plot(time_list,counter[:,jj],color =colors[jj],linestyle='--')
    else:
        ax.plot(time_list,counter[:,jj],color =colors[jj])
ax.set(xlabel='Time', ylabel='Population %',title='Share of Total Population over Time')





