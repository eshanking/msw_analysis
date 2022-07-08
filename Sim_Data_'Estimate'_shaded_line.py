#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:25:51 2022

@author: beckettpierce
"""

import numpy as np
import matplotlib.pyplot as plt
from seascapes_figures.classes.population_class import Population
import seascapes_figures

options = {    'k_abs':.95,
    'k_elim':.00839,
    'max_dose':5,
    'n_timestep':4000,
    'timestep_scale':.25,
    'fitness_data':'estimate',
    'curve_type':'pulsed',
    'prob_drop':.6
    }

p = Population(**options)
p.plot_fitness_curves()


x=np.arange(1000)
u=np.zeros(1000)


#need number/ to find the conc space
u = p.drug_curve
uhat = u
most_fit_at_conc = np.zeros(np.size(u))
conc_space = u

for z in range(np.size(u)):
    conc = conc_space[z]
    p_fit_list = p.gen_fit_land(conc)
    most_fit_at_conc[z] = (np.argmax(p_fit_list))
    
cc = p.gen_color_cycler()
cc_dict = cc.by_key()
colors = cc_dict['color']

fig, ax = plt.subplots()
#ax.plot(x, uhat,color='k')
ax.set(xlabel='Time', ylabel='drug concentration (ug/mL)',
       title='Dominant Type over Time')
plt.yscale("log")

# for l in range(np.size(x)):
#     if 0 == np.all(most_fit_at_conc == most_fit_at_conc[l]):
#         optimal_at_x = most_fit_at_conc[l]
#         optimal_at_x = int(optimal_at_x)
#         color_of_optimal = colors[optimal_at_x]
#         for c in range(np.size(x)):
#             if c+l > (np.size(x)-1):
#                 break
#             else:
#                 tester = int(most_fit_at_conc[c+l])
#                 if optimal_at_x != tester:
#                     terminal_point = c+l-1
        
#         # if optimal_at_x > 8:
#         #      ax.fill_between(x,0,uhat,color=color_of_optimal)
#         # else:
#         #      ax.fill_between(x,0,uhat,color=color_of_optimal,alpha=.5)
#         ax.fill_between(x[l:terminal_point],0,uhat[l:terminal_point],color = color_of_optimal)
#     else:
#          optimal_at_x = most_fit_at_conc[l]
#          optimal_at_x = int(optimal_at_x)
#          color_of_optimal = colors[optimal_at_x]
#          # if optimal_at_x > 8:
#          #     ax.fill_between(x,0,uhat,color=color_of_optimal)
         
#          ax.fill_between(x,0,uhat,color=color_of_optimal,alpha=.5)
for l in range(np.size(x)):
    if 0 == np.all(most_fit_at_conc == most_fit_at_conc[l]):
        if l == 0:
            the_tester = 19
        else:
            the_tester = most_fit_at_conc[l-1]
        optimal_at_x = most_fit_at_conc[l]
        no_terminal_counter = 0
        optimal_at_x = int(optimal_at_x)
        color_of_optimal = colors[optimal_at_x]
        if most_fit_at_conc[l] != the_tester:
            for c in range(np.size(x)):
                    if c+l > (np.size(x)-1):
                        break
                    else:
                        tester = int(most_fit_at_conc[c+l])
                        if optimal_at_x != tester:
                            terminal_point = c+l-1

                            break
                        else:
                            no_terminal_counter=c
                            
            if no_terminal_counter<np.size(x)+1:
                if optimal_at_x > 8:
                        ax.plot(x[l:terminal_point],uhat[l:terminal_point],color = color_of_optimal,linestyle='--')
                else:
                        ax.plot(x[l:terminal_point],uhat[l:terminal_point],color = color_of_optimal)
            if most_fit_at_conc[l]==most_fit_at_conc[np.size(x)-1] :
                if np.all(most_fit_at_conc[l:np.size(x)-1]==most_fit_at_conc[l]):
                    if optimal_at_x > 8:
                            ax.plot(x[l:np.size(x)-1],uhat[l:np.size(x)-1],color = color_of_optimal,linestyle='--')
                    else:
                            ax.plot(x[l:np.size(x)-1],uhat[l:np.size(x)-1],color = color_of_optimal)
            
    
       # ax.plot(x[l:terminal_point],uhat[l:terminal_point],color = color_of_optimal)
    else:
         optimal_at_x = most_fit_at_conc[l]
         optimal_at_x = int(optimal_at_x)
         color_of_optimal = colors[optimal_at_x]
         # if optimal_at_x > 8:
         #     ax.fill_between(x,0,uhat,color=color_of_optimal)
         if optimal_at_x > 8:
                 ax.plot(x,uhat,color=color_of_optimal,alpha=.5,linestyle='--')
         else:
                 ax.plot(x,uhat,color=color_of_optimal,alpha=.5)

        
         
plt.show()


