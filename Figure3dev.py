#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 14 13:30:19 2022

@author: beckettpierce
"""

import numpy as np
import matplotlib.pyplot as plt
from seascapes_figures.classes.population_class import Population
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import random


seed = 9345
np.random.seed(seed)
random.seed(seed)
x = np.arange(100)
u = np.zeros(100)
def steadystate(x,k,D,r):
    for j in range(len(x)):
        u[j] = (k/np.sqrt(4*D*r))*(np.e**(-1*np.abs(x[j])*np.sqrt(r/D)))
    return u

uhat = steadystate(x=x,k=100,D=6.45,r=.1)
uhat =  uhat/np.max(uhat)
#k in ug/ML
#D in 10^-6 cm^2/s
fig, ax = plt.subplots(1,2,figsize=(14, 5))
# ax[0].semilogy(x, uhat)
# ax[0].set(xlabel='x (10^-3 cm)', ylabel='% of Max Dose')


options = {    'k_abs':.95,
    'k_elim':.00839,
    'max_dose':5,
    'n_timestep':4000,
    'timestep_scale':.25,
    'fitness_data':'random',
    'curve_type':'pulsed',
    'prob_drop':.5,
    'n_allele':2
    }


p = Population(**options)


# p.ic50 = [-3.78766218, -4.02237296, -5.22221995, -3.09765892]
# p.drugless_rates = [1.22802236, 1.14399848, 1.28949852, 0.93619847]

x=np.arange(1000)
u=np.zeros(4000)


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
#ax.plot(x, uhat,color='k')


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
                        ax[0].plot(x[l:terminal_point],uhat[l:terminal_point],color = color_of_optimal,linestyle='--',linewidth=3)
                else:
                        ax[0].plot(x[l:terminal_point],uhat[l:terminal_point],color = color_of_optimal,linewidth=3)
            if most_fit_at_conc[l]==most_fit_at_conc[np.size(x)-1] :
                if np.all(most_fit_at_conc[l:np.size(x)-1]==most_fit_at_conc[l]):
                    if optimal_at_x > 8:
                            ax[0].plot(x[l:np.size(x)-1],uhat[l:np.size(x)-1],color = color_of_optimal,linestyle='--',linewidth=3)
                    else:
                            ax[0].plot(x[l:np.size(x)-1],uhat[l:np.size(x)-1],color = color_of_optimal,linewidth=3)
            
    
       # ax.plot(x[l:terminal_point],uhat[l:terminal_point],color = color_of_optimal)
    else:
         optimal_at_x = most_fit_at_conc[l]
         optimal_at_x = int(optimal_at_x)
         color_of_optimal = colors[optimal_at_x]
         # if optimal_at_x > 8:
         #     ax.fill_between(x,0,uhat,color=color_of_optimal)
         if optimal_at_x > 8:

                 ax[0].plot(x,uhat,color=color_of_optimal,alpha=.5,linestyle='--',linewidth=3)
         else:
                 ax[0].plot(x,uhat,color=color_of_optimal,alpha=.5,linewidth=3)


ax[0].set(xlabel='Time')
ax[0].set(ylabel='Drug concentration (ug/mL)')


x = np.logspace(-2,2,num=1000)
u = np.zeros(1000)
def steadystate(x,k,D,r):
    for j in range(np.size(x)):
        u[j] = (k/np.sqrt(4*D*r))*(np.e**(-1*np.abs(x[j])*np.sqrt(r/D)))
    return u

uhat = steadystate(x=x,k=100,D=6.45,r=.1)
#k in ug/ML
#D in 10^-6 cm^2/s

np.random.seed(seed)
random.seed(seed)

#now goal is to look at where conc space dominance and what is dominant where
most_fit_at_conc = np.zeros(np.size(u))
conc_space = u
options = {    'k_abs':.95,
    'k_elim':.00839,
    'max_dose':5,
    'n_timestep':4000,
    'timestep_scale':.25,
    'fitness_data':'random',
    'curve_type':'pulsed',
    'prob_drop':.5,
    'n_allele':2
    }
# p.ic50 = [-3.78766218, -4.02237296, -5.22221995, -3.09765892]
# p.drugless_rates = [1.22802236, 1.14399848, 1.28949852, 0.93619847]


p = Population(**options)
p.plot_fitness_curves()
for z in range(np.size(u)):
    conc = conc_space[z]
    p_fit_list = p.gen_fit_land(conc)
    most_fit_at_conc[z] = (np.argmax(p_fit_list))
##now most_fit_at_conc is our list of what is most fit at each point, now we must associate colors with each of these
cc = p.gen_color_cycler()
cc_dict = cc.by_key()
colors = cc_dict['color']



#we have list of colors for each thing, where we want each color
#now need to chop up x & uhat into different arrays based on where there is this optimal conc
# uhat = uhat/np.max(uhat)
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
                            
            if no_terminal_counter<1001:
                if optimal_at_x > 8:
                        ax[1].semilogy(x[l:terminal_point],uhat[l:terminal_point],color = color_of_optimal,linestyle='--',linewidth=3)
                else:
                        ax[1].semilogy(x[l:terminal_point],uhat[l:terminal_point],color = color_of_optimal,linewidth=3)
            if most_fit_at_conc[l]==most_fit_at_conc[np.size(x)-1] :
                if np.all(most_fit_at_conc[l:np.size(x)-1]==most_fit_at_conc[l]):
                    if optimal_at_x > 8:
                            ax[1].semilogy(x[l:np.size(x)-1],uhat[l:np.size(x)-1],color = color_of_optimal,linestyle='--',linewidth=3)
                    else:
                            ax[1].semilogy(x[l:np.size(x)-1],uhat[l:np.size(x)-1],color = color_of_optimal,linewidth=3)
    
       # ax.plot(x[l:terminal_point],uhat[l:terminal_point],color = color_of_optimal)
    else:
         optimal_at_x = most_fit_at_conc[l]
         optimal_at_x = int(optimal_at_x)
         color_of_optimal = colors[optimal_at_x]
         # if optimal_at_x > 8:
         #     ax.fill_between(x,0,uhat,color=color_of_optimal)
         if optimal_at_x > 8:
                 ax[1].semilogy(x,uhat,color=color_of_optimal,alpha=.5,linestyle='--',linewidth=3)
         else:
                 ax[1].semilogy(x,uhat,color=color_of_optimal,alpha=.5,linewidth=3)
                 
                 
legend_elements = [Line2D([0], [0], color='blue', lw=2, label='solid',linestyle='-'), 
                   Line2D([0], [0],  color='blue',lw=2, label='dot',linestyle='--'),
                   Line2D([0], [0],  color='orange',lw=2, label='dot',linestyle='-'),
                   Line2D([0], [0],  color='orange',lw=2, label='dot',linestyle='--'),
                   Line2D([0], [0],  color='green',lw=2, label='dot',linestyle='-'),
                   Line2D([0], [0],  color='green',lw=2, label='dot',linestyle='--'),
                   Line2D([0], [0],  color='red',lw=2, label='dot',linestyle='-'),
                   Line2D([0], [0],  color='red',lw=2, label='dot',linestyle='--'),
                   Line2D([0], [0],  color='purple',lw=2, label='dot',linestyle='-'),
                   Line2D([0], [0],  color='purple',lw=2, label='dot',linestyle='--'),
                   Line2D([0], [0],  color='brown',lw=2, label='dot',linestyle='-'),
                   Line2D([0], [0],  color='brown',lw=2, label='dot',linestyle='--'),
                   Line2D([0], [0],  color='pink',lw=2, label='dot',linestyle='-'),
                   Line2D([0], [0],  color='pink',lw=2, label='dot',linestyle='--'),
                   Line2D([0], [0],  color='yellow',lw=2, label='dot',linestyle='-'),
                   Line2D([0], [0],  color='yellow',lw=2, label='dot',linestyle='--'),
                   Line2D([0], [0],  color='gray',lw=2, label='dot',linestyle='-'),
                   Line2D([0], [0],  color='gray',lw=2, label='dot',linestyle='--')]
                          
# Create the figure
ax[1].legend(handles=legend_elements, loc='upper right',frameon=False)



ax[1].set(xlabel='x (10^-3 cm)', ylabel='Drug Concentration (ug/ml)')