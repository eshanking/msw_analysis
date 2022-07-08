#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 13:49:21 2022

@author: beckettpierce
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 12:18:41 2022

@author: beckettpierce
"""

import numpy as np
import matplotlib.pyplot as plt
from seascapes_figures.classes.population_class import Population
import seascapes_figures


x = np.logspace(-2,2,num=1000)
u = np.zeros(1000)
def steadystate(x,k,D,r):
    for j in range(np.size(x)):
        u[j] = (k/np.sqrt(4*D*r))*(np.e**(-1*np.abs(x[j])*np.sqrt(r/D)))
    return u

uhat = steadystate(x=x,k=100,D=6.45,r=.1)
#k in ug/ML
#D in 10^-6 cm^2/s


#now goal is to look at where conc space dominance and what is dominant where
most_fit_at_conc = np.zeros(np.size(u))
conc_space = u
options = {
           'doubling_time':1.5,
           'death_rate':0.016,
           'mut_rate':10**-4,
           'max_cells':10**6,
           'plot':True,
           'k_abs':0.001,
           'k_elim':10**-4,
           'max_dose':10**1,
           'curve_type':'pharm',
           'timestep_scale':3,
           'n_sims':1,
           'debug':False,
           'fitness_data':'random',
           'n_allele':4,
           }

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


fig, ax = plt.subplots()
ax.plot(x, uhat,color='k')
ax.set(xlabel='x (10^-3 cm)', ylabel='drug concentration (ug/mL)',
       title='Steady State Drug Diffusion w/Constant Source & Reaction')
plt.yscale("log")

#we have list of colors for each thing, where we want each color
#now need to chop up x & uhat into different arrays based on where there is this optimal conc
for l in range(np.size(x)):
    if 0 == np.all(most_fit_at_conc == most_fit_at_conc[l]):
        optimal_at_x = most_fit_at_conc[l]
        optimal_at_x = int(optimal_at_x)
        color_of_optimal = colors[optimal_at_x]
        for c in range(np.size(x)):
            if c+l > (np.size(x)-1):
                break
            else:
                tester = int(most_fit_at_conc[c+l])
                if optimal_at_x != tester:
                    terminal_point = c+l-1
        
        # if optimal_at_x > 8:
        #      ax.fill_between(x,0,uhat,color=color_of_optimal)
        # else:
        #      ax.fill_between(x,0,uhat,color=color_of_optimal,alpha=.5)
        ax.fill_between(x[l:terminal_point],0,uhat[l:terminal_point],color = color_of_optimal)
    else:
         optimal_at_x = most_fit_at_conc[l]
         optimal_at_x = int(optimal_at_x)
         color_of_optimal = colors[optimal_at_x]
         # if optimal_at_x > 8:
         #     ax.fill_between(x,0,uhat,color=color_of_optimal)
         
         ax.fill_between(x,0,uhat,color=color_of_optimal,alpha=.5)

            
#ax.fill_between(x, 0, uhat)
plt.show()


