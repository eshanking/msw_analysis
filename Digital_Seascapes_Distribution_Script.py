#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 11:42:33 2022

@author: beckettpierce
"""




from seascapes_figures.classes.population_class import Population
import numpy as np
import matplotlib.pyplot as plt

#k,r = p.simulate()
#initial fitness curve

#NEW GOAL -> arvitratily try it a t 10^-2, 10^-1, 1, 10,100,1000,10000
conc_space = np.array([10**-4,10**-3,10**-2,10**-1,10**0,10**1,10**2,10**3,10**4,10**5,10**6]) 
counts = np.zeros(np.size(conc_space))

for j in conc_space:
    conc = conc_space[np.where(conc_space==j)]
    for k in range(500):
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
        #ic50,initial rates, and
        true_ic50 = (10**((p.ic50)))
        # ic50_rates = p.drugless_rates/2 
    
        ic50 = p.ic50
        g_drugless = p.drugless_rates
    
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
                   #'fitness_data':'random',
                   'n_allele':4,
                   }
        p1 = Population(ic50=ic50,drugless_rates=g_drugless,digital_seascape=True,**options)
        p1.drugless_rates = p1.drugless_rates/max(p1.drugless_rates)
    
        p_fit_list = p.gen_fit_land(conc)
        p1_fit_list = p1.gen_fit_land(conc)
        np.where((p_fit_list) == max(p.gen_fit_land(conc)))
        np.where((p1_fit_list) == max(p1.gen_fit_land(conc)))
        if np.where((p_fit_list) == max(p.gen_fit_land(conc))) == np.where((p1_fit_list) == max(p1.gen_fit_land(conc))):
            counts[np.where(conc_space==j)] = counts[np.where(conc_space==j)]+1
## just need to actually plot these things now on a distribution
percent_counts = counts / 500
fig, ax = plt.subplots()
plt.xscale("log")
ax.plot(conc_space, percent_counts,color='r')
ax.set(xlabel='Concentration uM', ylabel='Accuracy Rate',
       title='Digital Seascape Accuracy at Different Drug Concentrations')

##

