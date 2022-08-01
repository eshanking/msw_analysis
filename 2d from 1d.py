#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 20 13:24:24 2022

@author: beckettpierce
"""


import numpy as np
import matplotlib.pyplot as plt
from seascapes_figures.classes.population_class import Population
import seascapes_figures
import math


xdim = np.linspace(-100,100,num=200)
ydim = np.linspace(-100,100,num=200)
uhat = np.zeros((np.size(xdim),np.size(xdim)))
def steadystate(x,k,D,r):
    uhat = (k/np.sqrt(4*D*r))*(np.e**(-1*np.abs(x)*np.sqrt(r/D)))
    return uhat

for o in range(np.size(xdim)):
    for q in range(np.size(xdim)):
        uhat[o,q] = steadystate((np.sqrt((xdim[o]**2) + ydim[q]**2)),k=100,D=6.45,r=.1)

#k in ug/ML
#D in 10^-6 cm^2/s


#now goal is to look at where conc space dominance and what is dominant where
most_fit_at_conc =np.zeros((np.size(xdim),np.size(xdim)))
conc_space = uhat
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
for z in range(np.size(xdim)):
    for j in range(np.size(xdim)):
        conc = conc_space[z,j]
        p_fit_list = p.gen_fit_land(conc)
        most_fit_at_conc[z,j] = (np.argmax(p_fit_list))
##now most_fit_at_conc is our list of what is most fit at each point, now we must associate colors with each of these
cc = p.gen_color_cycler()
cc_dict = cc.by_key()
colors = cc_dict['color']


fig, ax = plt.subplots(1,3,figsize=(15, 7))
ax[1].set(xlabel='x (10^-3 cm)', ylabel='y (10^-3 cm)')
# plt.yscale("log")


counts = np.zeros(16)
#we have list of colors for each thing, where we want each color
#now need to chop up x & uhat into different arrays based on where there is this optimal conc
for l in range(np.size(xdim)):
    for w in range(np.size(xdim)):
            optimal_at_x = most_fit_at_conc[l,w]
            optimal_at_x = int(optimal_at_x)
            counts[optimal_at_x] =   counts[optimal_at_x] + 1
            color_of_optimal = colors[optimal_at_x]
            ax[1].scatter(xdim[l],ydim[w], color=colors[optimal_at_x])



list = []
for i in range(16):
    list.append(p.int_to_binary(i))
    



# bar_list=[(f'{binary[0]}',counts[0]),('{binary[1]}',counts[1]),('{binary[2]}',counts[2]),
#           ('{binary[3]}',counts[3]),('{binary[4]}',counts[4]),('{binary[5]}',counts[5]),
#           ('{binary[6]}',counts[6]),('{binary[7]}',counts[7]),('{binary[8]}',counts[8]),
#           ('{binary[9]}',counts[9]),('{binary[10]}',counts[10]),
#           ('{binary[11]}',counts[11],)('{binary[12]}',counts[12]),
#           ('{binary[13]}',counts[13])('{binary[14]}',counts[14]),('{binary[15]}',counts[15])]

# labels, ys = zip(*bar_list)
# xs = np.arange(len(labels)) 
# width = 1
counts = counts / np.sum(counts)

ax[2].bar(list, counts, align='center')

# pad = .75
# pos = ax[2].get_position()
# pos.y1 = pos.y1-pad
# ax[2].set_position(pos)


# pos = ax[2].get_position()
# pos.y0 = pos.y0+pad
# pos.y1 = pos.y1+pad
# ax[2].set_position(pos)

# ax[2].xticks(xs, labels) #Replace default x-ticks with xs, then replace xs with labels
# ax[2].yticks(ys)




ax[0].imshow(uhat,cmap='hot')



# pad = .75
# pos = ax[1].get_position()
# pos.y1 = pos.y1-pad
# ax[1].set_position(pos)


# pos = ax[1].get_position()
# pos.y0 = pos.y0+pad
# pos.y1 = pos.y1+pad
# ax[1].set_position(pos)


# pad = .75
# pos = ax[0].get_position()
# pos.y0 = pos.y0+pad
# pos.y1 = pos.y1+pad
# ax[0].set_position(pos)


ax[0].set_xticks(np.arange(5),('-100','-50','0','50','100'))
ax[0].set_yticks(np.arange(5),('-100','-50','0','50','100'))



 
ax[1].set_xlim((-100, 100))
ax[1].set_ylim((-100, 100))