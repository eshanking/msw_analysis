#%%
"""
Created on Thu Jul 14 13:30:19 2022

@author: beckettpierce
"""

import numpy as np
import matplotlib.pyplot as plt
from fears.population import Population
from fears.utils import plotter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import random

def steadystate(x,k,D,r):

    u = np.zeros(len(x))

    for j in range(len(x)):
        u[j] = (k/np.sqrt(4*D*r))*(np.e**(-1*np.abs(x[j])*np.sqrt(r/D)))
    return u

def most_fit_at_conc(dc,pop):
    
    # for each point in the drug concentration curve, what is the most fit genotype
    mf = np.zeros(len(dc))

    for c in range(len(dc)):

        conc = dc[c]
        p_fit_list = pop.gen_fit_land(conc)
        mf[c] = (np.argmax(p_fit_list))
    
    return mf

def detect_changes(mf):
    
    x = 0
    chunks = []

    last_most_fit = mf[0]
    new_chunk = False
    left_indx = 0

    while x < len(mf):
        cur_most_fit = mf[x]
        if last_most_fit != cur_most_fit:

            right_indx = x
            chunks.append((left_indx,right_indx))
            left_indx = x

        x+=1
        last_most_fit = cur_most_fit
    
    chunks.append((right_indx,len(mf)))
    
    return chunks

def plot_drug_curve(ax,dc,mf,chunks,colors,pop,**kwargs):
    
    n = 0
    for chunk in chunks:
        x = np.arange(chunk[0],chunk[1])
        dc_t = dc[chunk[0]:chunk[1]]

        most_fit = mf[chunk[0]]
        color = colors[int(most_fit)]
        l = pop.int_to_binary(int(most_fit))
        ax.plot(x,dc_t,color=color,label=l,**kwargs)
        n+=1
    
    return ax

#%%
seed = 2022
np.random.seed(seed)
random.seed(seed)

fig, ax = plt.subplots(1,2,figsize=(8, 3))

options = {'k_abs':.95,
    'k_elim':.00839,
    'max_dose':5,
    'n_timestep':10000,
    'timestep_scale':0.05,
    'fitness_data':'random',
    'curve_type':'pulsed',
    'prob_drop':.5,
    'n_allele':2
    }

drug_conc_range = [-4,4]
p = Population(
               death_rate=0.1,
               drug_conc_range = drug_conc_range,
               ic50_limits=[-2.5,3],
               drugless_limits=[0.8,1.5],
               **options)

p = Population(**options)
p.drugless_rates = [1.28949852, 1.14399848, 1.22802236, 0.93619847]
p.ic50 = [-0.49205992, 1.76224515,  1.39341393,  2.84653598]

#need number/ to find the conc space
u = p.drug_curve

mf = most_fit_at_conc(u,p)
chunks = detect_changes(mf)
cc = plotter.gen_color_cycler()
cc_dict = cc.by_key()
colors = cc_dict['color']

ax[0] = plot_drug_curve(ax[0],u,mf,chunks,colors,p,linewidth=2)

ax[0] = plotter.x_ticks_to_days(p,ax[0])

ax[0].set_xlabel('Time (days)',fontsize=15)
ax[0].set_ylabel('Drug concentration (ug/mL)',fontsize=15)
# ax[0].set_yscale('log')

#%%

x = np.linspace(0,100,num=100)

diff = steadystate(x=x,k=100,D=6.45,r=.1)
#k in ug/ML
#D in 10^-6 cm^2/s

#now goal is to look at where conc space dominance and what is dominant where

mf = most_fit_at_conc(diff,p)
chunks = detect_changes(mf)
cc = plotter.gen_color_cycler()
cc_dict = cc.by_key()
colors = cc_dict['color']

ax[1] = plot_drug_curve(ax[1],diff,mf,chunks,colors,p,linewidth=5)

# ax[1].set(xlabel='x ($10^{-3}$ cm)', ylabel='Drug Concentration (ug/ml)')
ax[1].set_xlabel(xlabel='x ($10^{-3}$ cm)',fontsize=15)
ax[1].set_ylabel(ylabel='Drug Concentration (ug/mL)',fontsize=15)
ax[1].set_yscale('log')

ax[1] = plotter.shiftx(ax[1],0.1)

for a in ax:
    a.tick_params(axis='both', which='major', labelsize=12)

handles, labels = ax[1].get_legend_handles_labels()

unique_labels = sorted(set(labels))
labels = np.array(labels)
unique_handles = []

for lab in unique_labels:
    indx = np.argwhere(labels==lab)
    indx = indx[0][0]
    unique_handles.append(handles[indx])

ax[1].legend(unique_handles,unique_labels,loc = (-1,-0.4),frameon=False,
             fontsize=12,ncol=4)

fig.savefig('figures/1d_time_and_diff.pdf',bbox_inches='tight')