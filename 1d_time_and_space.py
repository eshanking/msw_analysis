#%%
import numpy as np
import matplotlib.pyplot as plt
from fears.population import Population
from fears.utils import plotter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import random

# def steadystate(x,k,D,r):

#     u = np.zeros(len(x))

#     for j in range(len(x)):
#         u[j] = (k/np.sqrt(4*D*r))*(np.e**(-1*np.abs(x[j])*np.sqrt(r/D)))
#     return u

def oneD_eqn(umax,x):
    return umax*np.exp(-x)

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

def plot_drug_curve(ax,dc,mf,chunks,colors,pop,x=None,**kwargs):
    
    n = 0
    for chunk in chunks:
        if x is None:
            x_t = np.arange(chunk[0],chunk[1])
        else:
            x_t = x[chunk[0]:chunk[1]]
        dc_t = dc[chunk[0]:chunk[1]]

        most_fit = mf[chunk[0]]
        color = colors[int(most_fit)]
        l = pop.int_to_binary(int(most_fit))
        ax.plot(x_t,dc_t,color=color,label=l,**kwargs)
        n+=1
    
    return ax

def plot_msw_vspan(ax,mf,chunks,colors,pop,x=None,**kwargs):
    
    n = 0
    for chunk in chunks:
        if x is None:
            x_t = np.arange(chunk[0],chunk[1])
        else:
            x_t = x[chunk[0]:chunk[1]]

        most_fit = mf[chunk[0]]
        color = colors[int(most_fit)]
        l = pop.int_to_binary(int(most_fit))
        ax.axvspan(x_t[0],x_t[-1],color=color,label=l,**kwargs)
        n+=1
    
    return ax

#%%
seed = 2022
np.random.seed(seed)
random.seed(seed)

fig, ax = plt.subplots(1,2,figsize=(8, 3))

options = {'death_model':None,
           'k_abs':.95,
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
cc = plotter.gen_color_cycler(style='solid',n_colors=4,palette='colorblind')
cc_dict = cc.by_key()
colors = cc_dict['color']

ax[0] = plot_drug_curve(ax[0],u,mf,chunks,colors,p,linewidth=4)

ax[0] = plotter.x_ticks_to_days(p,ax[0])

ax[0].set_xlabel('Time (days)',fontsize=15)
ax[0].set_ylabel('Drug concentration ($\mu$g/mL)',fontsize=15)
# ax[0].set_yscale('log')

#%%

x = np.linspace(0,15,num=1000)

ic50_ranked = np.sort(p.ic50)
umax = 10**np.mean((ic50_ranked[-1],ic50_ranked[-2]))
dc = oneD_eqn(umax,x)
# diff = steadystate(x=x,k=100,D=6.45,r=.1)
#k in ug/ML
#D in 10^-6 cm^2/s

#now goal is to look at where conc space dominance and what is dominant where

mf = most_fit_at_conc(dc,p)
chunks = detect_changes(mf)

ax[1] = plot_msw_vspan(ax[1],mf,chunks,colors,p,x=x,linewidth=5)
ax[1].plot(x,dc,linewidth=6,color='black')

# ax[1].set(xlabel='x ($10^{-3}$ cm)', ylabel='Drug Concentration (ug/ml)')
ax[1].set_xlabel(xlabel='x (mm)',fontsize=15)
ax[1].set_ylabel(ylabel='Drug Concentration ($\mu$g/mL)',fontsize=15)
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

ax[1].legend(unique_handles,unique_labels,loc = (-1,-0.45),frameon=False,
             fontsize=12,ncol=4)

fig.savefig('figures/1d_time_and_diff.png',bbox_inches='tight')
# %%
plotter.plot_msw(p,0,figsize=(7,3))
# %%

def log_eqn(g,ic50,conc):
    hc = -.5
    res = []
    for c in conc:
        res.append(g/(1+np.exp((ic50-np.log10(c))/hc)))
    return res

conc = np.logspace(-3,3)

ic50_1 = -1.5
ic50_2 = 1.5

g1 = 1.1
g2 = 0.8

f1 = log_eqn(g1,ic50_1,conc)
f2 = log_eqn(g2,ic50_2,conc)

fig,ax = plt.subplots(figsize=(5,3.5))

ax.plot(conc,f1,color='black',linewidth=3,label='sensitive')
ax.plot(conc,f2,color='red',linewidth=3,label='resistant')
ax.set_xscale('log')
ax.tick_params(axis='both', which='major', labelsize=12)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_xlabel('Drug concentration (ug/mL)',fontsize=14)
ax.set_ylabel('Growth rate',fontsize=14)
ax.legend(fontsize=11,frameon=False)

fig.savefig('ex_msw.pdf',bbox_inches='tight')
# %%
