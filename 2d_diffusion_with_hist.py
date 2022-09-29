"""
Created on Wed Aug  3 15:08:29 2022

@author: beckettpierce
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
from fears.population import Population
from fears.utils import plotter
import math
import scipy.signal
import random
# from matplotlib.colors import ListedColormap
from matplotlib import colors

def indx_matrix_from_zero(indx,mat):

    ycenter = int(mat.shape[0]/2)
    xcenter = int(mat.shape[1]/2)

    y = indx[0] - ycenter
    x = indx[1] - xcenter
    val = mat[int(y),int(x)]

    return val

def steadystate(x,k,D,r):
    uhat = (k/np.sqrt(4*D*r))*(np.e**(-1*np.abs(x)*np.sqrt(r/D)))
    return uhat

#%%
xdim = np.linspace(-100,100,num=200)
ydim = np.linspace(-100,100,num=200)
uhat = np.zeros((np.size(xdim),np.size(xdim)))

# Compute steady-state solution in radial coordinates
for i in range(np.size(xdim)):
    for j in range(np.size(xdim)):
        uhat[i,j] = steadystate((np.sqrt((xdim[i]**2) + ydim[j]**2)),k=100,D=6.45,r=.1)

#k in ug/ML
#D in 10^-6 cm^2/s

# Define impulse matrix (location of vessels)
delta_array = np.zeros((200,200))
delta_array[50,100]=1
delta_array[150,100]=1

# Convolve the steady state solution with the impulse matrix
r = scipy.signal.convolve2d(uhat,delta_array)


convolxdim = np.linspace(-200,200,num=400)
convolydim = np.linspace(-200,200,num=400)

tester = np.log(r)
rep_val = 0
tester[np.isinf(tester)] = -13.5
#plt.imshow(tester, cmap = 'hot')

#now goal is to look at where conc space dominance and what is dominant where
most_fit_at_conc = np.zeros((np.size(xdim)*2,np.size(xdim)*2))
conc_space = r
# seed=9345
#%%
seed = 109
np.random.seed(seed)
random.seed(seed)

drug_conc_range = [-4,4]
p = Population(fitness_data='random',
               n_allele=2,
               death_rate=0.1,
               drug_conc_range = drug_conc_range,
               ic50_limits=[-2.5,3],
               drugless_limits=[0.8,1.5])

p.drugless_rates = [1.28949852, 1.14399848, 1.22802236, 0.93619847]
p.ic50 = [-0.49205992, 1.76224515,  1.39341393,  2.84653598]
p.plot_fitness_curves()
#%%
for z in range(np.size(convolydim)-1):
    for j in range(np.size(convolydim)-1):
        conc = conc_space[z,j]
        p_fit_list = p.gen_fit_land(conc)
        most_fit_at_conc[z,j] = int((np.argmax(p_fit_list)))

##now most_fit_at_conc is our list of what is most fit at each point, now we must associate colors with each of these
cc = plotter.gen_color_cycler()
cc_dict = cc.by_key()
c = cc_dict['color']

c = c[:int(np.max(most_fit_at_conc))+1]

cmap = colors.ListedColormap(c)
bounds=[0,1,2,3,4]
norm = colors.BoundaryNorm(bounds, cmap.N)
#%% plotting

fig, ax = plt.subplots(1,3,figsize=(15, 7))

counts = np.zeros(p.n_genotype)
#we have list of colors for each thing, where we want each color
#now need to chop up x & uhat into different arrays based on where there is this optimal conc

final_range = (-100,100)

# Get counts only in the range we are plotting
for l in np.arange(final_range[0],final_range[1]):#range(final_range[1] - final_range[0]):
    for w in np.arange(final_range[0],final_range[1]):
            optimal_at_x = int(indx_matrix_from_zero((w,l),most_fit_at_conc))
            counts[optimal_at_x] =   counts[optimal_at_x] + 1


# resolution_scale = 10
# for l in range(int((np.size(convolydim))/resolution_scale)):
#     for w in range(int((np.size(convolxdim))/resolution_scale)):
#             optimal_at_x = most_fit_at_conc[w*resolution_scale,l*resolution_scale]
#             optimal_at_x = int(optimal_at_x)
#             # counts[optimal_at_x] =   counts[optimal_at_x] + 1
#             label = p.int_to_binary(optimal_at_x)
#             color_of_optimal = colors[optimal_at_x]
#             ax[1].scatter(convolydim[l*resolution_scale],
#                           convolxdim[w*resolution_scale], 
#                           color=colors[optimal_at_x],
#                           label=label)

imfit = ax[1].imshow(most_fit_at_conc,cmap=cmap,extent = [-200,200,-200,200],
                     interpolation='gaussian',interpolation_stage='rgba',norm=norm)

counts = counts / np.sum(counts)

x = np.arange(p.n_genotype)
xlab = [p.int_to_binary(xx) for xx in x]
ax[2].bar(x, counts, .8, align='center',color='coral')
ax[2].set_xticks(x)
ax[2].set_xticklabels(xlab)
# ax[2].set(xlabel='Genotype',ylabel='Proportion of Mutant Selection Windows')
ax[2].set_xlabel('Genotype',fontsize=15)
ax[2].set_ylabel('Proportion',fontsize=15)

# Plot 2D diffision
im = ax[0].imshow(tester,cmap='hot',extent = [-200,200,-200,200],vmin=-7,vmax=3)
ax[0].set_ylim((-100, 100))
ax[0].set_xlim((-100, 100))

ax[0].set_ylabel('y ($10^{-3}$ cm)',fontsize=15)
ax[0].set_xlabel('x ($10^{-3}$ cm)',fontsize=15)
ax[1].set_ylabel('y ($10^{-3}$ cm)',fontsize=15)
ax[1].set_xlabel('x ($10^{-3}$ cm)',fontsize=15)

for a in ax:
    a.tick_params(axis='both', which='major', labelsize=12)

# Display only this range
ax[1].set_xlim((-100, 100))
ax[1].set_ylim((-100, 100))

for a in ax[0:2]:
    a.set_aspect('equal')

pad = 0.03

# Increase x spacing of axes

cb = fig.colorbar(im, ax=ax[0],location='bottom')
cb.ax.tick_params(labelsize=12) 
cb.set_label('Drug concentration (log(ug/mL))',fontsize=12)

cbfit = fig.colorbar(imfit,ax=ax[1],cmap=cmap,location='bottom',
    ticks = [0.5,1.5,2.5,3.5])

cbfit.ax.set_xticklabels(['00','01','10','11'])
cbfit.ax.tick_params(labelsize=12) 

pos0 = ax[0].get_position()

for a in ax[1:]:
    pos = a.get_position()

    pos.y0 = pos0.y0
    pos.y1 = pos0.y1

    a.set_position(pos)

ax[1] = plotter.shiftx(ax[1],pad)
ax[2] = plotter.shiftx(ax[2],pad*2)
cbfit.ax = plotter.shiftx(cbfit.ax,pad)
# Add legend to middle axes
# handles, labels = ax[1].get_legend_handles_labels()

# unique_labels = sorted(set(labels))
# labels = np.array(labels)
# unique_handles = []

# for lab in unique_labels:
#     indx = np.argwhere(labels==lab)
#     indx = indx[0][0]
#     unique_handles.append(handles[indx])

# ax[1].legend(unique_handles,unique_labels,loc = (0,-0.3),frameon=True,
#              fontsize=12,ncol=4,fancybox=False,framealpha=1,
#              columnspacing=1)
# ax[1].colorbar()

fig.savefig('figures/2d_diffusion_with_hist.pdf',bbox_inches='tight')
# %%
