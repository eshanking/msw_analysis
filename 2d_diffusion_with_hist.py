#%%
import numpy as np
import matplotlib.pyplot as plt
from fears.population import Population
from fears.utils import plotter
import scipy
# from matplotlib.colors import ListedColormap
from matplotlib import colors

def indx_matrix_from_zero(indx,mat):

    ycenter = int(mat.shape[0]/2)
    xcenter = int(mat.shape[1]/2)

    y = indx[0] - ycenter
    x = indx[1] - xcenter
    val = mat[int(y),int(x)]

    return val

def twoD_eqn(umax,x,D_vessel=0.01):

    if x < D_vessel:
        u = umax
    else:
        u = umax*scipy.special.kn(0,x)
    return u

#%% Initialize population

seed = 2022
np.random.seed(seed)

# fig, ax = plt.subplots(1,2,figsize=(8, 3))

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

ic50_rank = np.sort(p.ic50)
umax = 20

#%%
s = 2 # mm, i.e. a 2 mm patch of tissue
step = 0.01 # each pixel is 10 um x 10 um
D_vessel = 0.01 # vessel diameter is 10 um

xdim = np.arange(-s,s,step)
ydim = np.arange(-s,s,step)
uhat = np.zeros((np.size(xdim),np.size(xdim)))

u_distances = np.zeros((np.size(xdim),np.size(xdim))) # a matrix of the distance of each pixel from the origin

# Compute steady-state solution in radial coordinates
for i in range(np.size(xdim)):
    for j in range(np.size(xdim)):
        gamma = np.sqrt((xdim[i]**2) + ydim[j]**2)
        uhat[i,j] = twoD_eqn(umax, gamma,D_vessel=D_vessel)
        u_distances[i,j] = gamma
        # uhat[i,j] = steadystate((np.sqrt((xdim[i]**2) + ydim[j]**2)),k=100,D=6.45,r=.1)

uhat = umax*uhat/np.max(uhat) # normalize to umax

uhat[u_distances < D_vessel] = umax # make sure the interior of the vessel is umax
#%%

s_d = s

delta_array = np.zeros((np.size(xdim),np.size(xdim)))

delta_array[200,160] = 1
delta_array[200,240] = 1

# Convolve the steady state solution with the impulse matrix
conv = scipy.signal.convolve2d(delta_array,uhat,mode='same')#,mode='same',boundary='wrap')

conv = conv[1:,:] # weird boundary conditions
log10_conv = np.log10(conv)
#%%
fig,ax = plt.subplots()

im = ax.imshow(log10_conv,extent=[-s_d,s_d,-s_d,s_d])
# im = ax.imshow(log10_conv)
fig.colorbar(im)


#%%
#now goal is to look at where conc space dominance and what is dominant where

most_fit_at_conc = np.zeros(log10_conv.shape)

for z in range(log10_conv.shape[0]):
    for j in range(log10_conv.shape[1]):
        conc = 10**log10_conv[z,j]
        p_fit_list = p.gen_fit_land(conc)
        most_fit_at_conc[z,j] = int((np.argmax(p_fit_list)))
        # most_fit.append(int((np.argmax(p_fit_list))))
        # conc_vect.append(conc)
        

##now most_fit_at_conc is our list of what is most fit at each point, now we must associate colors with each of these
cc = plotter.gen_color_cycler(style='solid',n_colors=4,palette='colorblind')
cc_dict = cc.by_key()
c = cc_dict['color']

indx = list(set(most_fit_at_conc.flatten()))
indx = [int(i) for i in indx]
c = [c[i] for i in indx]

cmap = colors.ListedColormap(c)
bounds=indx
norm = colors.BoundaryNorm(bounds, cmap.N)
#%% plotting

fig, ax = plt.subplots(1,3,figsize=(15, 7))

counts = np.zeros(p.n_genotype)
#we have list of colors for each thing, where we want each color
#now need to chop up x & uhat into different arrays based on where there is this optimal conc


final_range_x = (-150,150)
final_range_y = (-200,200)

# Get counts only in the range we are plotting
for l in np.arange(final_range_x[0],final_range_x[1]):#range(final_range[1] - final_range[0]):
    for w in np.arange(final_range_y[0],final_range_y[1]):
            optimal_at_x = int(indx_matrix_from_zero((w,l),most_fit_at_conc))
            counts[optimal_at_x] =   counts[optimal_at_x] + 1

imfit = ax[1].imshow(most_fit_at_conc,cmap=cmap,extent = [-s_d,s_d,-s_d,s_d],
                     interpolation='gaussian',interpolation_stage='rgba')

# imfit = ax[1].imshow(most_fit_at_conc,cmap=cmap,extent = [-200,200,-200,200])

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
im = ax[0].imshow(log10_conv,cmap='inferno',extent = [-s_d,s_d,-s_d,s_d])
# ax[0].set_ylim((-100, 100))
# ax[0].set_xlim((-100, 100))

# ax[1].set_ylim((-100, 100))
# ax[1].set_xlim((-100, 100))

ax[0].set_ylabel('y (mm)',fontsize=15)
ax[0].set_xlabel('x (mm)',fontsize=15)
ax[1].set_ylabel('y (mm)',fontsize=15)
ax[1].set_xlabel('x (mm)',fontsize=15)

for a in ax:
    a.tick_params(axis='both', which='major', labelsize=12)

# Display only this range
# ax[1].set_xlim((-100, 100))
# ax[1].set_ylim((-100, 100))

for a in ax[0:2]:
    a.set_aspect('equal')

pad = 0.03

# Increase x spacing of axes

cb = fig.colorbar(im, ax=ax[0],location='bottom')
cb.ax.tick_params(labelsize=12) 
cb.set_label('Drug concentration (log($\mu$g/mL))',fontsize=12)

cbfit = fig.colorbar(imfit,ax=ax[1],cmap=cmap,location='bottom',
    ticks = [1.3,2,2.7])

cbfit.ax.set_xticklabels(['01','10','11'])
cbfit.ax.tick_params(labelsize=12) 

cbfit.set_label('Genotype',fontsize=12)

ax[0].set_xlim(-1.5,1.5)
ax[1].set_xlim(-1.5,1.5)

pos0 = ax[0].get_position()

for a in ax[1:]:
    pos = a.get_position()

    pos.y0 = pos0.y0
    pos.y1 = pos0.y1

    a.set_position(pos)

# ax[1] = plotter.shiftx(ax[1],pad)
# ax[2] = plotter.shiftx(ax[2],pad*2)
# cbfit.ax = plotter.shiftx(cbfit.ax,pad)

fig.savefig('figures/2d_diffusion_with_hist.pdf',bbox_inches='tight')
# %%
