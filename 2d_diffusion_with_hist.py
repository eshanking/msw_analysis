#%%
import numpy as np
import matplotlib.pyplot as plt
from fears.population import Population
from fears.utils import plotter
import scipy
from matplotlib import colors
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

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

# delta_array[200,160] = 1
# delta_array[200,240] = 1
delta_array[160,200] = 1
delta_array[240,200] = 1

# Convolve the steady state solution with the impulse matrix
conv = scipy.signal.convolve2d(delta_array,uhat,mode='same')#,mode='same',boundary='wrap')

conv = conv[:,1:] # weird boundary conditions
log10_conv = np.log10(conv)
# #%%
# fig,ax = plt.subplots()

# im = ax.imshow(log10_conv,extent=[-s_d,s_d,-s_d,s_d])
# # im = ax.imshow(log10_conv)
# fig.colorbar(im)


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

fig, ax_list = plt.subplots(2,2,figsize=(8, 6))

counts = np.zeros(p.n_genotype)
#we have list of colors for each thing, where we want each color
#now need to chop up x & uhat into different arrays based on where there is this optimal conc


# final_range_x = (-150,150)
# final_range_y = (-200,200)
final_range_y = (-150,150)
final_range_x = (-200,200)

# Get counts only in the range we are plotting
# for l in np.arange(final_range_x[0],final_range_x[1]):#range(final_range[1] - final_range[0]):
#     for w in np.arange(final_range_y[0],final_range_y[1]):
#             optimal_at_x = int(indx_matrix_from_zero((w,l),most_fit_at_conc))
#             counts[optimal_at_x] =   counts[optimal_at_x] + 1

ax = ax_list[1,:]

imfit = ax[1].imshow(most_fit_at_conc,cmap=cmap,extent = [-s_d,s_d,-s_d,s_d],
                     interpolation='gaussian',interpolation_stage='rgba')

# Plot 2D diffision
im = ax[0].imshow(log10_conv,cmap='inferno',extent = [-s_d,s_d,-s_d,s_d])

# imfit = ax[1].imshow(most_fit_at_conc,cmap=cmap,extent = [-200,200,-200,200])

# counts = counts / np.sum(counts)

# x = np.arange(p.n_genotype)
# xlab = [p.int_to_binary(xx) for xx in x]
# ax[2].bar(x, counts, .8, align='center',color='coral')
# ax[2].set_xticks(x)
# ax[2].set_xticklabels(xlab)
# # ax[2].set(xlabel='Genotype',ylabel='Proportion of Mutant Selection Windows')
# ax[2].set_xlabel('Genotype',fontsize=15)
# ax[2].set_ylabel('Proportion',fontsize=15)


# ax[0].set_ylim((-100, 100))
# ax[0].set_xlim((-100, 100))

# ax[1].set_ylim((-100, 100))
# ax[1].set_xlim((-100, 100))

ax[0].set_ylabel('y (mm)',fontsize=15)
ax[0].set_xlabel('x (mm)',fontsize=15)
ax[1].set_ylabel('y (mm)',fontsize=15)
ax[1].set_xlabel('x (mm)',fontsize=15)

for a in ax:
    a.tick_params(axis='both', which='major', labelsize=14)

# Display only this range
# ax[1].set_xlim((-100, 100))
# ax[1].set_ylim((-100, 100))

for a in ax:
    a.set_aspect('equal')


# Increase x spacing of axes

cb = fig.colorbar(im, ax=ax_list[1,0],location='bottom',pad=0.2)
cb.ax.tick_params(labelsize=12) 
cb.set_label('Drug concentration (log($\mu$g/mL))',fontsize=14)

# cbfit = fig.colorbar(imfit,ax=ax[1],cmap=cmap,location='bottom',
#     ticks = [1.3,2,2.7],pad=0.1)

# cbfit.ax.set_xticklabels(['01','10','11'])
# cbfit.ax.tick_params(labelsize=12) 

# cbfit.set_label('Genotype',fontsize=14)

ax[0].set_ylim(-1.5,1.5)
ax[1].set_ylim(-1.5,1.5)

pos0 = ax[0].get_position()

# for a in ax[1:]:
#     pos = a.get_position()

#     pos.y0 = pos0.y0
#     pos.y1 = pos0.y1

#     pos.x0 = pos.x0 + 0.06
#     pos.x1 = pos.x1 + 0.06
#     a.set_position(pos)

#     pos = cbfit.ax.get_position()

#     pos.x0 = pos.x0 + 0.06
#     pos.x1 = pos.x1 + 0.06
#     cbfit.ax.set_position(pos)

# ax[1] = plotter.shiftx(ax[1],pad)
# ax[2] = plotter.shiftx(ax[2],pad*2)
# cbfit.ax = plotter.shiftx(cbfit.ax,pad)

# fig.savefig('figures/2d_diffusion_with_hist.pdf',bbox_inches='tight')
#%% 1D time and space

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

#%%
# seed = 2022
# np.random.seed(seed)
# random.seed(seed)

# fig, ax = plt.subplots(1,2,figsize=(8, 3))

# options = {'k_abs':.95,
#     'k_elim':.00839,
#     'max_dose':5,
#     'n_timestep':10000,
#     'timestep_scale':0.05,
#     'fitness_data':'random',
#     'curve_type':'pulsed',
#     'prob_drop':.5,
#     'n_allele':2
#     }

# drug_conc_range = [-4,4]
# p = Population(
#                death_rate=0.1,
#                drug_conc_range = drug_conc_range,
#                ic50_limits=[-2.5,3],
#                drugless_limits=[0.8,1.5],
#                **options)

# p = Population(**options)
# p.drugless_rates = [1.28949852, 1.14399848, 1.22802236, 0.93619847]
# p.ic50 = [-0.49205992, 1.76224515,  1.39341393,  2.84653598]

#need number/ to find the conc space

ax = ax_list[0,:]

u = p.drug_curve

mf = most_fit_at_conc(u,p)
chunks = detect_changes(mf)
cc = plotter.gen_color_cycler(style='solid',n_colors=4,palette='colorblind')
cc_dict = cc.by_key()
colors = cc_dict['color']

ax[0] = plot_drug_curve(ax[0],u,mf,chunks,colors,p,linewidth=2)

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

ax[1] = plot_drug_curve(ax[1],dc,mf,chunks,colors,p,x=x,linewidth=5)

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

ax[1].legend(unique_handles,unique_labels,loc = (-1,-0.5),frameon=False,
             fontsize=12,ncol=4)

# fig.savefig('figures/1d_time_and_diff.pdf',bbox_inches='tight')
# %% Adjust axes

pos0 = ax_list[0,0].get_position()
pos = ax_list[0,1].get_position()
pos.x0 = pos0.x0
pos.x1 = pos0.x1
ax_list[1,0].set_position(pos)

pos1 = ax_list[0,1].get_position()
pos = ax_list[1,1].get_position()
pos.x0 = pos1.x0
pos.x1 = pos1.x1
ax_list[1,1].set_position(pos)

#shift lower row down

pos = ax_list[1,1].get_position()
pos.y0 = pos.y0 - 0.15
pos.y1 = pos.y1 - 0.15

ax_list[1,1].set_position(pos)

pos_t = ax_list[1,0].get_position()
pos_t.y0 = pos.y0
pos_t.y1 = pos.y1
ax_list[1,0].set_position(pos_t)

# adjust colorbars

cb_pos = cb.ax.get_position()
h = cb_pos.height
cb_pos.y0 = pos_t.y0 - 0.15
cb_pos.y1 = cb_pos.y0 + h

cb.ax.set_position(cb_pos)

labels = ['A','B','C','D']
indx = 0

al = np.reshape(ax_list,4)

for ax in al:

    ax.annotate(labels[indx],(0,1.05),xycoords='axes fraction',fontsize=15)
    indx += 1


fig.savefig('figures/1d_2d_combined.pdf',bbox_inches='tight')

# %%
