#%%
from fears.population import Population
from fears.utils import plotter, fitness
import numpy as np
import matplotlib.pyplot as plt

np.random.seed(109)
drug_conc_range = [-3,5]
p = Population(fitness_data='random',
               n_allele=2,
               death_rate=0.1,
               drug_conc_range = drug_conc_range,
               ic50_limits=[-2.5,3],
               drugless_limits=[0.8,1.5])

p.drugless_rates = [1.28949852, 1.14399848, 1.22802236, 0.93619847]
p.ic50 = [-0.49205992, 1.76224515,  1.39341393,  2.84653598]
#%%
fig,ax = plt.subplots() 

genotypes = np.array([0,1,2,3])
ax = plotter.msw_grid(p,genotypes,ax=ax,
                      labelsize=12,
                      ticklabelsize=12,
                      comp_annotate_pos=10**5.1,
                      legendloc=(0,0.95))

ax.set_xlabel('Drug concentration ($\mathrm{\mu}$g/mL)',fontsize=14)
fig.savefig('figures/2_allele_grid.pdf')
#%%
# fig2 = plotter.plot_msw(p,2,ncols=1,figsize=(2.5,4))
fig2, ax = plotter.plot_msw(p,0,ncols=2,figsize=(6,2))

# for a in ax:
#     a.yaxis.tick_right()
ax[0].set_ylabel('Replication rate ($hr^{-1}$)',fontsize=14)
ax[0].yaxis.set_label_position("left")
ax[1].set_ylabel('')
fig2.savefig('figures/2_allele_msws.pdf',bbox_inches='tight')

fig3,ax = plt.subplots(figsize=(5,3))

fig3,ax = plotter.plot_fitness_curves(p,ax=ax,fig=fig3,
                                      color_kwargs={'style':'solid',
                                                    'n_colors':4,
                                                    'palette':'colorblind'})
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.legend(loc='best',frameon=False,fontsize=14)
fig3.savefig('figures/random_seascape.pdf')
# %%
