    #!/usr/bin/env python3
    # -*- coding: utf-8 -*-
    """
    Created on Wed Aug  3 15:08:29 2022
    
    @author: beckettpierce
    """
    
    import numpy as np
    import matplotlib.pyplot as plt
    from seascapes_figures.classes.population_class import Population
    import seascapes_figures
    import math
    import scipy.signal
    import random
    
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
    
    delta_array = np.zeros((200,200))
    delta_array[50,100]=1
    delta_array[150,100]=1
    
    
    r= scipy.signal.convolve2d(uhat,delta_array)
    
    convolxdim = np.linspace(-200,200,num=400)
    convolydim = np.linspace(-200,200,num=400)
    tester = np.log(r)
    rep_val = 0
    tester[np.isinf(tester)] = -13.5
    #plt.imshow(tester, cmap = 'hot')
    
    #now goal is to look at where conc space dominance and what is dominant where
    most_fit_at_conc =np.zeros((np.size(xdim)*2,np.size(xdim)*2))
    conc_space = r
    seed=9345
    np.random.seed(seed)
    random.seed(seed)
    options = {    'k_abs':.95,
        'k_elim':.00839,
        'max_dose':5,
        'n_timestep':4000,
        'timestep_scale':.25,
        'fitness_data':'random',
        'curve_type':'pulsed',
        'prob_drop':.4,
        'n_allele':2
        }
    
    p = Population(**options)
    p.plot_fitness_curves()
    for z in range(np.size(convolydim)-1):
        for j in range(np.size(convolydim)-1):
            conc = conc_space[z,j]
            p_fit_list = p.gen_fit_land(conc)
            most_fit_at_conc[z,j] = (np.argmax(p_fit_list))
    ##now most_fit_at_conc is our list of what is most fit at each point, now we must associate colors with each of these
    cc = p.gen_color_cycler()
    cc_dict = cc.by_key()
    colors = cc_dict['color']
    
    
    fig, ax = plt.subplots(1,3,figsize=(15, 7))
    ax[1].set(xlabel='x (10^-3 cm)')
    # plt.yscale("log")
    
    
    counts = np.zeros(16)
    #we have list of colors for each thing, where we want each color
    #now need to chop up x & uhat into different arrays based on where there is this optimal conc
    
    scalar = 20
    
    ###need to change this
    for l in range(int((np.size(convolydim))/scalar)):
        for w in range(int((np.size(convolxdim))/scalar)):
                optimal_at_x = most_fit_at_conc[w*scalar,l*scalar]
                optimal_at_x = int(optimal_at_x)
                counts[optimal_at_x] =   counts[optimal_at_x] + 1
                color_of_optimal = colors[optimal_at_x]
                ax[1].scatter(convolydim[l*scalar],convolxdim[w*scalar], color=colors[optimal_at_x])
    
    
    
    list = []
    for i in range(16):
        if counts[i] != 0:
            list.append(p.int_to_binary(i))
        else:
            list.append("")
            
    plotting_list = []
    counting_list = []
    list_array = np.array(list)
    for i in range(16):
        if list_array[i] != '':
           plotting_list.append(i)
        if counts[i] != 0:
            counting_list.append(counts[i])
    
    
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
    counting_list = counting_list / np.sum(counting_list)
    
    keys, values = zip(*sorted(zip(plotting_list,counting_list)))
    x = range(len(keys))
    ax[2].bar(x, values, .8, align='center')
    ax[2].set_xticks(x)
    ax[2].set(xlabel='Genotype',ylabel='Proportion of Total This Genotype is Dominant')
#    ax[2].tick_params(bottom=False)
    #ax[2].xticks(x, label='')
    
    
    
    
    pad = .325
    pos = ax[2].get_position()
    pos.y1 = pos.y1-pad
    ax[2].set_position(pos)
    
    pad= .2
    pos = ax[2].get_position()
    pos.y0 = pos.y0+pad
    pos.y1 = pos.y1+pad
    ax[2].set_position(pos)
    
    
    
    
    
    
    # ax[2].xticks(xs, labels) #Replace default x-ticks with xs, then replace xs with labels
    # ax[2].yticks(ys)
    
    
    #locs, labels = ax[2].xticks()
    

    
    ax[0].imshow(tester,cmap='hot',extent = [-200,200,-200,200])
    ax[0].set_ylim((-100, 100))
    ax[0].set_xlim((-100, 100))

    ax[0].set(xlabel='x (10^-3 cm)', ylabel='y (10^-3 cm)')
    
    
    
    pad = .325
    pos = ax[1].get_position()
    pos.y1 = pos.y1-pad
    ax[1].set_position(pos)
    
    
    
    pad = .2
    pos = ax[1].get_position()
    pos.y0 = pos.y0+pad
    pos.y1 = pos.y1+pad
    ax[1].set_position(pos)
    

        
    
    # pad = .05
    # pos = ax[0].get_position()
    # pos.y1 = pos.y1-pad
    # ax[0].set_position(pos)

    
    
    # pad = .04
    # pos = ax[0].get_position()
    # pos.y0 = pos.y0+pad
    # pos.y1 = pos.y1+pad
    # ax[0].set_position(pos)
    
    
    
    
     
    ax[1].set_xlim((-100, 100))
    ax[1].set_ylim((-100, 100))
    
