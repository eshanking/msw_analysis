#%%
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy
import os
import scipy.optimize as sciopt
import seaborn as sns
import matplotlib.patches as mpl_patches

class CancerGame:
    def __init__(self,folder_path,
                 dc=None,
                 cmap=None,
                 num_replicates=3,
                 proportions=None,
                 slope_est_options=None,
                 green_key=None,
                 red_key=None,
                 smooth_counts=False,
                 rolling_average_window=5,
                 fitness_estimate='slope',
                 exclude_wells=[],
                 parental_linear_fit=False,
                 mutant_linear_fit=False,
                 single_plate=False,
                 dt=4,
                 growth_rate_type='logistic'):
        
        if green_key is None:
            self.green_key = 'Count_gfp_objects'
        else:
            self.green_key = green_key

        if red_key is None:
            self.red_key = 'Count_rfp_objects' 
        else:
            self.red_key = red_key

        self.smooth_counts = smooth_counts
        self.rolling_average_window = rolling_average_window
        self.folder_path = folder_path
        if single_plate is False:
            self.plate_paths = self.get_plate_paths()
        else:
            self.plate_paths = [folder_path]

        if cmap is None:
            self.cmap = mpl.colormaps['viridis']
        else:
            self.cmap = cmap

        if dc is None:
            dc = [1]
            for i in range(1,9):
                dc.append(dc[i-1]/(2))

            dc.reverse()
            self.dc = dc
        else:
            self.dc = dc

        if proportions is None:

            self.proportions = [1,0.9,0.7,0.5,0.3,0.1,0.05,0]
        else:
            self.proportions = proportions

        self.data_list = []
        for plate_path in self.plate_paths:
            self.data_list.append(self.get_data(plate_path))

        if slope_est_options is None:
            self.slope_est_options = {'exclude':3,'window':5,'step':2,'thresh':0.5}
        else:
            self.slope_est_options = slope_est_options
        
        self.fitness_resistant = None
        self.fitness_parental = None
        self.resistant_dr_params = None
        self.parental_dr_params = None
        self.fitness_estimate = fitness_estimate
        self.num_replicates = num_replicates
        self.exclude_wells = exclude_wells # tuple of (plate,well) to exclude from analysis
        # self.dt = dt # per hour
        self.mutant_linear_fit = mutant_linear_fit
        self.parental_linear_fit = parental_linear_fit
        self.dt = dt # sampling frequency in hours
        self.growth_rate_type = growth_rate_type # 'logistic' or 'linear'

    def execute(self):
        # self.fitness_resistant,self.fitness_parental = self.get_fitness_lib()
        self.set_fitness_params()
        # self.plot_dose_response_curves()
        self.set_dose_response_params()

    def set_fitness_params(self):
        self.fitness_resistant,self.fitness_parental = self.get_fitness_lib()

    def get_plate_paths(self,folder_path=None):
        """Gets plate data paths
        Returns:
            list: list of plate data paths
        """
        if folder_path is None:
            folder_path = self.folder_path

        plate_files = os.listdir(path=folder_path)

        #Need to make sure we are only attempting to load .csv or .xlsx data
        plate_files = [i for i in plate_files if i[-4:] == '.csv' or i[-5:] == '.xlsx' and i[0] != '~']

        plate_files.sort()

        plate_data_paths = []

        for pf in plate_files:
            if pf[0] != '.':
                plate_path = folder_path + os.sep + pf
                plate_data_paths.append(plate_path)

        # sorts plates based on file name in a way a human would
        plate_data_paths.sort(key=self.natural_keys)
        return plate_data_paths

    def atoi(self,text):
        return int(text) if text.isdigit() else text

    def natural_keys(self,text):
        '''
        alist.sort(key=natural_keys) sorts in human order
        http://nedbatchelder.com/blog/200712/human_sorting.html
        (See Toothy's implementation in the comments)
        '''
        return [ self.atoi(c) for c in re.split(r'(\d+)', text) ]
    
    def rolling_average(self,x,N):
        """Rolling average function

        Args:
            x (array-like): array to be smoothed 
            N (int): window size

        Returns:
            list: smoothed array
        """
        indx = 0
        res = []
        while indx+N < len(x):
            xt = x[indx:indx+N]
            res.append(np.nanmean(xt))
            indx+=1
        
        for i in range(len(x[indx:])):
            res.append(np.nanmean(x[indx+i:]))

        res = np.array(res)
        x = np.array(x)
        res[x == 0] = 0

        return res
    
    def get_data(self,plate_path):
        # print(plate_path)
        if plate_path[-4:] == '.csv':
            df = pd.read_csv(plate_path)
        elif plate_path[-5:] == '.xlsx':
            df = pd.read_excel(plate_path)
        index_column = 'FileName_gfp'
        data = {}

        green_key = self.green_key
        red_key = self.red_key

        for i in df.index:
            cur_row  = df.iloc[i]
            name = cur_row[index_column]
            well_indx = name.find('_')
            well = cur_row[index_column][0:well_indx]
            if well in data.keys():
                # ts_green = cur_row[green_key]
                # ts_red = cur_row[red_key]
                # if self.smooth_counts:
                #     ts_green = self.rolling_average(ts_green,self.rolling_average_window)
                #     ts_red = self.rolling_average(ts_red,self.rolling_average_window)
                data[well]['green'].append(cur_row[green_key])
                data[well]['red'].append(cur_row[red_key])
            else:
                dict_t = {'green':[cur_row[green_key]],
                        'red':[cur_row[red_key]]}
                data[well] = dict_t
        
        if self.smooth_counts:
            for key in data.keys():
                data[key]['green'] = self.rolling_average(data[key]['green'],self.rolling_average_window)
                data[key]['red'] = self.rolling_average(data[key]['red'],self.rolling_average_window)

        return data
    
    def plot_plate(self,data=None,plate_num=0,plot_fit=False,ylim=[3,10],
                   show_title=False):
        
        if data is None:
            data = self.data_list[plate_num]

        fig,ax_list = plt.subplots(nrows=6,ncols=9,figsize=(12,8),
                                   sharex=True,sharey=True)

        fig.tight_layout()
        
        key0 = list(data.keys())[0]
        time = np.arange(len(data[key0]['green']))*self.dt

        row_list = ['B','C','D','E','F','G']
        col_list = np.arange(9) + 2
        col_list = [str(c) for c in col_list]

        row_indx = 0
        for row in row_list:
            col_indx = 0
            
            for col in col_list:
                ax = ax_list[row_indx,col_indx]
                key = row+col
                green_data = data[key]['green']
                red_data = data[key]['red']
                ax.plot(time,np.log(green_data),color='tab:cyan',linewidth=2.5)
                ax.plot(time,np.log(red_data),color='tab:brown',linewidth=2.5)
                if plot_fit:

                    if np.max(green_data) > 10**3:
                        fit = self.est_linear_slope(green_data,return_fit=True,**self.slope_est_options)
                        
                        if type(fit) == tuple:
                            fit_t = time*fit[0] + fit[1]
                            ax.plot(time,fit_t,color='black',alpha=0.9)

                    if np.max(red_data) > 10**3:
                        fit = self.est_linear_slope(red_data,return_fit=True,**self.slope_est_options)
                        
                        if type(fit) == tuple:
                            fit_t = time*fit[0] + fit[1]
                            ax.plot(time,fit_t,'--',color='black',alpha=0.9)
                    # if not np.isnan(fit):
                    #     ax.plot(time,fit,color='black')
                    # fit = self.est_linear_slope(red_data,return_fit=True)
                    # if type(fit) == tuple:
                    #     ax.plot(fit[0],fit[1],'--',color='black')
                    # if not np.isnan(fit):
                    #     ax.plot(time,fit,'--',color='black')
                # ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
                # ax.set_yscale('log')
                # ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
                if show_title:
                    ax.set_title(key)
                ax.set_ylim(ylim)

                if (plate_num,key) in self.exclude_wells:
                    ax.set_facecolor('grey')

                col_indx += 1
            row_indx += 1

        return fig,ax
    
    def est_AUC(self,counts,debug=False,time=None):
        """Estimates area under the curve of a growth curve

        Args:
            counts (array-like): Raw optical density data
            time (array-like): Time data

        Returns:
            float: Area under the curve
        """
        if time is None:
            time = np.arange(len(counts))*self.dt
        AUC = np.trapz(counts,time)
        return AUC
    
    def est_linear_slope(self,counts,
                        # dt=1, # per hour
                        dt=None,
                        time=None,
                        window=10,
                        step=5,
                        thresh=0.65,
                        debug=True,
                        title=None,
                        exclude=0,
                        cell_death_thresh=0.7,
                        return_fit=False):
        """A non-parametric method for estimating the slope of growth curves

        Args:
            OD (array-like): Raw optical density data
            window (float,optional): Window for rolling linear regression. Defaults to 10.
            step (float, optional): Step size for rolling lienar regression. Defaults to 5.
            thresh (float, optional): Percentage of max slope to include in final 
            linear regression. Defaults to 0.65.
            debug (bool, optional): If True, plots linear fit over data. Defaults to False.

        Returns:
            _type_: _description_
        """

        if dt is None:
            dt = self.dt

        n_measurements = len(counts)
        if type(counts) is list:
            counts = np.array(counts)

        if time is None:
            time = np.arange(len(counts))

        # remove zero values
        time = np.delete(time,np.argwhere(counts<=0))
        counts = np.delete(counts,np.argwhere(counts<=0))
        if len(counts) < n_measurements/2:
            return np.nan
        
        lnCount = np.log(counts) # normalize OD data

        # check if the population is decreasing
        if np.mean(lnCount[-3:]) < np.max(lnCount)*cell_death_thresh:
            return 0

        slopes = []
        x_pos = []

        # calculate piecewise slopes
        for indx in range(exclude,len(lnCount)-window+step,step):
            time_t = time[indx:indx+window]
            subset = lnCount[indx:indx+window]
            fit = scipy.stats.linregress(time_t,subset)
            slopes.append(fit.slope)
            x_pos.append(indx)

        # print(slopes)
        if all(np.array(slopes)<=0):
            return 0
        # find the max slope
        
        lb = thresh*np.nanmax(slopes)
        # print(lb)
        # print(slopes)
        use_indices = np.argwhere(slopes>=lb)[:,0]

        if len(use_indices) == 1:

            return np.nanmax(slopes)
        # print(use_indices)
        # print(x_pos)
        if len(use_indices) > 1:
            lin_range = x_pos[np.min(use_indices):np.max(use_indices)+1]

        else:
            lin_range = x_pos[use_indices[0]:use_indices[0]+1]
        # print(lin_range)

        # compute slope and plot result
        time_t = time[lin_range]
        lin_seg = lnCount[lin_range]
        # print(lin_seg)
        fit = scipy.stats.linregress(time_t,lin_seg)
        slope = fit.slope

        count_fit = time_t*fit.slope + fit.intercept
        
        if return_fit:
            return (fit.slope/dt,fit.intercept)

        if np.isnan(slope):
            # raise Warning('Slope is NaN, adjust parameters')
            return slope

        if self.growth_rate_type == 'linear':
            rate = np.exp(slope/dt)
        else:
            rate = slope/dt

        # plot the linear regression
        if debug:
            fig,ax_list = plt.subplots(ncols=2,figsize=(8,3))
            ax = ax_list[0]
            ax.plot(time,lnCount,linewidth=2,label='lnCount')
            ax.plot(time_t,count_fit,linewidth=2,label='Linear fit')
            # ax.plot(time[x_pos],(10000*np.array(slopes))+np.min(lnCount),label='slope')
            # ax.plot(cutoff_time,cutoff,label='Threshold')
            ax.legend(frameon=False)

            # title = title + ' ' + 'err = ' + str(round(100000*fit.stderr,2))

            ax.set_title(title)

            ax = ax_list[1]
            ax.plot(time,counts)
            ax.set_title('Raw count')
            # ax.set_ylim(0,1.1)
            fig.suptitle('Slope = ' + str(round(slope,3)), fontsize=14)
            AUC = np.trapz(counts,time)
            max_count = np.max(counts)
            AUC_ratio = AUC/np.trapz(max_count*np.ones(len(time)),time)
            ax.annotate('AUC ratio = ' + str(round(AUC_ratio,3)),xy=(0.5,0.5),xycoords='axes fraction')

        return rate # per hour
    
    def plate_fitness(self,data=None,plate_path=None,row_list=None,col_list=None):
        """Calculate fitness for each condition in a plate
        """
        if data is None:
            if plate_path is not None:
                data = self.get_data(plate_path)
            else:
                raise ValueError('Must provide data or plate_path')
        
        if row_list is None:
            row_list = ['B','C','D','E','F','G'] # assume exclude outer wells
        if col_list is None:
            col_list = ['2','3','4','5','6','7','8','9','10','11']
        
        res = {}

        for row in row_list:
            for col in col_list:
                well = row+col
                res[well] = self.get_fitness(data,well)
        
        

    def get_fitness_lib(self, plate_paths=None, proportions=None, 
        num_replicates=None, slope_kwargs=None, debug=False):
        """Generate library of fitness values for each condition

        Args:
            plate_paths (list, optional): List of plate paths. Defaults to None.
            proportions (list, optional): List of proportions. Defaults to None.
            num_replicates (int, optional): Number of replicates. Defaults to 3.
            slope_kwargs (dict, optional): Slope estimation parameters. Defaults to None.
            debug (bool, optional): If True, plots linear fit over data. Defaults to False.

        Returns:
            tuple: (fitness_res_avg, fitness_parent_avg)
        """
        if proportions is None:
            proportions = self.proportions
        if plate_paths is None:
            plate_paths = self.plate_paths
        if slope_kwargs is None:
            slope_kwargs = self.slope_est_options
        if num_replicates is None:
            num_replicates = self.num_replicates

        col_list = np.arange(9) + 2
        col_list = [str(c) for c in col_list]
        num_plates = len(plate_paths)

        if len(proportions) == num_plates * 2:
            proportions = proportions * num_replicates
            proportions.sort()
            proportions.reverse()

        plate_indx = 0
        fitness_values = {}
        prop_indx = 0  # index for proportion vector
        
        if self.fitness_estimate == 'slope':
            fitness_func = self.est_linear_slope
        elif self.fitness_estimate == 'AUC' or self.fitness_estimate == 'auc':
            fitness_func = self.est_AUC
            slope_kwargs = {}
        row_list = ['B', 'C', 'D', 'E', 'F', 'G']

        for plate in plate_paths:
            data = self.get_data(plate)
            
            for row in row_list:
                prop = proportions[prop_indx]
                fitness_resistant = []
                fitness_parental = []

                for col in col_list:
                    key = row + col
                    if not (plate_indx, key) in self.exclude_wells:
                        count_res = data[key]['red']
                        count_parent = data[key]['green']

                        if prop == 0:
                            if np.sum(count_parent) > 0:
                                # s = self.est_linear_slope(count_parent, debug=debug, **slope_kwargs)
                                s = fitness_func(count_parent, debug=debug, **slope_kwargs)

                                if self.fitness_estimate == 'AUC' or self.fitness_estimate == 'auc':
                                    s = s / np.mean(count_parent[0:5])   
                                    s = s / (self.dt*len(count_parent)) # convert to per hour
                                fitness_parental.append(s)
                            else:
                                fitness_parental.append(np.nan)
                        elif prop == 1:
                            if np.sum(count_res) > 0:
                                # s = self.est_linear_slope(count_res, debug=debug, **slope_kwargs)
                                s = fitness_func(count_res, debug=debug, **slope_kwargs)

                                if self.fitness_estimate == 'AUC' or self.fitness_estimate == 'auc':
                                    s = s / np.mean(count_res[0:5])
                                    s = s / (self.dt*len(count_res)) # convert to per hour

                                fitness_resistant.append(s)
                            else:
                                fitness_resistant.append(np.nan)
                        else:
                            if np.sum(count_parent) > 0:
                                # s_par = self.est_linear_slope(count_parent, debug=debug, **slope_kwargs)
                                s_par = fitness_func(count_parent, debug=debug, **slope_kwargs)
                                
                                if self.fitness_estimate == 'AUC' or self.fitness_estimate == 'auc':
                                    s_par = s_par / np.mean(count_parent[0:5])
                                    s_par = s_par / (self.dt*len(count_parent)) # convert to per hour

                                fitness_parental.append(s_par)
                            else:
                                fitness_parental.append(np.nan)

                            if np.sum(count_res) > 0:
                                # s_res = self.est_linear_slope(count_res, debug=debug, **slope_kwargs)
                                s_res = fitness_func(count_res, debug=debug, **slope_kwargs)
                                
                                if self.fitness_estimate == 'AUC' or self.fitness_estimate == 'auc':
                                    s_res = s_res / np.mean(count_res[0:5])
                                    s_res = s_res / (self.dt*len(count_res)) # convert to per hour
                                
                                fitness_resistant.append(s_res)
                            else:
                                fitness_resistant.append(np.nan)
                    else:
                        fitness_resistant.append(np.nan)
                        fitness_parental.append(np.nan)
                dict_t = {'resistant': fitness_resistant, 'parental': fitness_parental}
                fitness_values[prop_indx] = dict_t
                prop_indx += 1

            plate_indx += 1

        fitness_res_avg = {}
        fitness_parent_avg = {}

        for i in range(0, len(proportions), num_replicates):
            prop = proportions[i]

            if prop == 1:
                fitness_t = np.zeros((len(col_list), num_replicates))

                for k in range(0, num_replicates):
                    key = i + k
                    fitness_t[:, k] = np.array(fitness_values[key]['resistant'])

                fitness_avg = np.nanmean(fitness_t, axis=1)
                fitness_std = np.nanstd(fitness_t, axis=1)
                fitness_err = fitness_std / np.sqrt(num_replicates)
                dict_t = {'avg': fitness_avg, 'err': fitness_err}
                fitness_res_avg[prop] = dict_t

            elif prop == 0:
                fitness_t = np.zeros((len(col_list), num_replicates))

                for k in range(0, num_replicates):
                    key = i + k
                    fitness_t[:, k] = np.array(fitness_values[key]['parental'])

                fitness_avg = np.nanmean(fitness_t, axis=1)
                fitness_std = np.nanstd(fitness_t, axis=1)
                fitness_err = fitness_std / np.sqrt(num_replicates)
                dict_t = {'avg': fitness_avg, 'err': fitness_err}
                fitness_parent_avg[prop] = dict_t

            else:
                fitness_t_resistant = np.zeros((len(col_list), num_replicates))
                fitness_t_parental = np.zeros((len(col_list), num_replicates))

                for k in range(0, num_replicates):
                    key = i + k
                    fitness_t_resistant[:, k] = np.array(fitness_values[key]['resistant'])
                    fitness_t_parental[:, k] = np.array(fitness_values[key]['parental'])

                fitness_avg_resistant = np.nanmean(fitness_t_resistant, axis=1)
                fitness_std_resistant = np.nanstd(fitness_t_resistant, axis=1)
                fitness_err_resistant = fitness_std_resistant / np.sqrt(num_replicates)
                dict_t_resistant = {'avg': fitness_avg_resistant, 'err': fitness_err_resistant}
                fitness_res_avg[prop] = dict_t_resistant

                fitness_avg_parental = np.nanmean(fitness_t_parental, axis=1)
                fitness_std_parental = np.nanstd(fitness_t_parental, axis=1)
                fitness_err_parental = fitness_std_parental / np.sqrt(num_replicates)
                dict_t_parental = {'avg': fitness_avg_parental, 'err': fitness_err_parental}
                fitness_parent_avg[prop] = dict_t_parental

        return fitness_res_avg, fitness_parent_avg


    def plot_dose_response_curves(self,gr_parent=None,gr_res=None):
        
        cmap = mpl.colormaps['viridis']
        
        if gr_parent is None:
            if self.fitness_parental is None:
                raise TypeError('No parental growth rate data to plot.')
            else:
                gr_parent = self.fitness_parental
        if gr_res is None:
            if self.fitness_resistant is None:
                raise TypeError('No resistant growth rate data to plot.')
            else:
                gr_res = self.fitness_resistant
        
        if self.fitness_estimate == 'slope':
            ylabel = 'Growth rate (hr$^{-1}$)'
        elif self.fitness_estimate == 'AUC' or self.fitness_estimate == 'auc':
            ylabel = 'Relative AUC'
        fig,ax_list = plt.subplots(ncols=2,figsize=(8,3))

        ax = ax_list[0]
        indx = 0

        for key in gr_parent.keys():
            yt = np.array(gr_parent[key]['avg'])
            yerr = np.array(gr_parent[key]['err'])
            x = np.arange(len(yt))
            ax.errorbar(x,yt,color=cmap(indx/6),yerr=yerr,label=self.proportions[indx+1])
            indx+=1



        ax.legend(frameon=False,title='proportion \nmutant',ncol=1,loc=(1,0.2))
        # ax.set_title('BRAF vs parental')
        ax.set_xticks(np.arange(9))
        dc_log = [np.log2(c) for c in self.dc]
        dc_log[0] = 'nd'
        ax.set_xticklabels(dc_log)

        ax.set_xlabel('Log$_{2}$ drug concentration (relative to max)',fontsize=12)
        # ax.set_ylabel('Parental growth rate',fontsize=12)
        ax.set_ylabel(ylabel,fontsize=12)
        ax.set_title('Parental',fontsize=12)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # ax.set_ylim(0.05,0.14)

        # fig,ax = plt.subplots()
        ax = ax_list[1]
        indx = 0
        for key in gr_res.keys():
            yt = np.array(gr_res[key]['avg'])
            yerr = np.array(gr_res[key]['err'])
            x = np.arange(len(yt))
            ax.errorbar(x,yt,color=cmap(indx/6),yerr=yerr,label=self.proportions[indx])
            indx+=1

        yl_min = np.inf
        yl_max = 0
        for ax in ax_list:
            yl = ax.get_ylim()
            if np.min(yl) < yl_min:
                yl_min = np.min(yl)
            if np.max(yl) > yl_max:
                yl_max = np.max(yl)

        for ax in ax_list:
            ax.set_ylim(yl_min,yl_max)

        # put the legend between the two axes
        ax.legend(frameon=False,title='proportion \nmutant',ncol=1,loc=(1,0.2))
        # ax.set_title('BRAF vs parental')
        ax.set_xticks(np.arange(9))
        dc_log = [np.log2(c) for c in self.dc]
        dc_log[0] = 'nd'
        ax.set_xticklabels(dc_log)

        ax.set_xlabel('Log$_{2}$ drug concentration (relative to max)',fontsize=12)
        ax.set_ylabel(ylabel,fontsize=12)
        # ax.set_ylabel('Mutant growth rate',fontsize=12)
        ax.set_title('Mutant',fontsize=12)
        # ax.set_ylim(0.05,0.14)
        # ax.set_yscale('log')

        pos = ax.get_position()
        pos.x0 = pos.x0 + 0.15
        pos.x1 = pos.x1 + 0.15
        ax.set_position(pos)

        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        return fig,ax
    
    def hill_fn(self,conc,gmax, gmin, hc, ic_50):
        
        y = gmax + ((gmin - gmax) * conc**hc) / (ic_50**hc + conc**hc)
        return y

    def fds_model(self,frac, g_100, alpha, beta):
        g_100 = g_100 # Growth rate at 100% fraction
        g_0 = (g_100*(1+alpha) + beta) # Maximally moderated growth ( 0% frac)
        g_frac = frac*g_100 + (1-frac)*g_0 # Growth rate (fraction)
        return g_frac
    
    def fds_model_curve_fit(self,X,alpha,beta):
        """FDS model for curve fitting

        Args:
            X (list): [g_100,frac]
            alpha (float): alpha
            beta (float): beta

        Returns:
            float: growth rate
        """
        g_100 = X[0]
        frac = X[1]
        g_0 = (g_100*(1+alpha) + beta) # Maximally moderated growth ( 0% frac)
        g_frac = frac*g_100 + (1-frac)*g_0 # Growth rate (fraction)
        return g_frac
    
    def plot_slope_vs_g100(self,which='parental'):
        if which == 'mutant':
            data = self.fitness_resistant
            g_100_indx = 1
            proportions = self.proportions[0:-1]
            # prop_mutant = self.proportions[0:-1]
        elif which == 'parental':
            data = self.fitness_parental
            g_100_indx = 0
            proportions = np.array(self.proportions[1:])
            # prop_mutant = self.proportions[1:]
            proportions = np.flip(proportions)
        
        fig,ax = plt.subplots()

        fig2,ax_list = plt.subplots(ncols=3,nrows=3,figsize=(12,12),sharex=True,sharey=True)
        ax_list = ax_list.flatten()

        g_100 = []
        slope = []
        slope_err = []
        for c_indx in range(len(self.dc)):
            g_100.append(data[g_100_indx]['avg'][c_indx])
            g = []
            g_err = []
            for p in proportions:
                g.append(data[p]['avg'][c_indx])
                g_err.append(data[p]['err'][c_indx])
            # linear regression
            fit = scipy.stats.linregress(proportions,g)
            slope.append(fit.slope)
            slope_err.append(fit.stderr)

            ax_t = ax_list[c_indx]
            ax_t.errorbar(proportions,g,yerr=g_err,fmt='o',color='black')
            x = np.linspace(0,1,100)
            y = fit.slope*x + fit.intercept
            ax_t.plot(x,y,color='r')

        ax.errorbar(g_100,slope,yerr=slope_err,fmt='o')
        ax.set_xlabel('Monotypic growth rate',fontsize=12)
        ax.set_ylabel('Frequency-dependent interaction slope',fontsize=12)
        ax.set_title(which,fontsize=12)

        return fig,ax
    
    def plot_slope_vs_conc(self,which='parental',ylim=[-0.007,0.004],ax=None,fig=None):
        if which == 'mutant':
            data = self.fitness_resistant
            g_100_indx = 1
            proportions = self.proportions[0:-1]
            # prop_mutant = self.proportions[0:-1]
        elif which == 'parental':
            data = self.fitness_parental
            g_100_indx = 0
            proportions = np.array(self.proportions[1:])
            # prop_mutant = self.proportions[1:]
            proportions = np.flip(proportions)

        if ax is None:
            fig,ax = plt.subplots()

        # fig2,ax_list = plt.subplots(ncols=3,nrows=3,figsize=(12,12),sharex=True,sharey=True)
        # ax_list = ax_list.flatten()

        g_100 = []
        slope = []
        slope_err = []
        for c_indx in range(len(self.dc)):
            g_100.append(data[g_100_indx]['avg'][c_indx])
            g = []
            g_err = []
            for p in proportions:
                g.append(data[p]['avg'][c_indx])
                g_err.append(data[p]['err'][c_indx])
            # linear regression
            fit = scipy.stats.linregress(proportions,g)
            slope.append(fit.slope)
            slope_err.append(fit.stderr)

            # ax_t = ax_list[c_indx]
            # ax_t.errorbar(proportions,g,yerr=g_err,fmt='o',color='black')
            # x = np.linspace(0,1,100)
            # y = fit.slope*x + fit.intercept
            # ax_t.plot(x,y,color='r')

        ax.errorbar(self.dc,slope,yerr=slope_err,fmt='o')
        ax.set_xscale('log')
        ax.set_xlabel('Drug concentration',fontsize=14)
        ax.set_ylabel('Frequency-dependent interaction slope',fontsize=14)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_title(which,fontsize=14)
        ax.set_ylim(ylim)
        ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))

        return fig,ax


    def game_objective_function(self,data,proportions,g_100_indx,alpha_est,beta_est):
        
        loss = 0        

        for c_indx in range(len(self.dc)):
            
            g_100 = data[g_100_indx]['avg'][c_indx]
            for prop in proportions:
                gr_t = data[prop]['avg'][c_indx]
                gr_est = self.fds_model(prop,g_100,alpha_est,beta_est)
                err = (gr_t - gr_est)**2
                loss += err

        return loss

    def estimate_game_params(self,which='mutant',debug=False,save_fig=False,
                             save_path=None,bounds=None,p0=None):
        
        if which == 'mutant':
            data = self.fitness_resistant
            g_100_indx = 1
            proportions = self.proportions[0:-1]
            prop_mutant = self.proportions[0:-1]
        elif which == 'parental':
            data = self.fitness_parental
            g_100_indx = 0
            proportions = np.array(self.proportions[1:])
            prop_mutant = self.proportions[1:]
            proportions = np.flip(proportions)

        # Fit FDS model
        if p0 is None:
            p0 = [0,0]
        # bounds = ([-1,-.1],[1,.1])
        if bounds is None:
            bounds = ([-1,-1],[1,1])
        # res = sciopt.minimize(lambda param_est: 
        #                       self.game_objective_function(data,prop_mutant,g_100_indx,param_est[0],param_est[1]), 
        #                       p0, bounds=bounds)
        
        # try with scipy curve fit
        g_100 = data[g_100_indx]['avg']
        # generate X and y data
        X = [[],[]]
        y = []
        for conc_indx in range(len(self.dc)):
            for frac in proportions:
                X[0].append(g_100[conc_indx])
                X[1].append(frac)
                y.append(data[frac]['avg'][conc_indx])

        popt,pcov = sciopt.curve_fit(self.fds_model_curve_fit, X, y, p0=p0, bounds=bounds)

        alpha_est = popt[0]
        beta_est = popt[1]
        alpha_err = np.sqrt(pcov[0,0])
        beta_err = np.sqrt(pcov[1,1])
        
        # alpha_est = res.x[0]
        # beta_est = res.x[1]
        # alpha_err = np.sqrt(res.hess_inv.todense()[0,0])
        # beta_err = np.sqrt(res.hess_inv.todense()[1,1])

        res = {'alpha':alpha_est,'beta':beta_est,'alpha_err':alpha_err,'beta_err':beta_err}

        if debug:
            fig,ax = plt.subplots()
            for c_indx in range(len(self.dc)):
                
                gr_t = []
                gr_err = []
                for prop in prop_mutant:
                    gr_t.append(data[prop]['avg'][c_indx])
                    gr_err.append(data[prop]['err'][c_indx])
                gr_t = np.array(gr_t)
                gr_err = np.array(gr_err)
                if which == 'mutant':
                    g_100 = gr_t[0]
                elif which == 'parental':
                    g_100 = gr_t[-1]
                # g_100 = gr_t[0]
                ax.errorbar(proportions,gr_t,yerr=gr_err,fmt='o',color=self.cmap(c_indx/len(self.dc)))
                yfit = []
                for frac in proportions:
                    yfit.append(self.fds_model(frac,g_100,alpha_est,beta_est))
                ax.plot(proportions,yfit,color=self.cmap(c_indx/len(self.dc)),label=np.log2(self.dc[c_indx]))
                # ax.set_title(str(round(popt[0],3)))
            if which == 'mutant':
                ax.set_xlabel('Fraction of mutant cells',fontsize=14)
            elif which == 'parental':
                ax.set_xlabel('Fraction of parental cells',fontsize=14)
            # ax.set_xlabel('Fraction of resistant cells',fontsize=14)
            loss = self.game_objective_function(data,prop_mutant,g_100_indx,alpha_est,beta_est)
            ax.set_ylabel('Growth rate',fontsize=14)
            ax.set_title('alpha = ' + str(round(alpha_est,3)) + ', beta = ' + str(round(beta_est,3)) + ', loss = ' + str(round(loss,3)),fontsize=14)
            ax.legend(frameon=False,fontsize=12,loc=(1.05,.5),title='drug concentration')
            ax.tick_params(labelsize=12)
            fig.suptitle(which,fontsize=14,y=1.02)
            if save_fig:
                if save_path is None:
                    save_path = 'figures/'
                fig.savefig(save_path + which + '_game_fit.png',dpi=300,bbox_inches='tight')
        return res
    
    def plot_linear_game(self,game_params=None,hill_params=None,which='parental',
                      style='heatmap',cmap='magma',save_fig=False,save_path=None):
        if game_params is None:
            if which == 'parental':
                game_params = self.parental_game_params
            elif which == 'mutant':
                game_params = self.mutant_game_params
        if hill_params is None:
            if which == 'parental':
                hill_params = self.parental_dr_params[0]
            elif which == 'mutant':
                hill_params = self.mutant_dr_params[1]
            
        if which == 'parental':
            proportions = np.array(self.proportions[1:])
            proportions = np.flip(proportions)
            data = self.fitness_parental
            xlabel = 'Fraction of parental cells'
        elif which == 'mutant':
            proportions = np.array(self.proportions[0:-1])
            data = self.fitness_resistant
            xlabel = 'Fraction of mutant cells'
        
        fig,ax_list = plt.subplots(ncols=3,nrows=3,figsize=(12,12),sharex=True,sharey=True)
        ax_list = ax_list.flatten()
        for c_indx in range(len(self.dc)):
            g_100 = data[proportions[0]]['avg'][c_indx]
            ax = ax_list[c_indx]
            gr_t = []
            gr_err = []
            # gr_fit = []
            for prop in proportions:
                gr_t.append(data[prop]['avg'][c_indx])
                gr_err.append(data[prop]['err'][c_indx])
                if which == 'parental':
                    prop_fit = 1 - prop
                else:
                    prop_fit = prop
                # gr_fit.append(self.fds_model(prop_fit,g_100,
                #                                 game_params['alpha'],
                #                                 game_params['beta']))

            res = scipy.stats.linregress(proportions,gr_t)
            gr_fit = res.intercept + res.slope*proportions
            ax.errorbar(proportions,gr_t,yerr=gr_err,fmt='o',color='k',capsize=5)
            # plot game fit
            ax.plot(proportions,gr_fit,'-',color='r')
            ax.set_title('Drug conc. = ' + str(self.dc[c_indx]),fontsize=12)
            # annotate with R^2
            # convert lists to arrays
            gr_t = np.array(gr_t)
            gr_fit = np.array(gr_fit)

            # try r_sqr from correlation coefficient
            # normalize data

            r,p = scipy.stats.pearsonr(gr_t,gr_fit)
            r_sqr = r**2

            # ax.annotate('$R^2$ = ' + str(np.round(r_sqr,2)),xy=(0.05,0.85),xycoords='axes fraction',fontsize=12)
            # use legend to place the annotate in the best position
            handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                lw=0, alpha=0)]*2
            ax.legend(handles,['$R^2$ = ' + str(np.round(r_sqr,2)),'p = ' + str(np.round(p,5))],
                        loc='best',fontsize=12,frameon=False)
        # add labels to the left column
        for ax in ax_list[0::3]:
            if self.fitness_estimate == 'AUC':
                yl = 'AUC'
            else:
                yl = 'Growth rate'
            ax.set_ylabel(yl,fontsize=12)
        # add labels to the bottom row
        for ax in ax_list[-3:]:
            ax.set_xlabel(xlabel,fontsize=12)
        fig.tight_layout()
        fig.suptitle(which + ' game fit',fontsize=14,y=1.02)
        return fig,ax_list
    
    def plot_game_fit(self,game_params=None,hill_params=None,which='parental',
                      style='heatmap',cmap='magma',save_fig=False,save_path=None):
        if game_params is None:
            if which == 'parental':
                game_params = self.parental_game_params
            elif which == 'mutant':
                game_params = self.mutant_game_params
        if hill_params is None:
            if which == 'parental':
                hill_params = self.parental_dr_params[0]
            elif which == 'mutant':
                hill_params = self.mutant_dr_params[1]
            
        if which == 'parental':
            proportions = np.array(self.proportions[1:])
            proportions = np.flip(proportions)
            data = self.fitness_parental
            xlabel = 'Fraction of parental cells'
        elif which == 'mutant':
            proportions = self.proportions[0:-1]
            data = self.fitness_resistant
            xlabel = 'Fraction of mutant cells'

        if style == 'heatmap':
            fig,ax_list = plt.subplots(ncols=3,figsize=(12,4))
            data_matrix = np.zeros((len(self.dc),len(proportions)))
            fit_matrix = np.zeros((len(self.dc),len(proportions)))
            for c_indx in range(len(self.dc)):
                g_100 = data[proportions[0]]['avg'][c_indx]
                for p_indx in range(len(proportions)):
                    prop = proportions[p_indx]
                    data_matrix[c_indx,p_indx] = data[prop]['avg'][c_indx]
                    fit_matrix[c_indx,p_indx] = self.fds_model(prop,g_100,game_params['alpha'],game_params['beta'])

            residuals = (data_matrix - fit_matrix)
            vmin = np.min((data_matrix,fit_matrix))
            vmax = np.max((data_matrix,fit_matrix))
            sns.heatmap(residuals,ax=ax_list[2],cmap=cmap)
            sns.heatmap(data_matrix,ax=ax_list[0],cmap=cmap,vmin=vmin,vmax=vmax)
            sns.heatmap(fit_matrix,ax=ax_list[1],cmap=cmap,vmin=vmin,vmax=vmax)

            ax_list[0].set_title('Data',fontsize=14)
            ax_list[1].set_title('Fit',fontsize=14)
            ax_list[2].set_title('Residuals',fontsize=14)

            for ax in ax_list:
                ax.set_xticks(np.arange(len(proportions))+.5)
                ax.set_xticklabels(np.round(proportions,2),rotation=45)
                ax.set_yticks(np.arange(len(self.dc))+.5)
                ax.set_yticklabels(np.log2(self.dc),rotation=45)
                ax.set_xlabel(xlabel,fontsize=12)
                ax.set_ylabel('Drug concentration',fontsize=12)
                ax.tick_params(labelsize=12)
            
            fig.tight_layout()
            fig.suptitle(which + ' game fit',fontsize=14,y=1.05)
        
        elif style == 'cartesian':
            fig,ax_list = plt.subplots(ncols=3,nrows=3,figsize=(12,12),sharex=True,sharey=True)
            ax_list = ax_list.flatten()
            for c_indx in range(len(self.dc)):
                g_100 = data[proportions[0]]['avg'][c_indx]
                ax = ax_list[c_indx]
                gr_t = []
                gr_err = []
                gr_fit = []
                for prop in proportions:
                    gr_t.append(data[prop]['avg'][c_indx])
                    gr_err.append(data[prop]['err'][c_indx])
                    if which == 'parental':
                        prop_fit = 1 - prop
                    else:
                        prop_fit = prop
                    gr_fit.append(self.fds_model(prop_fit,g_100,
                                                 game_params['alpha'],
                                                 game_params['beta']))
                ax.errorbar(proportions,gr_t,yerr=gr_err,fmt='o',color='k',capsize=5)
                # plot game fit
                ax.plot(proportions,gr_fit,'-',color='r')
                ax.set_title('Drug conc. = ' + str(self.dc[c_indx]),fontsize=12)
                # annotate with R^2
                # convert lists to arrays
                gr_t = np.array(gr_t)
                gr_fit = np.array(gr_fit)

                # try r_sqr from correlation coefficient
                # normalize data

                r,p = scipy.stats.pearsonr(gr_t,gr_fit)
                r_sqr = r**2

                # ax.annotate('$R^2$ = ' + str(np.round(r_sqr,2)),xy=(0.05,0.85),xycoords='axes fraction',fontsize=12)
                # use legend to place the annotate in the best position
                handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", 
                                 lw=0, alpha=0)]*2
                ax.legend(handles,['$R^2$ = ' + str(np.round(r_sqr,2)),'p = ' + str(np.round(p,5))],
                          loc='best',fontsize=12,frameon=False)
            # add labels to the left column
            for ax in ax_list[0::3]:
                if self.fitness_estimate == 'AUC':
                    yl = 'AUC'
                else:
                    yl = 'Growth rate'
                ax.set_ylabel(yl,fontsize=12)
            # add labels to the bottom row
            for ax in ax_list[-3:]:
                ax.set_xlabel(xlabel,fontsize=12)
            fig.tight_layout()
            fig.suptitle(which + ' game fit',fontsize=14,y=1.02)
        return fig,ax_list

    def fit_dr_fn(self,which='mutant',debug=False,mutant_frac=None):

        if which == 'mutant':
            if mutant_frac is None:
                mutant_frac = 1
            data = self.fitness_resistant
            gr_0 = data[mutant_frac]['avg']
            err = data[mutant_frac]['err']

            if self.mutant_linear_fit:
                linear_fit = True
            else:
                linear_fit = False

        elif which == 'parental':
            if mutant_frac is None:
                mutant_frac = 0
            data = self.fitness_parental
            gr_0 = data[mutant_frac]['avg']
            err = data[mutant_frac]['err']

            if self.parental_linear_fit:
                linear_fit = True
            else:
                linear_fit = False
        
        if not linear_fit:
            gmax_est= np.max(gr_0)
            gmin_est = np.min(gr_0)
            hc_est = 1
            ic_50_est = 1
            p0 = [gmax_est,gmin_est,hc_est,ic_50_est]
            popt,pcov = sciopt.curve_fit(self.hill_fn,self.dc,gr_0,p0=p0,sigma=err,
                                        maxfev=10000)
            # check if the values are reasonable
            
            # if popt[0] < popt[1]:
            #     raise ValueError('gmax < gmin')
            # if popt[2] < 0:
            #     raise ValueError('hc < 0')
            # if popt[3] < 0 or popt[3] > np.max(self.dc):
            #     raise ValueError('ic_50 < 0 or ic_50 > max(dc)')        
        
        if linear_fit:
            # log-linear fit to data
            res,pcov = sciopt.curve_fit(self.linear_fn,self.dc,gr_0,sigma=err,
                            maxfev=10000)
            popt = [res[0],res[1],0,0]

        if debug:
            xfit = np.logspace(np.min(self.dc),np.max(self.dc),100)
            # xfit = np.logspace(np.log10(np.min(self.dc)),np.log10(np.max(self.dc)),100)
            # xfit = np.lin
            if linear_fit:
                yfit = self.linear_fn(xfit,*popt)
            else:
                yfit = self.hill_fn(xfit,*popt)
                
            fig,ax = plt.subplots()
            ax.plot(xfit,yfit)
            ax.errorbar(self.dc,gr_0,yerr=err,fmt='o')
            ax.set_xscale('log')
            ax.set_title(which + ' fit,' + ' mutant prop = ' + str(mutant_frac))
            # annotate estimated params with legend trick
            handles = [mpl_patches.Rectangle((0, 0), 1, 1, fc="white", ec="white",
                                                lw=0, alpha=0)]*4
            
            ax.legend(handles,['gmax = ' + str(np.round(popt[0],2)),
                                 'gmin = ' + str(np.round(popt[1],2)),
                                 'hc = ' + str(np.round(popt[2],2)),
                                 'ic50 = ' + str(np.round(popt[3],2))],
                         loc='best',fontsize=10,frameon=False)
        return popt
    
    def linear_fn(self,x,m,b,c=None,d=None):
        x = np.log10(x)
        return m*x + b
    
    def fit_all_dr_curves(self,which='mutant',debug=False,separate_plots=False):
        if which == 'mutant':
            data = self.fitness_resistant
            if self.mutant_linear_fit:
                dr_func = self.linear_fn
            else:
                dr_func = self.hill_fn
        elif which == 'parental':
            data = self.fitness_parental
            if self.parental_linear_fit:
                dr_func = self.linear_fn
            else:
                dr_func = self.hill_fn

        popt_dict = {}
        for key in data.keys():
            popt_dict[key] = self.fit_dr_fn(which=which,mutant_frac=key)

        if debug:
            if separate_plots:
                for key in popt_dict.keys():
                    xfit = np.logspace(np.min(np.log10(self.dc)),
                                       np.max(np.log10(self.dc)),100)
                    yfit = dr_func(xfit,*popt_dict[key])
                    fig,ax = plt.subplots()
                    ax.plot(xfit,yfit)
                    ax.errorbar(self.dc,data[key]['avg'],yerr=data[key]['err'],
                                fmt='o')
                    ax.set_xscale('log')
            else:
                fig,ax = plt.subplots()
                for key in popt_dict.keys():
                    xfit = np.logspace(np.min(np.log10(self.dc)),
                                       np.max(np.log10(self.dc)),100)
                    yfit = dr_func(xfit,*popt_dict[key])

                    if which == 'parental':
                        cmap_val = key/0.9
                    else:
                        cmap_val = (key-0.05)/0.9

                    ax.plot(xfit,yfit,label=key,color=self.cmap(cmap_val))
                    ax.errorbar(self.dc,data[key]['avg'],yerr=data[key]['err'],
                                fmt='o',color=self.cmap(cmap_val))
                    ax.legend()
                    ax.set_xscale('log')
        return popt_dict
    
    def set_dose_response_params(self):
        self.mutant_dr_params = self.fit_all_dr_curves(which='mutant')
        self.parental_dr_params = self.fit_all_dr_curves(which='parental')

    def set_game_params(self):
        self.mutant_game_params = self.estimate_game_params(which='mutant')
        self.parental_game_params = self.estimate_game_params(which='parental')

    def unified_model(self,conc,f,game_params,hill_params,dr_curve='hill'):
        """Unified frequency-drug interaction model

        Args:
            conc (float): drug concentration
            f (float): fraction of cells of interest
            game_params (list): list of game parameters (alpha,beta)
            hill_params (list): 

        Returns:
            _type_: _description_
        """
        # get g_100 at conc
        if dr_curve == 'hill':
            g_100 = self.hill_fn(conc,*hill_params)
        elif dr_curve == 'linear':
            g_100 = self.linear_fn(conc,*hill_params)

        return self.fds_model(f,g_100,*game_params)

    def plot_unified_model(self,game_params=None,hill_params=None,conc_range=None,
                           ax=None,fig=None,cmap=None,legend=True,which='parental'):
        if game_params is None:
            if which == 'mutant':
                game_params = [self.mutant_game_params['alpha'],self.mutant_game_params['beta']]
            elif which == 'parental':
                game_params = [self.parental_game_params['alpha'],self.parental_game_params['beta']]
        if hill_params is None:
            if which == 'mutant':
                hill_params = self.mutant_dr_params[1]
            elif which == 'parental':
                hill_params = self.parental_dr_params[0]

        dr_curve = 'hill'
        if which == 'mutant':
            if self.mutant_linear_fit:
                dr_curve = 'linear'
        elif which == 'parental':
            if self.parental_linear_fit:
               dr_curve = 'linear'
        
        if conc_range is None:
            conc_range = np.logspace(-3,0,100)
        if ax is None:
            fig,ax = plt.subplots()
        if cmap is None:
            cmap = self.cmap
        for f in np.linspace(0,1,11):
            ax.plot(conc_range,self.unified_model(conc_range,f,game_params,
                                                  hill_params,dr_curve=dr_curve),
                    label=np.round(f,2),color=cmap(f),linewidth=2)
        if legend:
            ax.legend(frameon=False,fontsize=10)
        ax.set_xlabel('Drug concentration',fontsize=14)
        if self.fitness_estimate == 'AUC':
            yl = 'AUC'
        else:
            yl = 'Growth rate'
        ax.set_ylabel(yl,fontsize=14)
        ax.set_xscale('log')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        return fig,ax

    def plot_pharmacodynamic_curves(self,which='parental'):
        """Plot pharmacodynamic curves

        Args:
            which (str, optional): 'parental' or 'mutant'. Defaults to 'parental'.

        Returns:
            figure, axes tuple: figure and axes objects 
        """
        if which == 'parental':
            dr_params = self.parental_dr_params
            data = self.fitness_parental
            color_pad = 0.9
            if self.parental_linear_fit:
                dr_func = self.linear_fn
            else:
                dr_func = self.hill_fn
        elif which == 'mutant':
            dr_params = self.mutant_dr_params
            data = self.fitness_resistant
            color_pad = 1
            if self.mutant_linear_fit:
                dr_func = self.linear_fn
            else:
                dr_func = self.hill_fn


        fig,ax = plt.subplots()
        for key in dr_params.keys():
            # xfit = np.linspace(np.min(self.dc),np.max(self.dc),100)
            xfit = np.logspace(np.min(np.log10(self.dc)),
                               np.max(np.log10(self.dc)),100)
            yfit = dr_func(xfit,*dr_params[key])
            ax.plot(xfit,yfit,label=key,color=self.cmap(key/color_pad))
            ax.errorbar(self.dc,data[key]['avg'],yerr=data[key]['err'],
                        fmt='x',color=self.cmap(key/color_pad))
            ax.legend(frameon=False,title='proportion mutant',ncol=2)   
            ax.set_xscale('log')
            ax.spines['right'].set_visible(False)
            ax.spines['top'].set_visible(False)
            ax.tick_params(labelsize=14)
            ax.set_xlabel('Drug concentration',fontsize=14)
            if self.fitness_estimate == 'AUC':
                ax.set_ylabel('AUC',fontsize=14)
            else:
                ax.set_ylabel('Growth rate',fontsize=14)
            ax.set_title(which,fontsize=14)

        return fig,ax
    
    def spot_check_growth_curves(self,wells=None,plate_num=0,which='parental',prop=None,conc=None):
        """Spot check growth curves. Must specify either specific wells or a specific 
        propotion and drug concentration.

        Args:
            wells (list): list of well names
            which (str, optional): 'parental' or 'mutant'. Defaults to 'parental'.

        Returns:
            figure, axes tuple: figure and axes objects 
        """
        fig,ax = plt.subplots()
        if which == 'parental':
            color_key = 'green'
        elif which == 'mutant':
            color_key = 'red'

        if wells is None:
            # figure out which wells to plot from prop and conc
            if prop is None:
                raise TypeError('Must specify proportion of mutant cells if not specifying specific wells.')
            elif conc is None:
                raise TypeError('Must specify drug concentration if not specifying specific wells.')
                
        for well in wells:
            data = self.get_data(self.plate_paths[plate_num])
            ax.plot(data[well][color_key],label=well)
            ax.legend(frameon=False)
        return fig,ax
    
    def growth_rate_to_div_prob(self,growth_rate,dt=None):
        """Convert growth rate to division probability

        Args:
            growth_rate (float): growth rate

        Returns:
            float: division probability
        """
        if dt is None:
            dt = self.dt
        

        return 1 - np.exp(-np.exp(growth_rate)*dt)

