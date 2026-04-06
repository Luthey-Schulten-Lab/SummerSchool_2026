import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
from matplotlib.colors import to_rgba
import numpy as np
import h5py
import WCM_diagnosis as diagnosis


def get_sim_data(traces, reps, filename):
    """
    Description: serialize the trajectories from the output LM file into a 3D numpy array
    """


    f=h5py.File(filename, "r")

    for r in range(reps):

        traces[:,:,r]=f['Simulations'][str(r+1).zfill(7)]['SpeciesCounts'][()]

    f.close()


    return traces


def plot_hists(fig_dir, fig_name, fig_size,
               data_list, legends, colors, xlabel, ylabel, title, bins,
               mean_median=[False, False],
               title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', range=None, 
               tick_setting=[4.0, 1.5, 5, 'out'], line_widths = [1,1.5], legend_pos='best'):
    """
    Description:
    fonts_sizes: xlable, ylabel, title, legend
    tick_setting: tick_length, tick_width, tick label fontsize, direction
    """

    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    
    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
    line_width, spine_width = line_widths

    # colors = ['limegreen', 'royalblue', 'darkorange', 'purple', 'red', 'cyan']  # Predefined color list
    
    ax = plt.gca()
    
    for i, data in enumerate(data_list):
        color = colors[i]  # Cycle through colors if needed

        # print the necessary statistics of data

        print(f"Hist {fig_path} {legends[i]}: Min {np.min(data):.2f}, Mean {np.mean(data):.2f}, Median {np.median(data):.2f}, Max {np.max(data):.2f}")

        if range is None:
            plt.hist(data, bins=bins, alpha=0.7, color=color, edgecolor='black', linewidth=line_width, label=f'{legends[i]}', histtype='stepfilled')
        else:
            plt.hist(data, bins=bins, range=range, alpha=0.7, color=color, edgecolor='black', linewidth=line_width, label=f'{legends[i]}', histtype='stepfilled')
        
        
        # Add mean, median lines per dataset
        
        if mean_median[0]:
            mean = np.nanmean(data)
            plt.axvline(mean, color=color, linestyle='solid', linewidth=1.5*line_width)

            ax.text(mean, 0, f'{mean:.2f}', color=color, rotation=-45,
            ha='left', va='top', fontsize=tick_setting[2],
            transform=ax.get_xaxis_transform())  # align with x-axis

        if mean_median[1]:
            median = np.nanmedian(data)
            plt.axvline(median, color=color, linestyle='dotted', linewidth=1.5*line_width)

            ax.text(median, 0, f'{median:.2f}', color=color, rotation=45,
            ha='right', va='top', fontsize=tick_setting[2],
            transform=ax.get_xaxis_transform())  # align with x-axis

    xlabel = xlabel.replace('_', '\_')
    ax.set_xlabel(r'{0}'.format(xlabel), fontsize=xlabel_fontsize, labelpad=1.5)
    
    ylabel = ylabel.replace('_', '\_')
    ax.set_ylabel(r'{0}'.format(ylabel), fontsize=ylabel_fontsize, labelpad=1.5)
    
    if title_set:
        title = title.replace('_', '\_')
        ax.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)
    
    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    ax.tick_params(labelsize=tick_setting[2], length=tick_length, width=tick_width, direction=tick_setting[3],
                left=True, right=False, bottom=True, top=False, which='major')
    
    ax.spines['left'].set_linewidth(spine_width)
    ax.spines['bottom'].set_linewidth(spine_width)
    ax.spines['right'].set_linewidth(spine_width)
    ax.spines['top'].set_linewidth(spine_width)
    
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=legend_fontsize, loc=legend_pos, frameon=False)

    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()
    plt.close()

    return None


def plot_time_ranges(fig_dir, fig_name, fig_size,
               time, data_list, legends, colors, xlabel, ylabel, title,
               percentile=[10,90], plot_avg=True, plot_range=True, xlimit=[0,100],
               title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 'out'], line_widths = [1,1.5], legend_pos='best',
               linestyles=None, ylimit=None):
    """
    Description:
        Plot the 
    Input:
        fonts_sizes: xlable, ylabel, title, legend
        tick_setting: tick_length, tick_width, tick label fontsize, direction
    """

    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    
    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
    line_width, spine_width = line_widths

    ax = plt.gca()
    # ax.set_facecolor('white')

    ax.set_xlim(xlimit[0], xlimit[1])

    if ylimit != None:
        ax.set_ylim(ylimit[0], ylimit[1])
        
    if linestyles == None:
        linestyles = len(legends)*['-']

    for y, legend, color, ls in zip(data_list, legends, colors, linestyles):
        # print(y.shape)
        # print(f"{legend}: any all-NaN in y? {np.isnan(y).all(axis=1).any()}")
        if plot_avg:
            mean_y = np.nanmean(y, axis=1)
            ax.plot(time, mean_y, alpha=0.75, linewidth=line_width, color=color, linestyle=ls, label=f"{legend}")
            # print(mean_y, legend)

        if plot_range:
            lower_bound = np.percentile(y, percentile[0], axis=1)
            upper_bound = np.percentile(y, percentile[1], axis=1)
            # ax.fill_between(self.t / 60, lower_bound, upper_bound, color=color, alpha=0.3, label=f'{range[0]}th-{range[1]}th Percentile ({label})')

            ax.fill_between(time, lower_bound, upper_bound, color=color, alpha=0.3)

 
    xlabel = xlabel.replace('_', '\_')
    ax.set_xlabel(r'{0}'.format(xlabel), fontsize=xlabel_fontsize, labelpad=1.5)
    
    ylabel = ylabel.replace('_', '\_')
    ax.set_ylabel(r'{0}'.format(ylabel), fontsize=ylabel_fontsize, labelpad=1.5)
    
    if title_set:
        title = title.replace('_', '\_')
        ax.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)
    
    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    ax.tick_params(labelsize=tick_setting[2], length=tick_length, width=tick_width, direction=tick_setting[3],
                left=True, right=False, bottom=True, top=False, which='major')
    
    ax.spines['left'].set_linewidth(spine_width)
    ax.spines['bottom'].set_linewidth(spine_width)
    ax.spines['right'].set_linewidth(spine_width)
    ax.spines['top'].set_linewidth(spine_width)
    
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=legend_fontsize, loc=legend_pos, frameon=False)


    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()
    plt.close()

    return None

def plot_doubling(fig_dir, fig_name, fig_size,
               pkl, DNA_content,
               ylabel, title, reps, colors, legends, xlimit=[0,100],
               title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 'out'], line_widths = [1,1.5],
               legend_pos='best'):
    """
    Plot the scaled SA, Volume, and Chromosome
    """
    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    ax = plt.gca()

    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
    line_width, spine_width = line_widths

    reps = [rep -1 for rep in reps]

    vol = pkl.volumes[:,reps] # times by reps
    vol_S = vol/vol[0,:]

    SA = pkl.get_specie_trace('SA_nm2')[:,reps]
    SA_S = SA/SA[0,:]

    DNA_content = DNA_content[:,reps]
    DNA_S = DNA_content/DNA_content[0,:]

    ax.set_xlim(xlimit[0], xlimit[1])
    ax.autoscale(False) # Fix x limits

    y_low = 0.9; y_high = 2.2
    ax.set_ylim(y_low,y_high)


    quantities = [vol_S, SA_S, DNA_S]
        
    # cmap, c_space = pkl.choose_colormap(len(quantities))
    
    temp_handles = []

    for ith, quant in enumerate(quantities):
        
        mean_y = np.nanmean(quant, axis=1)
        low_y = np.percentile(quant, 10, axis=1)
        high_y = np.percentile(quant, 90, axis=1)

        # Mask out regions where y > 2
        mean_y[mean_y > 2] = np.nan
        low_y[low_y > 2] = 2
        high_y[high_y > 2] = 2

        p, = ax.plot(pkl.t/60, mean_y,
                        alpha=0.75,
                        linewidth=line_width,
                        color=colors[ith],
                        label=legends[ith])
        
        temp_handles.append(p)

        ax.fill_between(pkl.t/60, low_y, high_y,
                        color=colors[ith], alpha=0.5)
        
        # plot the ranges of doubling times on time axis
        doubling_times, not_doubled_reps = diagnosis.get_doubling_time(pkl, legends[ith], quant, reps, print_flag=False)
        
        median_doubling_time = int(np.nanmedian(doubling_times)/60)

        print(f"Doubling of {legends[ith]}: Min {np.min(doubling_times)/60:.2f} Mean {np.mean(doubling_times)/60:.2f} Median {np.median(doubling_times)/60:.2f} Max {np.max(doubling_times)/60:.2f}")

        # ax.axvline(x=median_doubling_time, color=colors[ith],
        #             linewidth=0.5, linestyle='--')

        # ax.text(median_doubling_time, -0.05, f'{median_doubling_time:.2f}', color=colors[ith], rotation=45,
        #     ha='right', va='top', fontsize=tick_setting[2], weight='heavy',
        #     transform=ax.get_xaxis_transform())  # align with x-axis

        # ax.axvline(x=max(doubling_times)/60, color=cmap(c_space[ith]),
        #            linewidth=0.5, linestyle='--')


    ax.legend(handles=temp_handles,fontsize=legend_fontsize, loc=legend_pos, frameon=False)

    if title_set:
        title = title.replace('_', '\_')
        ax.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)

    ax.set_xlabel('Time [Min]',
                    fontsize=xlabel_fontsize,
                    labelpad=1.5)
    
    ylabel = ylabel.replace('_', '\_')
    ax.set_ylabel(ylabel,fontsize=ylabel_fontsize,
                    labelpad=1.5)
    
    ax.set_xlim(xlimit[0], xlimit[1])
    ax.autoscale(False) # Fix x limits

    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    ax.tick_params(labelsize=tick_setting[2], length=tick_length, width=tick_width, direction=tick_setting[3],
                left=True, right=False, bottom=True, top=False, which='major')
    
    ax.spines['left'].set_linewidth(spine_width)
    ax.spines['bottom'].set_linewidth(spine_width)
    ax.spines['right'].set_linewidth(spine_width)
    ax.spines['top'].set_linewidth(spine_width)

    print('plot_doubling',ax.get_xlim())
    
    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()
    plt.close()

    return None

def plot_time_dualAxes(fig_dir, fig_name, fig_size,
               time, xlabel, title, percentile,
               left_data, left_legends, left_colors, left_ylabel, left_plots, left_ylabel_color,
               right_data, right_legends, right_colors, right_ylabel, right_plots, right_ylabel_color,
               xlimit=[0,100], title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 'out'], line_widths = [1,1.5],legend_pos='best'):

    fig_path = fig_dir + fig_name + extension
    fig, ax1 = plt.subplots(figsize=(fig_size[0], fig_size[1]))
    
    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
    line_width, spine_width = line_widths
    
    ax1.set_xlim(xlimit[0], xlimit[1])
    
    xlabel = xlabel.replace('_', '\_')
    ax1.set_xlabel(r'{0}'.format(xlabel), fontsize=xlabel_fontsize, labelpad=1.5)
    
    left_ylabel = left_ylabel.replace('_', '\_')
    ax1.set_ylabel(r'{0}'.format(left_ylabel), 
                   fontsize=ylabel_fontsize, color=left_ylabel_color,
                   labelpad=1.5)
    
    if title_set:
        title = title.replace('_', '\_')
        ax1.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)

    for y, legend, color, left_plot in zip(left_data, left_legends, left_colors, left_plots):
        # print(y.shape)
        # print(f"{legend}: any all-NaN in y? {np.isnan(y).all(axis=1).any()}")
        if left_plot == 'single':
            ax1.plot(time, y, alpha=1, linewidth=line_width, color=color, label=f"{legend}")

        elif left_plot == 'range':
            lower_bound = np.percentile(y, percentile[0], axis=1)
            upper_bound = np.percentile(y, percentile[1], axis=1)
            ax1.fill_between(time, lower_bound, upper_bound, color=color, alpha=0.3)

        elif left_plot == 'range_avg':
            mean_y = np.nanmean(y, axis=1)
            ax1.plot(time, mean_y, alpha=1, linewidth=line_width, color=color, label=f"{legend}")
            lower_bound = np.percentile(y, percentile[0], axis=1)
            upper_bound = np.percentile(y, percentile[1], axis=1)
            ax1.fill_between(time, lower_bound, upper_bound, color=color, alpha=0.3)

        else:
            print('Plot Method Not Matched')

    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    ax1.tick_params(labelsize=tick_setting[2], length=tick_length, width=tick_width, direction=tick_setting[3],
                left=True, right=False, bottom=True, top=False, which='major')
    
    ax1.spines['left'].set_linewidth(spine_width)
    ax1.spines['bottom'].set_linewidth(spine_width)
    ax1.spines['right'].set_linewidth(spine_width)
    ax1.spines['top'].set_linewidth(spine_width)
    
    ax2 = ax1.twinx()
    
    ax2.set_ylabel(r''+right_ylabel.replace('_','\_'),
                      fontsize=ylabel_fontsize, color=right_ylabel_color,
                      labelpad=1.5)
    
    for y, legend, color, right_plot in zip(right_data, right_legends, right_colors, right_plots):
        # print(y.shape)
        # print(f"{legend}: any all-NaN in y? {np.isnan(y).all(axis=1).any()}")
        if right_plot == 'single':
            ax2.plot(time, y, alpha=1, linewidth=line_width, color=color, label=f"{legend}")

        elif right_plot == 'range':
            lower_bound = np.percentile(y, percentile[0], axis=1)
            upper_bound = np.percentile(y, percentile[1], axis=1)
            ax2.fill_between(time, lower_bound, upper_bound, color=color, alpha=0.3)

        elif right_plot == 'range_avg':
            mean_y = np.nanmean(y, axis=1)
            ax2.plot(time, mean_y, alpha=1, linewidth=line_width, color=color, label=f"{legend}")
            lower_bound = np.percentile(y, percentile[0], axis=1)
            upper_bound = np.percentile(y, percentile[1], axis=1)
            ax2.fill_between(time, lower_bound, upper_bound, color=color, alpha=0.3)

        else:
            print('Plot Method Not Matched')

    ax2.tick_params(labelsize=tick_setting[2], length=tick_length, width=tick_width, direction=tick_setting[3],
                left=False, right=True, bottom=False, top=False, which='major')
    
    # Merge legends
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()

    if handles1+handles2:
        # Set the legend with a higher zorder
        legend = ax2.legend(handles1 + handles2, labels1 + labels2, 
                            loc=legend_pos, 
                            fontsize=legend_fontsize)
        legend.set_zorder(20)  # Ensure legend is on top


    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()
    plt.close()

    return None
