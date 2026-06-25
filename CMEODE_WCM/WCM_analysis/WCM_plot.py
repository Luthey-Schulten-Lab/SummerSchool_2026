"""
Custom plot functions
"""

# 87 mm equals to 3.5 inch
# JPCB accepts 4.157 inch to 7 inch wide and max 9.167 inch long
# mm = 1/25.4
# l_w_ratio = 1/1.618
# fig_size = [87,87*l_w_ratio]

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as colormaps
from matplotlib.colors import to_rgba
import numpy as np
import WCM_diagnosis as diagnosis

def plot_hists(fig_dir, fig_name, fig_size,
               data_list, legends, colors, xlabel, ylabel, title, bins,
               mean_median=[False, False],
               title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', range=None, 
               tick_setting=[4.0, 1.5, 5, 'out'], legend_pos='best'):
    """
    Description:
    fonts_sizes: xlable, ylabel, title, legend
    tick_setting: tick_length, tick_width, tick label fontsize, direction
    """

    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    
    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes

    # colors = ['limegreen', 'royalblue', 'darkorange', 'purple', 'red', 'cyan']  # Predefined color list
    
    ax = plt.gca()
    
    for i, data in enumerate(data_list):
        color = colors[i]  # Cycle through colors if needed

        # print the necessary statistics of data

        print(f"Hist {fig_path} {legends[i]}: Min {np.min(data):.2E}, Mean {np.mean(data):.2E}, Median {np.median(data):.2E}, Max {np.max(data):.2E}")

        if range is None:
            plt.hist(data, bins=bins, alpha=0.7, color=color, edgecolor='black', linewidth=1, label=f'{legends[i]}', histtype='stepfilled')
        else:
            plt.hist(data, bins=bins, range=range, alpha=0.7, color=color, edgecolor='black', linewidth=1, label=f'{legends[i]}', histtype='stepfilled')
        
        
        # Add mean, median lines per dataset
        
        if mean_median[0]:
            mean = np.nanmean(data)
            plt.axvline(mean, color=color, linestyle='solid', linewidth=1.5)

            ax.text(mean, 0, f'{mean:.2E}', color=color, rotation=-45,
            ha='left', va='top', fontsize=tick_setting[2],
            transform=ax.get_xaxis_transform())  # align with x-axis

        if mean_median[1]:
            median = np.nanmedian(data)
            plt.axvline(median, color=color, linestyle='dotted', linewidth=1.5)

            ax.text(median, 0, f'{median:.2E}', color=color, rotation=45,
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
    
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=legend_fontsize, loc=legend_pos, frameon=False)

    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()

    return None


def plot_time_ranges(fig_dir, fig_name, fig_size,
               time, data_list, legends, colors, xlabel, ylabel, title,
               percentile=[10,90], plot_avg=True, plot_range=True, xlimit=[0,100],
               title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 'out'], legend_pos='best',
               linestyles=None, ylimit=None):
    """
    Description:
    fonts_sizes: xlable, ylabel, title, legend
    tick_setting: tick_length, tick_width, tick label fontsize, direction
    """

    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    
    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
        
    ax = plt.gca()
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
            ax.plot(time, mean_y, alpha=0.75, linewidth=1, color=color, linestyle=ls, label=f"{legend}")
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
    
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=legend_fontsize, loc=legend_pos, frameon=False)


    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()

    return None


def plot_multigens_time_ranges(fig_dir, fig_name, fig_size,
               times, xtick_time_interval, data_list, legends, colors, xlabel, ylabel, title,
               percentile=[10,90], plot_avg=True, plot_range=True, xlimit=[0,100],
               title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 'out'], legend_pos='best',
               linestyles=None, ylimit=None, grid=[True, 'major', 'gray', '--', 0.5, 0.6]):
    
    """
    Input:
    times: list of time series
    data_list: size of number of quantities by number of generation
    fonts_sizes: xlable, ylabel, title, legend
    tick_setting: tick_length, tick_width, tick label fontsize, direction
    """

    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    
    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
    
    ax = plt.gca()
    ax.set_xlim(xlimit[0], xlimit[1])

    if ylimit != None:
        ax.set_ylim(ylimit[0], ylimit[1])
        
    if linestyles == None:
        linestyles = len(legends)*['-']

    time_ends = [int(time[-1] - time[0]) for time in times]

    if len(times) != len(data_list[0]):
        raise ValueError(f"# Generations in times {len(times)} not math # generations {len(data_list[0])} in data_list")
    
    for y_multigen, legend, color, ls in zip(data_list, legends, colors, linestyles):
        # print(y.shape)
        # print(f"{legend}: any all-NaN in y? {np.isnan(y).all(axis=1).any()}")
        ax.plot([], [], color=color, label=f"{legend}") # Add a dummy point to create a legend

        if plot_avg:
            for i_gen, (time, y) in enumerate(zip(times, y_multigen)):
                forwared_time = time + sum(time_ends[0:i_gen])
                mean_y = np.nanmean(y, axis=1)

                ax.plot(forwared_time, mean_y, alpha=0.75, linewidth=1, color=color, linestyle=ls)

            # print(mean_y, legend)

        if plot_range:
            for i_gen, (time, y) in enumerate(zip(times, y_multigen)):
                forwared_time = time + sum(time_ends[0:i_gen])

                lower_bound = np.percentile(y, percentile[0], axis=1)
                upper_bound = np.percentile(y, percentile[1], axis=1)
                # ax.fill_between(self.t / 60, lower_bound, upper_bound, color=color, alpha=0.3, label=f'{range[0]}th-{range[1]}th Percentile ({label})')

                ax.fill_between(forwared_time, lower_bound, upper_bound, color=color, alpha=0.3)

 
    xlabel = xlabel.replace('_', '\_')
    ax.set_xlabel(r'{0}'.format(xlabel), fontsize=xlabel_fontsize, labelpad=1.5)
    
    ylabel = ylabel.replace('_', '\_')
    ax.set_ylabel(r'{0}'.format(ylabel), fontsize=ylabel_fontsize, labelpad=1.5)
    
    if title_set:
        title = title.replace('_', '\_')
        ax.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)
    
    def prepare_xticks(time_ends, xtick_time_inteval):
    
        xticks, xticklabels = [0], ["0"]

        for i_time, time_end in enumerate(time_ends):
            tick_position = sum(time_ends[0:i_time]) + np.arange(0, time_end, xtick_time_inteval)[1:]
            xticks.extend(tick_position)
            xticks.append(sum(time_ends[0:i_time])+time_end)

            tick_labels = [f"{_}" for _ in np.arange(0, time_end, xtick_time_inteval)[1:]]
            xticklabels.extend(tick_labels)
            xticklabels.append(f"{time_end}")

        return xticks, xticklabels
    
    x_positions, xticklabels = prepare_xticks(time_ends, xtick_time_interval)
    ax.set_xticks(x_positions)
    ax.set_xticklabels(xticklabels)
    
    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    ax.tick_params(labelsize=tick_setting[2], length=tick_length, width=tick_width, direction=tick_setting[3],
                left=True, right=False, bottom=True, top=False, which='major')
    
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=legend_fontsize, loc=legend_pos, frameon=False)

    # Set grid
    grid_flag, grid_which, grid_color, grid_ls, grid_lw, grid_alpha = grid
    if grid_flag:
        ax.grid(True, which=grid_which, color=grid_color, linestyle=grid_ls, linewidth=grid_lw, alpha=grid_alpha)

    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()


    return None



def plot_scatter(fig_dir, fig_name, fig_size,
               x, data_list, legends, colors, markers,  
               xlabel, ylabel, title, regression=False,
               title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 5, 'out'], legend_pos='best',
               xtick_labels=[]):
    """
    Description:
    fonts_sizes: xlable, ylabel, title, legend
    tick_setting: tick_length, tick_width, tick label fontsize, direction
    """

    from scipy.stats import linregress

    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    
    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
        
    ax = plt.gca()

    for y, legend, color, marker in zip(data_list, legends, colors, markers):
        ax.scatter(x,y, color=color, label=legend, marker=marker, alpha=0.7)

        if regression:
            slope, intercept, r_value, _, _ = linregress(x, y)
            regression_line = slope * x + intercept
            ax.plot(x, regression_line, color=color, label=f"Fit: y={slope:.2E}x+{intercept:.2E} \n R: {r_value:.2f}")

 
    xlabel = xlabel.replace('_', '\_')
    ax.set_xlabel(r'{0}'.format(xlabel), fontsize=xlabel_fontsize, labelpad=1.5)
    
    ylabel = ylabel.replace('_', '\_')
    ax.set_ylabel(r'{0}'.format(ylabel), fontsize=ylabel_fontsize, labelpad=1.5)
    
    if title_set:
        title = title.replace('_', '\_')
        ax.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)
    
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)

    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=legend_fontsize, loc=legend_pos, frameon=False)

    if xtick_labels != []:
        ax.set_xticks(x)
        ax.set_xticklabels(xtick_labels, rotation=45, ha="right")
        ax.set_xlim(x[0]-1, x[-1]+1)  # Expands x-axis to add space around boxes
        ax.margins(x=0.05, y=0.1)  # Adds padding around both axes   

    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    xticklabel_size = tick_setting[2]
    yticklabel_size = tick_setting[3]
    tick_direction = tick_setting[4]

    # X-axis ticks
    ax.tick_params(axis='x', labelsize=xticklabel_size, length=tick_length, width=tick_width,
                direction=tick_direction, bottom=True, top=False)

    # Y-axis ticks
    ax.tick_params(axis='y', labelsize=yticklabel_size, length=tick_length, width=tick_width,
                direction=tick_direction, left=True, right=False)
    
    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()

    return None

def plot_time_dualAxes(fig_dir, fig_name, fig_size,
               time, xlabel, title, percentile,
               left_data, left_legends, left_colors, left_ylabel, left_plots, left_ylabel_color,
               right_data, right_legends, right_colors, right_ylabel, right_plots, right_ylabel_color,
               xlimit=[0,100], title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 'out'], legend_pos='best'):

    fig_path = fig_dir + fig_name + extension
    fig, ax1 = plt.subplots(figsize=(fig_size[0], fig_size[1]))
    
    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
        
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
            ax1.plot(time, y, alpha=1, linewidth=1, color=color, label=f"{legend}")

        elif left_plot == 'range':
            lower_bound = np.percentile(y, percentile[0], axis=1)
            upper_bound = np.percentile(y, percentile[1], axis=1)
            ax1.fill_between(time, lower_bound, upper_bound, color=color, alpha=0.3)

        elif left_plot == 'range_avg':
            mean_y = np.nanmean(y, axis=1)
            ax1.plot(time, mean_y, alpha=1, linewidth=1, color=color, label=f"{legend}")
            lower_bound = np.percentile(y, percentile[0], axis=1)
            upper_bound = np.percentile(y, percentile[1], axis=1)
            ax1.fill_between(time, lower_bound, upper_bound, color=color, alpha=0.3)

        else:
            print('Plot Method Not Matched')

    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    ax1.tick_params(labelsize=tick_setting[2], length=tick_length, width=tick_width, direction=tick_setting[3],
                left=True, right=False, bottom=True, top=False, which='major')
    
    ax1.spines['left'].set_linewidth(1.5)
    ax1.spines['bottom'].set_linewidth(1.5)
    ax1.spines['right'].set_linewidth(1.5)
    ax1.spines['top'].set_linewidth(1.5)
    
    ax2 = ax1.twinx()
    
    ax2.set_ylabel(r''+right_ylabel.replace('_','\_'),
                      fontsize=ylabel_fontsize, color=right_ylabel_color,
                      labelpad=1.5)
    
    for y, legend, color, right_plot in zip(right_data, right_legends, right_colors, right_plots):
        # print(y.shape)
        # print(f"{legend}: any all-NaN in y? {np.isnan(y).all(axis=1).any()}")
        if right_plot == 'single':
            ax2.plot(time, y, alpha=1, linewidth=1, color=color, label=f"{legend}")

        elif right_plot == 'range':
            lower_bound = np.percentile(y, percentile[0], axis=1)
            upper_bound = np.percentile(y, percentile[1], axis=1)
            ax2.fill_between(time, lower_bound, upper_bound, color=color, alpha=0.3)

        elif right_plot == 'range_avg':
            mean_y = np.nanmean(y, axis=1)
            ax2.plot(time, mean_y, alpha=1, linewidth=1, color=color, label=f"{legend}")
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

    return None

def plot_box(fig_dir, fig_name, fig_size,
               data_list, xlabels, colors, ylabel, title,
               flier_setting=['x',3,'red'], 
               highlight_positions=[], highlight_setting=['o',5,'green'],
               highlight_boxes=[],
               title_set=True, fonts_sizes=[7, 8],
               extension='.png', tick_setting=[4.0, 1.5, 7, 5, 'out']):
    """
    Description:
    fonts_sizes: xlable, ylabel, title, legend
    tick_setting: tick_length, tick_width, tick label fontsize, direction
    """

    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    
    ylabel_fontsize, title_fontsize = fonts_sizes
    flier_marker, flier_size, flier_color = flier_setting

    # colors = ['limegreen', 'royalblue', 'darkorange', 'purple', 'red', 'cyan']  # Predefined color list
    
    ax = plt.gca()
    box = ax.boxplot(data_list, labels=xlabels, patch_artist=True, 
                     flierprops=dict(marker=flier_marker, markersize=flier_size, color=flier_color))
    
    ylabel = ylabel.replace('_', '\_')
    ax.set_ylabel(r'{0}'.format(ylabel), fontsize=ylabel_fontsize, labelpad=1.5)
    
    if title_set:
        title = title.replace('_', '\_')
        ax.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)

    ax.set_xticks(np.arange(len(data_list))+1)
    ax.set_xticklabels(xlabels)

    # if len(xlabels) > 10:
    #     plt.setp(ax.get_xticklabels(), rotation=60, ha="right")
    #     ax.set_xlim(-1, len(xlabels) + 2)  # Expands x-axis to add space around boxes
    #     ax.margins(x=0.05, y=0.1)          # Adds padding around both axes
    #     xticklabel_size = 5
    # else:

    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    xticklabel_size = tick_setting[2]


    for i, patch in enumerate(box['boxes']):
        if i in highlight_boxes:
            patch.set(facecolor=colors[1])  # Highlighted color (red)
        else:
            patch.set(facecolor=colors[0])

    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    yticklabel_size = tick_setting[3]
    tick_direction = tick_setting[4]

    # X-axis ticks
    ax.tick_params(axis='x', labelsize=xticklabel_size, length=tick_length, width=tick_width,
                direction=tick_direction, bottom=True, top=False)

    # Y-axis ticks
    ax.tick_params(axis='y', labelsize=yticklabel_size, length=tick_length, width=tick_width,
                direction=tick_direction, left=True, right=False)
    
    if highlight_positions != []: # Add markers to each positions
        high_marker, high_size, high_color = highlight_setting
        ax.scatter(np.arange(1, len(data_list) + 1), highlight_positions, 
                   color=high_color, marker=high_marker, s=high_size, 
                   zorder=20)

    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    

    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()


    return None

def plot_bars(fig_dir, fig_name, fig_size,
              data_dict, color_dict, text_dict, xlabels, ylabel, title,
              bars_between_labels=2, width=0.8,
              title_set=True, fonts_sizes=[7,8,6], variance=True,
              extension='.png', tick_setting=[4.0, 1.5, 7, 5, 'out'],
              spine_setting=[1.5,1.5,1.5,1.5], text_color=None):
    
    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    ax = plt.gca()

    ylabel_fontsize, title_fontsize, text_fontsize = fonts_sizes
    def keys_match_in_order(*dicts):
        keys_list = [list(d.keys()) for d in dicts]
        return all(k == keys_list[0] for k in keys_list[1:])
    
    if not keys_match_in_order(data_dict, color_dict, text_dict):
        print(f"WARNING: Input Size Not Match")

    x_positions = np.arange(len(data_dict))

    max_height = 0

    # Iterate different xlables
    for i_x, (xlabel, data) in enumerate(data_dict.items()):
        color = color_dict[xlabel]
        text = text_dict[xlabel]
        spacing = 1 / (len(data) + bars_between_labels)
        bar_width = width*spacing
        # Iterate different bars
        for j, (d, c, t) in enumerate(zip(data, color, text)):

            x_pos = x_positions[i_x] + (j-len(data)/2)*spacing  # Offset bars within each label
            if np.isscalar(d):
                y_value = d
                y_err = None
            else:
                y_value = np.mean(d)
                if variance:
                    y_err = np.std(d)  # or use np.std(d, ddof=1) for sample std
                else:
                    y_err = None

            bars = ax.bar(
                x_pos, y_value,
                width=bar_width,
                color=c,
                edgecolor='black',
                linewidth=bar_width/10,
                yerr=y_err,
                capsize=1)  # adds "caps" on error bars; adjust as needed)
            
            bar_heights = max([bar.get_height() for bar in bars])

            if bar_heights > max_height:
                max_height = bar_heights

            if text_color != None:
                c = text_color
            ax.text(x_pos, y_value*1.1, f"{t}",
                    ha='center', va='bottom', 
                    fontsize=text_fontsize, rotation=45, 
                    color=c)
            


    # Adjust x-axis for many labels
    # ax.set_xticks(x_positions)
    # ax.set_xticklabels(xlabels,
    #                     ha='center',
    #                     rotation=30,
    #                     fontsize=xlabel_fontsize)  # Adjusted tick labels
    
    ax.set_xticks(x_positions)
    ax.set_xticklabels(xlabels)
    plt.setp(ax.get_xticklabels(),
             rotation=30, ha='right')

    ylabel = ylabel.replace('_', '\_')
    ax.set_ylabel(r'{0}'.format(ylabel), 
                  fontsize=ylabel_fontsize, 
                  labelpad=1.5)
    
    ax.set_ylim(top=max_height * 1.3)  # extra headroom

    if title_set:
        title = title.replace('_', '\_')
        ax.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)
    
    ax.spines['left'].set_linewidth(spine_setting[0])
    ax.spines['bottom'].set_linewidth(spine_setting[1])
    ax.spines['right'].set_linewidth(spine_setting[2])
    ax.spines['top'].set_linewidth(spine_setting[3])

    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    xticklabel_size = tick_setting[2]
    yticklabel_size = tick_setting[3]
    tick_direction = tick_setting[4]

    # X-axis ticks
    ax.tick_params(axis='x', labelsize=xticklabel_size, length=tick_length, width=tick_width,
                direction=tick_direction, bottom=True, top=False)

    # Y-axis ticks
    ax.tick_params(axis='y', labelsize=yticklabel_size, length=tick_length, width=tick_width,
                direction=tick_direction, left=True, right=False)

    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()

    return None


def plot_heatmap(fig_dir, fig_name, fig_size,
               data, cmap,
               index_ticks, index_labels,
               column_ticks, column_labels,
               xlabel, ylabel, title, colorbar=['Value', 7],
               title_set=True, fonts_sizes=[7, 7, 8, 4],
               extension='.png', tick_setting=[1, 0.25, 5, 5, 5], annot=False):
    """
    twoDarray: the first dimesion is index and second one is columns

    Description: plot the heatmap of 2D data
    """
    import seaborn as sns

    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    ax = plt.gca()

    xlabel_fontsize, ylabel_fontsize, title_fontsize, annot_fontsize = fonts_sizes
    colorbar_text, colorbar_fontsize = colorbar
    tick_length, tick_width, xtick_fontsize, ytick_fontsize, cbartick_fontsize = tick_setting

    # df = pd.DataFrame(array, index = index, columns= columns)

    heatmap = sns.heatmap(data, annot=annot, annot_kws={"size": annot_fontsize}, 
                          cmap = cmap, linewidths=0.1, linecolor='gray', 
                            xticklabels=True, yticklabels=True)

    ax.set_xlabel(xlabel.replace('_','\_'),
            fontsize=xlabel_fontsize,
            labelpad=1.5)
    
    ax.set_ylabel(ylabel.replace('_','\_'),
            fontsize=ylabel_fontsize,
            labelpad=1.5)
    
    if title_set:
        title = title.replace('_', '\_')
        ax.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)
    
    ax.set_yticks(index_ticks)
    ax.set_yticklabels(index_labels, 
                       fontsize=ytick_fontsize, 
                       rotation=0,
                       va='center')

    ax.set_xticks(column_ticks)
    ax.set_xticklabels(column_labels, 
                       fontsize=xtick_fontsize, 
                       rotation=45, 
                       ha='right')

    ax.tick_params(length=tick_length,
                    width=tick_width)

    # Set colorbar font size
    cbar = heatmap.collections[0].colorbar  # Get the colorbar object
    cbar.ax.tick_params(labelsize=cbartick_fontsize, 
                        length=tick_length,
                        width=tick_width)  # Set colorbar tick font size
    
    # Set colorbar label font size
    cbar.set_label(colorbar_text, 
                   fontsize=colorbar_fontsize)  # Adjust label font size


    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()

    return None


def plot_balance(fig_dir, fig_name, fig_size,
               time, xlabel, title, left_ylabel, right_ylabel,
               gen_data, gen_legends, gen_colors, gen_type,
               con_data, con_legends, con_colors, con_type,
               cum_data, cum_legend, cum_color, cum_linestyle=[1.5, 'dotted'],
               xlimit=[0,100], title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 'out'], 
               legend_pos='best'):
    
    """
    Description: Plot the consumption and generation of one quantity and show the net cumulative on the second axe
    """
    fig_path = fig_dir + fig_name + extension
    fig, ax1 = plt.subplots(figsize=(fig_size[0], fig_size[1]))

    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes

    ax1.set_xlim(xlimit[0], xlimit[1])

    xlabel = xlabel.replace('_', '\_')
    ax1.set_xlabel(r'{0}'.format(xlabel), fontsize=xlabel_fontsize, labelpad=1.5)

    left_ylabel = left_ylabel.replace('_', '\_')
    ax1.set_ylabel(r'{0}'.format(left_ylabel), 
                   fontsize=ylabel_fontsize,
                   labelpad=1.5)
    if title_set:
        title = title.replace('_', '\_')
        ax1.set_title(r'{0}'.format(title), fontsize=title_fontsize, pad=4)

    if gen_type == 'stack':
        ax1.stackplot(time, gen_data, 
                      colors=gen_colors,labels=gen_legends, 
                        alpha=0.6)
    elif gen_type == 'line':
        ax1.plot(time, gen_data, 
                        color=gen_colors,label=gen_legends,
                        linewidth=1.5)
    else:
        print('WRONG Plot Type')

    if con_type == 'stack':
        ax1.stackplot(time, con_data, 
                      colors=con_colors,labels=con_legends, 
                        alpha=0.6)
    elif con_type == 'line':
        ax1.plot(time, con_data, 
                        color=con_colors,label=con_legends,
                        linewidth=1.5)
    else:
        print('WRONG Plot Type')

    tick_length = tick_setting[0]
    tick_width = tick_setting[1]
    ax1.tick_params(labelsize=tick_setting[2], length=tick_length, width=tick_width, direction=tick_setting[3],
                left=True, right=False, bottom=True, top=False, which='major')
    
    ax1.spines['left'].set_linewidth(1.5)
    ax1.spines['bottom'].set_linewidth(1.5)
    ax1.spines['right'].set_linewidth(1.5)
    ax1.spines['top'].set_linewidth(1.5)

    
    ax2 = ax1.twinx()
    
    ax2.set_ylabel(r''+right_ylabel.replace('_','\_'),
                      fontsize=7,
                      labelpad=1.5)

    ax2.plot(time, cum_data, color=cum_color,
             linewidth=cum_linestyle[0],
             linestyle=cum_linestyle[1],
             label=cum_legend)
    
    ax2.tick_params(labelsize=tick_setting[2], length=tick_length, width=tick_width, direction=tick_setting[3],
                left=False, right=True, bottom=False, top=False, which='major')
    
    # Merge legends
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()

    if handles1+handles2:
        # Set the legend with a higher zorder
        legend = ax2.legend(handles1 + handles2, labels1 + labels2, 
                            loc=legend_pos, 
                            fontsize=legend_fontsize,
                            frameon=False)
        legend.set_zorder(20)  # Ensure legend is on top


    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()

    return None


def plot_doubling(fig_dir, fig_name, fig_size,
               pkl, DNA_content,
               ylabel, title, reps, colors, legends, xlimit=[0,100],
               title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 'out'], 
               legend_pos='best'):
    """
    Plot the scaled SA, Volume, and Chromosome
    """
    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    ax = plt.gca()

    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
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

        doubling_times, not_doubled_reps = diagnosis.get_doubling_time(pkl, legends[ith], quant, reps, print_flag=False)

        p, = ax.plot(pkl.t/60, mean_y,
                        alpha=0.75,
                        linewidth=1,
                        color=colors[ith],
                        label=f"{legends[ith]} {np.mean(doubling_times)/60:.1f} Min")
        
        temp_handles.append(p)

        ax.fill_between(pkl.t/60, low_y, high_y,
                        color=colors[ith], alpha=0.5)
        
        # plot the ranges of doubling times on time axis
        doubling_times, not_doubled_reps = diagnosis.get_doubling_time(pkl, legends[ith], quant, reps, print_flag=False)
        
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

    ax.set_xlabel(r'$t$ [Min]',
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
    
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)

    print('plot_doubling',ax.get_xlim())
    
    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    
    plt.show()

    return None

def plot_time_reps(fig_dir, fig_name, fig_size,
               time, data_list, legends, colors, xlabel, ylabel, title,
               plot_avg=True, plot_reps=True, alphas = [1,0.3], xlimit=[0,100],
               title_set=True, fonts_sizes=[7, 7, 8, 6],
               extension='.png', tick_setting=[4.0, 1.5, 5, 'out'], legend_pos='best',
               linestyles=None, ylimit=None, vertical_line=[True, 0, 'red', '--', 'xtime']):
    """
    Description:
    Plot the traces of each rep in transparent and avg trace in solid
    fonts_sizes: xlable, ylabel, title, legend
    tick_setting: tick_length, tick_width, tick label fontsize, direction
    """

    fig_path = fig_dir + fig_name + extension
    fig = plt.figure(figsize=(fig_size[0], fig_size[1]))
    
    xlabel_fontsize, ylabel_fontsize, title_fontsize, legend_fontsize = fonts_sizes
    alpha_avg, alpha_reps = alphas
    ax = plt.gca()
    ax.set_xlim(xlimit[0], xlimit[1])

    if ylimit != None:
        ax.set_ylim(ylimit[0], ylimit[1])
        
    if linestyles == None:
        linestyles = len(legends)*['-']

    for y, legend, color, ls in zip(data_list, legends, colors, linestyles):
        # print(y.shape)
        # print(f"{legend}: any all-NaN in y? {np.isnan(y).all(axis=1).any()}")
        if plot_reps:
            for i_rep in range(y.shape[1]):
                ax.plot(time, y[:,i_rep], alpha=alpha_reps, linewidth=0.5, color=color, linestyle=ls)

    for y, legend, color, ls in zip(data_list, legends, colors, linestyles):
        if plot_avg:
            mean_y = np.nanmean(y, axis=1)
            ax.plot(time, mean_y, alpha=alpha_avg, linewidth=1, color=color, linestyle=ls, label=f"{legend}")
            # print(mean_y, legend)



    if vertical_line[0]:
        flag, x_pos, x_color, x_ls, x_legend = vertical_line
        plt.axvline(x_pos, color=x_color, linestyle=x_ls, linewidth=1, label=f"{x_legend}")

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
    
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.spines['right'].set_linewidth(1.5)
    ax.spines['top'].set_linewidth(1.5)
    
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(fontsize=legend_fontsize, loc=legend_pos, frameon=False)


    plt.tight_layout()
    fig.savefig(fig_path, dpi=600, transparent=False)
    plt.show()

    return None