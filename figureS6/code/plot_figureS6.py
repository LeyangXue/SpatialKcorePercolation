# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 23:01:23 2024

@author: Leyang Xue
"""
import sys
import matplotlib.pyplot as plt
import numpy as np

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def PlotPictureAx(net_matrix, ax, step, title, xlabel, ylabel):
    '''
    Plot the network matrix on a given axis.

    Parameters
    ----------
    net_matrix : numpy.ndarray
        The network matrix to be plotted.
    ax : matplotlib.axes.Axes
        The axis object to plot on.
    step : int
        The step or timestamp associated with the network matrix.
    title : str
        The title of the plot.
    xlabel : str
        The label for the x-axis.
    ylabel : str
        The label for the y-axis.

    Returns
    -------
    None.

    '''    
    cmp = plt.get_cmap('viridis')
    ax.imshow(net_matrix, vmin=0, vmax=1, cmap=cmp, interpolation_stage='rgba')
    ax.set_xticks([])
    ax.set_yticks([])
    ax.text(750,950, 't='+str(step), color='white', size=8)
    kp.PlotAxes(ax, xlabel, ylabel, title)
 
def GenerateCircleCordinate(diameter, origin, L):
    '''
    Generate coordinates for plotting a circle.

    Parameters
    ----------
    diameter : float
        Diameter of the circle.
    origin : tuple
        Coordinates of the origin point of the circle.
    L : int
        Size of the lattice.

    Returns
    -------
    x_left : list
        the x coordinates.
    y_large_left : list
        y coordinates for the larger arc.
    y_small_left : list
        y coordinates for the smaller arc.

    '''  
    radius = diameter/2
    x0 = origin[0]
    y0 = origin[1]
    xseg = np.arange(0,radius,5)
    y_large = []
    y_small = []
    
    for delta_x in xseg:
        
        delta_y = np.sqrt(np.power(radius,2)-np.power(delta_x,2))
        s = y0 + delta_y
        t = y0 - delta_y
        if s > L:
            y_large.append(s%L)
        else:
            y_large.append(s)
        if t < 0:
            y_small.append((t+L)%L)
        else:
            y_small.append(t)
    
        
    y_large_left = [y_large[i] for i in np.arange(len(y_large)-1,-1,-1)]
    y_large_left.extend(y_large)
    y_small_left = [y_small[i] for i in np.arange(len(y_small)-1,-1,-1)]
    y_small_left.extend(y_small)
    x_left = list(reversed(x0-xseg))
    x_left.extend(list(x0+xseg))
    
    return x_left, y_large_left, y_small_left

def PlotSizeofHole(ax,networkpath, kcore, zeta, steps, L):
    '''
    Plot the size of holes on the given axes.

    Parameters
    ----------
    ax : list
        List of axes for plotting.
    networkpath : str
        Path to the network data.
    kcore : int
        K-core parameter.
    zeta : int
        Zeta parameter.
    steps : list
        List of steps for plotting.
    L : int
        Size of the lattice.

    Returns
    -------
    None.

    '''
    
    # Define plot parameters
    ms = 1.2
    bgms = 5
    textsize = 8
    mew = 1
    color = 'white'
    mfc = None
    titles = ['(a)', '(b)', '(c)']
    
    # Iterate over steps
    for i, step in enumerate(steps):
        
        # Load the network matrix
        net_matrix = kp.load(networkpath+f'/zeta{zeta}_step{step}_net_matrix.pkl')         
    
        # Identify the origin and diameter
        [origin, diameter]= kp.LocalizeCircleOrigin(net_matrix, L)
        
        # Plot the origin and radius
        radius = diameter/2        
        ax[i].plot(origin[0], origin[1], 'o', color=color, ms=bgms, mew = mew, mfc=mfc)
        ax[i].hlines(origin[1], origin[0], origin[0]+radius, color=color, lw= 1.2, linestyles='dashed')
        if i == 0:
            ax[i].text(origin[0]+radius/10, (origin[1]-radius/10), r'$r_h$='+str(int(radius)), color=color, size=textsize)
        else:
            ax[i].text(origin[0]+radius/6, (origin[1]-radius/10), r'$r_h$='+str(int(radius)), color=color, size=textsize)
        
        # Calculate the circle and plot it
        [xs, y_larges, y_smalls]= GenerateCircleCordinate(diameter, origin, L)
        ax[i].plot(xs, y_larges, '*', color=color, ms=ms)
        ax[i].plot(xs, y_smalls, '*', color=color, ms=ms)
        
        # Plot the picture on the axis
        PlotPictureAx(net_matrix, ax[i], step, titles[i], '', '')

def PlotSolidLines(ax, x, y, star_point, end_point, slope_line_color, slope):
    '''
    Plot lines representing slopes between two points.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The subplot axes.
    x : array_like
        The x-coordinates of the points.
    y : array_like
        The y-coordinates of the points.
    star_point : int
        The index of the starting point.
    end_point : int
        The index of the ending point.
    slope_line_color : str or tuple
        The color of the lines.
    slope : float
        The slope value.

    Returns
    -------
    None.

    '''
    # Get midpoint coordinates
    x_mid = np.convolve(x[star_point:end_point], [0.5, 0.5], mode='valid')
    x_tri = np.vstack((x_mid, x_mid + 80))
    y_tri = np.interp(x_tri, x, y)
    
    # Plot the lines
    ax.plot(x_tri+30, np.tile(y_tri[0, :], [2, 1]), color = slope_line_color)  # Horizontal line
    ax.plot(np.tile(x_tri[1, :], [2, 1])+30, y_tri, color = slope_line_color)  # Vertical line
    ax.text(x_tri[1, 0]+40, y_tri[:, 0][0]+20, f'speed = {round(slope,2)}', color = slope_line_color)#-\frac{1}{\zeta}
    
def PlotRh2time(networkpath, ax4, kcore, zetas):
    '''
    Plot hole propagation over time.

    Parameters
    ----------
    networkpath : str
        The path to the network data.
    ax4 : matplotlib.axes.Axes
        The subplot axes.
    kcore : int
        The kcore value.
    zetas : list
        A list of zeta values.

    Returns
    -------
    None.

    '''
    speed_zeta = kp.load(networkpath+'/hole_propogation/hole_speed.pkl')

    # Set the color
    color = plt.get_cmap('Paired')
    slope_line_color = plt.get_cmap("tab20c")(17)
    ms = 6
    mew = 0.8
    mfc = 'None'
    
    for i, zeta in enumerate(zetas):
        # Load the hole propagation over time
        rh_times = kp.load(networkpath+f'/hole_propogation/kcore{kcore}_zeta{zeta}_simu1_rh_time_dict.pkl')
        x = np.array(list(rh_times.keys()))
        y = np.array(list(rh_times.values()))
                
        # Plot the data
        ax4.plot(x[::5], y[::5], 'o-', mfc = mfc, color=color(2*i+1), ms = ms, mew =mew, label='$\zeta$='+str(zeta), lw=1.5)

        # Plot the hole propagation speed
        if zeta == 8:
            star_step = 150
            end_step = 450
            mid = int((star_step+end_step)/2)
            star_point = mid - 20
            end_point = mid -18
            PlotSolidLines(ax4, x, y, star_point, end_point, slope_line_color, round(speed_zeta[zeta][0],2))
            
    # Set y-axis ticks and labels
    ax4.set_yticks(np.arange(0,501,100))
    ax4.set_yticklabels(np.arange(0,501,100))
    
    # Set axis labels and title
    kp.PlotAxes(ax4, 'Step t', r'$r_h$', '(d)')
    # Add legend
    ax4.legend(loc='upper left',framealpha=0, fontsize=8.5, ncol=1)
    
def PlotSpeed(ax, networkpath):
    '''
    Plot the relationship between zeta values and hole propagation speeds.

    Parameters
    ----------
    ax5 : matplotlib.axes.Axes
        The subplot axes.
    networkpath : str
        The path to the network data.

    Returns
    -------
    None.

    '''
    speed_zeta = kp.load(networkpath+'/hole_propogation/hole_speed.pkl')
    
    slope_line_color = plt.get_cmap("tab20c")(17)
    mfc = 'None'
    
    # Extract zeta values and speeds
    x = np.array(list(speed_zeta.keys()))
    y = np.array([speed_zeta[zeta][0] for zeta in x])
    
    # Perform linear regression to find slope and intercept
    coefficients = np.polyfit(x, y, 1)
    slope, intercept = coefficients
    
    # Plot the data points and the linear fit line
    ax.plot(x, y, 'o', mfc=mfc, color=slope_line_color)
    ax.plot(x, np.polyval(coefficients, x), '-', color=slope_line_color,
        label=f'speed = {slope:.2f} $\zeta$ + {intercept:.2f}')
    
    # Add labels and legend
    kp.PlotAxes(ax, r'$\zeta$', 'speed', '(e)', mode=True)
    
def PlotRh(networkpath, figurepath, kcore, zetas):
    '''
    Plot the figures related to hole propagation analysis.

    Parameters
    ----------
    networkpath : str
        Path to the network data.
    figurepath : str
        Path to save the figures.
    kcore : int
        Value of the k-core.
    zetas : list
        List of values of zeta.

    Returns
    -------
    None.

    '''
    
    fig = plt.figure(figsize=(6, 6.5), constrained_layout = True)
    gs1 = fig.add_gridspec(nrows=3, ncols=3, hspace=0, wspace=0)#width_ratios=[1, 1, 1], height_ratios=[2, 2, 2], hspace=0.0, wspace=0.0
    ax1 = fig.add_subplot(gs1[0:1, 0:1])
    ax2 = fig.add_subplot(gs1[0:1, 1:2])
    ax3 = fig.add_subplot(gs1[0:1, 2:3])
    ax4 = fig.add_subplot(gs1[1:2, 0:3])
    ax5 = fig.add_subplot(gs1[2, 0:3])
    axes = [ax1,ax2,ax3]
    
    #plot the subfig1 
    L = 1000
    steps = [150, 240, 300]
    zeta = 10
    PlotSizeofHole(axes,networkpath, kcore, zeta, steps, L)

    #plot the subfig2
    PlotRh2time(networkpath, ax4, kcore, zetas)    
    
    #plot the subfig3
    PlotSpeed(ax5,networkpath)
    
    plt.savefig(figurepath+'/FigS6.png', dpi = 1000)
    plt.savefig(figurepath+'/FigS6.pdf')
    plt.savefig(figurepath+'/FigS6.eps')

if __name__ == "__main__":
    
    networkpath = '../kcorePercolation/figureS6/network'
    figurepath = '../kcorePercolation/figureS6/figure'
    
    kcore = 5 
    zetas = [8,10,12,14,20]
    PlotRh(networkpath, figurepath, kcore, zetas)
    