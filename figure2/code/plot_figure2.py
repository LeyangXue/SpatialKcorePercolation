# -*- coding: utf-8 -*-
"""
Created on Sat Mar 25 21:35:23 2023

@author: Leyang Xue
"""
import sys
import matplotlib.pyplot as plt
import numpy as np

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def AverageSimulation(simulation_results_array, ps, threshold):
    '''
    Calculate the average critical point (Pc) and pinfty from simulation results.

    Parameters
    ----------
    simulation_results_array : ndarray
        Array containing pinfty results from simulations.
    ps : list
        List of probability values.
    threshold : float
        Threshold for k-core percolation.

    Returns
    -------
    pc_avg : float
        the average critical point (Pc).
    pinfty_avg : array
        pinfty values.

    '''
    
    # Initialize array to store pinfty averages
    pinfty_avg = np.zeros(simulation_results_array.shape[1])
    
    # Identify the index of average Pc
    avg_pc_indexs = [] 
    pc_list = []
    for each in simulation_results_array:
       # Find the index of the first pinfty value above the threshold
       pc_index =  np.argwhere(each>threshold)[0,0]
       avg_pc_indexs.append(pc_index)
       pc_list.append(ps[pc_index])
    
    # Calculate the average Pc index and Pc value
    avg_pc_index = int(round(np.average(avg_pc_indexs),0))
    pc_avg = round(np.average(pc_list),6)
    
    # For pinfty values above the average Pc index, calculate the average pinfty
    pinfty = simulation_results_array[:, avg_pc_index:]
    count_nonzero = pinfty.copy()
    count_nonzero[count_nonzero>threshold] = 1
    count_nonzero[count_nonzero!=1] = 0
    nonzero_sum = np.sum(count_nonzero, axis=0)
    pinfty_avg[avg_pc_index:] =  np.sum(pinfty, axis=0)/nonzero_sum
        
    return pc_avg, pinfty_avg

def makeaverage(rootpath, resultpath, kcore, zeta, simutimes, threshold):
    '''
    Calculate the average critical point (Pc) and pinfty for the given k-core and zeta.

    Parameters
    ----------
    rootpath : str
        Path to the root directory containing network data.
    resultpath : str
        Path to the directory containing result data.
    kcore : int
        The k-core value.
    zeta : int
        The zeta value..
    simutimes : int
        Number of simulations.
    threshold : float
        Threshold for k-core percolation.

    Returns
    -------
    pc_avg : float
        the average critical point (Pc).
    pinfty_avg_dict : dict
        pinfty values.

    '''
    simulation_results =  []
    for simu in np.arange(simutimes):
        # Load pinfty results for each simulation
        zeta_result = kp.load(f'{rootpath}/figure1/result/kcore{kcore}/zeta{zeta}/kcore{kcore}_zeta{zeta}_simu{simu}_pinfty_dict.pkl')
        pinfty = [zeta_result[p] for p in sorted(zeta_result.keys())] 
        simulation_results.append(pinfty)
        
    # Sort ps (keys) and convert simulation_results to an array
    ps = sorted(zeta_result.keys())
    simulation_results_array = np.array(simulation_results)
    
    # Calculate average Pc and pinfty
    pc_avg, pinfty_avg = AverageSimulation(simulation_results_array, ps, threshold)
    pinfty_avg_dict = {p:pinfty for p, pinfty in zip(ps, pinfty_avg)}
    
    return pc_avg, pinfty_avg_dict

def calulatePc2Zeta(rootpath, resultpath, kcore, zetas, zetacs):
    '''
    Calculate the critical point (Pc) for each zeta value and save the results.

    Parameters
    ----------
    rootpath : str
        Path to the root directory containing network and result data.
    resultpath : str
        Path to save the results.
    kcore : int
        The k-core value.
    zetas : list
        List of zeta values to analyze..
    zetacs : dict
        Dictionary of critical zeta values for each k-core.

    Returns
    -------
    None.

    '''
    simutimes = 12 # Number of simulations used to make an average for continus phase transition 
    threshold = 0.001 #Threshold of continuous phase transition for k-core percolation 
    
    avg_zeta2pc = {}
    for zeta in zetas:
        if zeta < zetacs[kcore]:
            # For zeta values less than critical zeta for the given k-core
            pc_avg, _ = makeaverage(rootpath, resultpath, kcore, zeta, simutimes, threshold)
            avg_zeta2pc[zeta] = pc_avg
        else:
            # For zeta values greater than or equal to critical zeta
            zetapc_binarysearch = kp.load(resultpath+f'/kcore{kcore}/kcore{kcore}_zeta{zeta}_binarysearchpc.pkl')
            pc_avg = np.average(zetapc_binarysearch)
            avg_zeta2pc[zeta] = pc_avg
            
    # Save the average Pc values for each zeta
    kp.save(resultpath+f'/kcore{kcore}/kcore{kcore}_binarysearchpc_avg.pkl', avg_zeta2pc)

def PlotPc2Zeta(resultpath, figurepath, kcore):
    '''
    Plot Pc as a function of zeta and mark specific points on the plot.

    Parameters
    ----------
    resultpath : str
        Path to the directory containing result data.
    figurepath : str
        Path to save the generated figures.
    kcore : int
        The k-core value.

    Returns
    -------
    None.
    
    '''
    # Load the Pc results
    zeta_threshold_pc = kp.load(resultpath+f'/kcore{kcore}/kcore{kcore}_binarysearchpc_avg.pkl')

    # Create the figure
    fig = plt.figure(figsize=(6.5, 2.5), constrained_layout = True) # constrained_layout tight_layout
    gs1 = fig.add_gridspec(nrows=1, ncols=3, hspace=3, wspace=0)
    ax = fig.add_subplot(gs1[0, 0:3])

    #parameter setting 
    color = plt.get_cmap('Paired')
    mew = 1.2
    lw = 1.5
    ms = 9
    ms_big = 12
    #bg_color = 'black'
    fontsize = 10
    font_label = {'family': "Arial", 'size':fontsize}
    n_legend = 8
    
    # Load the datasets    
    x = np.array(list(zeta_threshold_pc.keys()))
    y = np.array(list(zeta_threshold_pc.values()))
    max_index = np.argmax(y)
    
    # Plot the figure
    ax.axvline(x[max_index], color = color(9), lw =lw, ls='--')
    ax.text(x[max_index]+0.5, min(y)+(max(y)-min(y))/10, r'$\zeta_c$='+str(x[max_index]),color= color(9), size=ms_big)
    ax.plot(x, y,'o', color=color(9), ms=ms, mfc='None', mew =mew, ls='-', lw =lw)#label='threshold'
    ax.set_xscale("log")
    ax.set_xticks([3, 10, 100, 1000])
    ax.set_xticklabels([r'$3$', r'$10^1$', r'$10^2$',  r'$10^3$'])
    
    # Mark A, B, C points
    zeta_values = [4, 10, 100]
    markers = ['(a)', '(d)', '(b)']
    for zeta, marker in zip(zeta_values, markers):
        ax.plot(zeta, zeta_threshold_pc[zeta], 's', color=color(5), lw=lw, ms=ms_big, mfc='None')
        #ax.text(zeta - 0.8, zeta_threshold_pc[zeta] + 0.004, marker, color=bg_color, size=ms_big)
    ax.set_ylim(0.675, 0.727)
    
    # Insert an inset axes
    zeta_c = {3:15, 4:9, 5:6, 6:4} 
    axinx = ax.inset_axes((0.75,0.45,0.2,0.48))
    axinx.plot(zeta_c.keys(), zeta_c.values(), 'o-',  lw = lw, ms=ms, mfc='None', color =color(9))
    axinx.set_ylim(2, 17)
    axinx.set_xlim(2.5,6.5)
    axinx.set_yticks([4, 8, 12, 16])
    axinx.set_yticklabels(['4', '8', '12', '16'])
    axinx.set_xticks([3, 4, 5, 6])
    axinx.set_xticklabels(['3', '4', '5', '6'])
   
    axinx.set_xlabel('k-core', fontdict = font_label)
    axinx.set_ylabel(r'$\zeta_c$', fontdict = font_label)
    axinx.tick_params(direction='out', which='both',length =2, width=0.5, pad=0.5,labelsize=n_legend)
    kp.PlotAxes(ax, r'$\zeta$', r'$p_c$','', mode=False)
    
    plt.savefig(f"{figurepath}/fig2_s1.png", dpi=1000)
    plt.savefig(f"{figurepath}/fig2_s1.pdf") 
    plt.savefig(f"{figurepath}/fig2_s1.eps") 

def shiftmatrix(zeta_net_matrix, move_cols):
    '''
    Shift the matrix to the left by removing columns from the right and appending them to the left.

    Parameters
    ----------
    zeta_net_matrix : numpy.ndarray
        The input matrix.
    move_cols : int
        Number of columns to move from the right.

    Returns
    -------
    shifted_matrix : numpy.ndarray
        Shifted matrix.

    '''
    
    # Shift the matrix to the left
    movedcols = zeta_net_matrix[:,:move_cols]
    shifted_matrix = np.concatenate((zeta_net_matrix[:,move_cols:], movedcols), axis=1)
    
    return shifted_matrix
        
def PlotPicture(networkpath, figurepath, picturename, zeta, step, title, xlabel, ylabel, figurename, remove_cols=0):
    '''
    Plot a network picture.

    Parameters
    ----------
    networkpath : str
        Path to the directory containing network data.
    picturename : str
        Name of the picture file.
    zeta : int
        The zeta value.
    step : int
        The step value.
    title : str
        Title of the plot.
    xlabel : str
        Label for the x-axis.
    ylabel : str
        Label for the y-axis.
    figurename : str
        Name of the figure file to save.
    remove_cols : int, optional
        Number of columns to move. Defaults to 0.

    Returns
    -------
    None.

    '''    
    
    # Transform the result into net_matrix
    zeta_net = kp.load(f"{networkpath}/NetID0_avgk10_zeta{zeta}_spatialNet.pkl")
    results = kp.load(f"{networkpath}/lcc_nodes/{picturename}")
    zeta_net_matrix = kp.transform2M(results, zeta_net)
    
    if zeta == 10:
        shiftedmatrix = shiftmatrix(zeta_net_matrix, remove_cols)
        kp.save(f"{networkpath}/matrix/zeta{zeta}_step{step}_net_matrix.pkl", shiftedmatrix)   
    else:
        kp.save(f"{networkpath}/matrix/zeta{zeta}_step{step}_net_matrix.pkl", zeta_net_matrix)   

    # Load the result 
    net_matrix = kp.load(f"{networkpath}/matrix/zeta{zeta}_step{step}_net_matrix.pkl")         
    
    # Plot the picture
    fig = plt.figure(figsize=(3.25, 3.25), constrained_layout = True)
    gs = fig.add_gridspec(nrows=1, ncols=1, hspace=-0.05, wspace=0)
    ax = fig.add_subplot(gs[0,0])
    
    cmp = plt.get_cmap('viridis')
    ax.imshow(net_matrix, vmin=0, vmax=1, cmap=cmp, interpolation_stage='rgba')
    ax.set_xticks([])
    ax.set_yticks([])
    kp.PlotAxes(ax, xlabel, ylabel, title)
    
    plt.savefig(f"{figurepath}/{figurename}", dpi=900)

def PlotPictureMatrix(figurepath, figurename, net_matrix, title, xlabel, ylabel):
    '''
    Plot a network matrix.

    Parameters
    ----------
    figurepath : str
        Path to save the figure.
    figurename : str
        Name of the figure file.
    net_matrix : numpy.ndarray
        Matrix to plot.
    title : str
        Title of the plot.
    xlabel : str
        Label for the x-axis.
    ylabel : str
        Label for the y-axis.

    Returns
    -------
    None.

    '''
    fig = plt.figure(figsize=(3.25, 3.25), constrained_layout = True)
    gs = fig.add_gridspec(nrows=1, ncols=1, hspace=-0.05, wspace=0)
    ax = fig.add_subplot(gs[0,0])
    
    cmp = plt.get_cmap('viridis')
    ax.imshow(net_matrix, vmin=0, vmax=1, cmap=cmp, interpolation_stage='rgba')
    ax.set_xticks([])
    ax.set_yticks([])
    kp.PlotAxes(ax, xlabel, ylabel, title)
    plt.savefig(f"{figurepath}/{figurename}", dpi=900)

def CalculateHoleSpeed(networkpath, kcore, zetas, propogate_steps):
    '''
    Calculate the speed of hole propagation for different values of zeta.

    Parameters
    ----------
    networkpath : str
        Path to the network data.
    kcore : int
        Value of the k-core.
    zetas : list
        List of values of zeta.
    propogate_steps : dict
        Dictionary containing propagation steps for each zeta.

    Returns
    -------
    speed_zeta : dict
        A dictionary containing the speed of hole propagation for each zeta.
        Each entry contains a tuple with slope, x values, and y values.
    '''
    
    speed_zeta = {}
    
    for i, zeta in enumerate(zetas):
        
        rh_times = kp.load(networkpath+f'/hole_propogation/kcore{kcore}_zeta{zeta}_simu1_rh_time_dict.pkl')
        x = np.array(list(rh_times.keys()))
        y = np.array(list(rh_times.values()))
        
        #calculate the speed
        (star_step, end_step) = propogate_steps[zeta]
        
        # Calculate the speed using linear regression
        ex = x[star_step:end_step]
        ey = y[star_step:end_step]
        [slope, beta]= np.polyfit(ex, ey, 1)
        rh_t = slope*ex +beta
        
        speed_zeta[zeta] = (slope, ex, rh_t)
    
    kp.save(networkpath+'/hole_propogation/hole_speed.pkl', speed_zeta)
    
    return speed_zeta

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
    kp.PlotAxes(ax4, 'Step t', r'$r_h$', '')
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
    
def PlotRhSpeed(figurepath, networkpath, kcore, zetas):
    '''
    Plot the hole propagation speed.

    Parameters
    ----------
    figurepath : str
        The path where the plot will be saved.
    networkpath : str
        The path to the network data.
    kcore : int
        The kcore value.
    zetas : list
        A list of zeta values.
    
    Returns
    -------
    None.

    '''
    fig = plt.figure(figsize=(6.5, 2.5), constrained_layout = True)
    gs1 = fig.add_gridspec(nrows=1, ncols=3, hspace=3, wspace=0)
    ax = fig.add_subplot(gs1[0, 0:3])
    PlotRh2time(networkpath, ax, kcore, zetas)    
    plt.savefig(figurepath+'/fig2_s3.png', dpi = 1000)
    
if __name__ == "__main__":
    
    # Set the root, result, network, and figure paths
    rootpath = '../kcorePercolation/'
    resultpath = '../kcorePercolation/figure2/result'
    networkpath = '../kcorePercolation/figure2/network'
    figurepath = '../kcorePercolation/figure2/figure'
    
    #calculate the pc as a function of zeta
    kcore = 5 
    zetas = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,50,100,500,1000]
    zetacs = {2:0, 3:15, 4:9, 5:6, 6:4}
    calulatePc2Zeta(rootpath, resultpath, kcore, zetas, zetacs)
    
    #FigS1: plot the pc as a function of zeta
    PlotPc2Zeta(resultpath, figurepath, kcore)
    
    #FigS2: plot the picuture for zeta = 4
    zeta = 4
    step = 58
    
    # Option 1: Plot the picture directly from the network matrix
    #picturename = f'kcore{kcore}_zeta{zeta}_p0.701501_simu1_step{step}_lcc_nodes.pkl'
    #PlotPicture(networkpath, figurepath, picturename, zeta, step, '', '', ' \n', f'fig2_zeta{zeta}_step{step}')
    
    # Option 2: Plot the picture using pre-saved network matrix
    zeta4_net_matrix = kp.load(networkpath+'/matrix/zeta4_net_matrix.pkl')         
    PlotPictureMatrix(figurepath, f'fig2_zeta{zeta}_step{step}', zeta4_net_matrix,'', '', ' \n')
    
    # Plot pictures for zeta = 10 at different steps
    zeta = 10
    remove_cols = 400
    steps = [140, 200, 280]
    for step in steps:
        picturename = f'kcore{kcore}_zeta{zeta}_p0.700501_simu1_step{step}_lcc_nodes.pkl'
        PlotPicture(networkpath, figurepath, picturename, zeta, step, '', '', ' \n', f'fig2_zeta{zeta}_step{step}', remove_cols)
    
    # Plot pictures for zeta = 100
    zeta = 100
    step = 199
    picturename = f'kcore{kcore}_zeta{zeta}_p0.680601_simu1_step{step}_lcc_nodes.pkl'
    PlotPicture(networkpath, figurepath, picturename, zeta, step, '', '', ' \n', f'fig2_zeta{zeta}_step{step}')
    
    # Calculate the hole propagation speed
    zetasspeed = [6,8,10,12,14,16,18,20]
    propogate_steps = {6:(200, 800),  8:(130, 450), 10:(100,300), 12:(70,200), 14:(70,190), 16:(70,180), 18:(70,150), 20:(70,150)}
    zetas_speed = CalculateHoleSpeed(networkpath, kcore, zetasspeed, propogate_steps)

    # Plot the size of holes 
    zetas = [8,10,12,14,20]
    PlotRhSpeed(figurepath, networkpath, kcore, zetas)
