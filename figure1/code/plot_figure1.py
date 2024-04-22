# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 19:44:46 2024

@author: LeyangXue
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def PlotAxes(ax, xlabel, ylabel, title, mode=False):
    '''
    Decorate the axes
    
    Parameters
    ----------
    ax : axes
        Axes object.
    xlabel : str
        Label for the x-axis.
    ylabel : str
        Label for the y-axis.
    title : str
        Title of the plot.
    mode : bool, optional
        Whether to show the legend. Default is False.
    
    Returns
    -------
    None
    
    '''
    fontsize = 12
    font_label = {'family': "Arial", 'size':fontsize}
    
    n_legend = 10
    ax.set_xlabel(xlabel,  fontdict=font_label)
    ax.set_ylabel(ylabel, fontdict=font_label)
    ax.set_title(title, loc='left', fontdict={'family': "Arial", 'size':fontsize})
    ax.tick_params(direction='out', which='both', length=4, width=1, pad=1, labelsize=n_legend)

    # ax.minorticks_on()
    if mode == True:
        ax.legend(loc='best', framealpha=0, fontsize=n_legend)

def GenerateNetwork(L, avg_k, zeta, networkpath, filename):
    '''
    Generate a spatial network.

    Parameters
    ----------
    L : int
        Length of the lattice.
    avg_k : float
        Average degree.
    zeta : float
        Characteristic length.
    networkpath : str
        Path to save the network.
    filename : str
        Name of the file to save the network.

    Returns
    -------
    G : networkx.Graph
        Generated network.

    '''
    N = L * L  # the number of nodes 
    M = int((N * avg_k) / 2)  # the number of edges
    range_coordinate = np.arange(0, L, 1)  # range of coordinates
    
    # Generate the sequence of links with length L             
    links = kp.expvariate_sequence(zeta, M, L)

    # Create the network 
    G = kp.NetworkCreate(range_coordinate, links, networkpath, filename)
    
    return G

def kcore_percolation_LCC_nodes(G, rm_nodes, k):
    '''
    Perform k-core percolation and return the number of stages and LCC nodes dictionary.

    Parameters
    ----------
    G : networkx.Graph
        Input graph.
    rm_nodes : list
        Nodes to be removed.
    k : int
        K-core value.

    Returns
    -------
    NOI : int
        Number of stages.
    lcc_step_nodes_dict : dict
        Dictionary containing LCC nodes for each stage.

    '''
    Gr = G.copy()
    Gr.remove_nodes_from(rm_nodes)
    degree_dict = {node: degree for node, degree in Gr.degree()}

    NOI = 0
    lcc_step_nodes_dict = {}
    lcc_step_nodes_dict[NOI] = sorted(nx.connected_components(Gr), key=len, reverse=True)

    while len(degree_dict) != 0 and min(degree_dict.values()) < k:
        nodes = list(Gr.nodes())
        for node in nodes:
            if degree_dict[node] < k:
                Gr.remove_node(node)
        
        degree_dict = {node: degree for node, degree in Gr.degree()}
        NOI += 1
        lcc_step_nodes_dict[NOI] = sorted(nx.connected_components(Gr), key=len, reverse=True)
    
    return NOI, lcc_step_nodes_dict
    
def PlotCurve(ax, x1, x2, y1, y2, r, linecolor, lw, dx, n, alpha=1):
    '''
    Plot a curve with specified characteristics.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axes object for plotting.
    x1 : float
        X-coordinate of the starting point.
    x2 : float
        X-coordinate of the ending point.
    y1 : float
        Y-coordinate of the starting point.
    y2 : float
        Y-coordinate of the ending point.
    r : float
        Radius of the curve.
    linecolor : str
        Color of the curve.
    lw : float
        Line width of the curve.
    dx : float
        X-offset for the vertical lines.
    n : int
        Number of points to skip from the edges.
    alpha : float, optional
        Transparency of the curve. Default is 1.

    Returns
    -------
    None

    '''
    x = (x1 + x2) / 2
    y = (y1 + y2) / 2 
    deltays = np.arange(-r, r, 0.001) 
    deltaxs = np.array([np.sqrt(np.power(r, 2) - np.power(deltay, 2)) for deltay in deltays])
    deltax_x = x + deltaxs
    deltay_y = y + deltays
    ax.plot(deltax_x[n:-n], deltay_y[n:-n], '-', color=linecolor, lw=lw, alpha=alpha)
    ax.plot([x + dx, x + dx], [y1, deltay_y[0]], '-', color=linecolor, lw=lw, alpha=alpha)
    ax.plot([x + dx, x + dx], [deltay_y[-1], y2], '-', color=linecolor, lw=lw, alpha=alpha)
    
def PlotZetaModelRmnodes(G, L, ax, lcc, has_remove_nodes, titles, initialstate):
    '''
    Plot the illustration of spatial percolation.

    Parameters
    ----------
    G : networkx.Graph
        Input graph.
    L : int
        Length of the lattice.
    ax : matplotlib.axes.Axes
        Axes object for plotting.
    lcc : dict
        Dictionary containing LCC nodes for each stage.
    has_remove_nodes : list
        Nodes that have been removed.
    titles : str
        Title for the plot.
    initialstate : int
        Initial state flag (0 or 1).

    Returns
    -------
    None

    '''
    connected_nodes = lcc[0]
    mew = 1.2
    lmew = 1.5
    ms = 12
    fontsize = 16
    dx = 0.11
    delay = 24
    r = 0.32
    alpha = 1
    mec = 'black'
    mfc = '#85C8DD' 
    line_color = mfc
    if initialstate == 0:
        rm_color = '#FFA62B'
    else:
        rm_color = '#EEEEEE'
    
    ax.grid(visible=False, which='major', axis='both')
    ax.set_xlim(-1, 5)
    ax.set_ylim(-0.6, 4.6)
    ax.set_xticks(np.arange(L))
    ax.set_xticks([])
    ax.set_yticks(np.arange(L))
    ax.set_yticks([])
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    
    for edge_str in G.edges():
        edge_node1 = eval(edge_str[0])
        edge_node2 = eval(edge_str[1])
        if edge_str[0] in has_remove_nodes or edge_str[1] in has_remove_nodes:
            ax.plot([edge_node1[0], edge_node2[0]], [edge_node1[1], edge_node2[1]], '-', color=rm_color, lw=lmew, alpha=alpha)
        else:
            ax.plot([edge_node1[0], edge_node2[0]], [edge_node1[1], edge_node2[1]], '-', color=line_color, lw=lmew, ms=2)
    
    x = [3, 3]
    y = [1, 3]
    PlotCurve(ax, x[0], x[1], y[0], y[1], r=r, linecolor=line_color, lw=lmew, dx=dx, n=delay, alpha=alpha)

    for node_str in G.nodes():
        node = eval(node_str)
        if node in has_remove_nodes: 
            ax.plot(node[0], node[1], 'o', ms=ms, color=rm_color, mec=mec, mew=mew, alpha=alpha)
        elif node_str in connected_nodes:
            ax.plot(node[0], node[1], 'o', ms=ms, mfc=mfc, mec=mec, mew=mew)
        else:
            ax.plot(node[0], node[1], 'o', ms=ms, mfc=rm_color, mec=mec, mew=mew)

    PlotAxes(ax, '', '', '')
    ax.set_title(titles, loc='center', fontdict={'family': "Arial", 'size':fontsize})
    
def SpatialPercoIllust(networkpath, ax1, ax2, ax3, ax4):
    """
    Perform spatial percolation illustration and plot the network state at different time steps.
    
    Parameters:
        networkpath (str): 
            Path to the network data.
        ax1, ax2, ax3, ax4 (matplotlib.axes.Axes): 
            Axes to plot the network state at different time steps.
    """

    L = 5
    kcore = 3
    rm_nodes = ['(0,1)','(4,4)', '(2,1)', '(3,0)', '(1,3)']

    G = nx.read_edgelist(networkpath+'/toy_model_network/edgelist.txt')
    [NOI, lcc_step_nodes_dict] = kcore_percolation_LCC_nodes(G, rm_nodes, kcore)
    
    t0_deleted_nodes = []
    t0_deleted_nodes.extend(rm_nodes)
    PlotZetaModelRmnodes(G, L, ax1,  lcc_step_nodes_dict[0], t0_deleted_nodes, '$t=0$', initialstate=0)
    t0_deleted_nodes.extend(list(lcc_step_nodes_dict[0][0]-lcc_step_nodes_dict[1][0]))
    PlotZetaModelRmnodes(G, L, ax2,  lcc_step_nodes_dict[1], t0_deleted_nodes, '$t=1$', initialstate=1)
    t0_deleted_nodes.extend(list(lcc_step_nodes_dict[1][0]-lcc_step_nodes_dict[2][0]))
    PlotZetaModelRmnodes(G, L, ax3,  lcc_step_nodes_dict[2], t0_deleted_nodes, '$t=2$', initialstate=1)
    t0_deleted_nodes.extend(list(lcc_step_nodes_dict[2][0]-lcc_step_nodes_dict[3][0]))
    PlotZetaModelRmnodes(G, L, ax4,  lcc_step_nodes_dict[3], t0_deleted_nodes, '$t=3$', initialstate=1)
    ax1.set_title('(a)', loc='left',fontdict = {'family': "Arial", 'size':16})

def AveragePinfty(simulation_results_array):
    """
    Calculate the average value of the percolation parameter (pinfty) across multiple simulation results.

    Parameters:
       simulation_results_array (numpy.ndarray): Array containing simulation results.

    Returns:
       numpy.ndarray: Average pinfty values.
   """
   
    #set the array to store the results
    pinfty_avg = np.zeros(simulation_results_array.shape[1])
    threshold = 0.0001
    
    #identify the index of avg. pc
    avg_pc_indexs = [] 
    for each in simulation_results_array:
       pc_index =  np.argwhere(each>threshold)[0,0]
       avg_pc_indexs.append(pc_index)
    avg_pc_index = int(round(np.average(avg_pc_indexs),0))
    
    #it will be zero below the index of avg.pc, make a average over cases that is larger than the threshold above the index of avg.pc
    pinfty = simulation_results_array[:, avg_pc_index:]
    count_nonzero = pinfty.copy()
    count_nonzero[count_nonzero>threshold] = 1
    count_nonzero[count_nonzero!=1] = 0
    nonzero_sum = np.sum(count_nonzero, axis=0)
    pinfty_avg[avg_pc_index:] =  np.sum(pinfty, axis=0)/nonzero_sum
        
    return pinfty_avg

def makeaverage(kcore, zeta, simulation):
    '''
    Compute the average pinfty values for a given kcore and zeta parameter by aggregating simulation results.

    Parameters
    ----------
    kcore : int
        Core value of the network.
    zeta : int
        Zeta parameter value.
    simulation : int
        Simulation times 
        
    Returns
    -------
    pinfty_avg_dict : dict
        Dictionary containing average pinfty values.
        
    '''
    # Initialize list to store simulation results
    simulation_results =  []
    
    # Perform simulations for the specified number of times
    for simu in np.arange(simulation):
        # Load simulation results for the current iteration
        zeta_result = kp.load(resultpath+f'/kcore{kcore}/zeta{zeta}/kcore{kcore}_zeta{zeta}_simu{simu}_pinfty_dict.pkl')
        # Extract pinfty values and append to the list
        pinfty = [zeta_result[p] for p in sorted(zeta_result.keys())] 
        simulation_results.append(pinfty)
    
    # Convert simulation results to array
    simulation_results_array = np.array(simulation_results)
    # Compute average pinfty values
    pinfty_avg = AveragePinfty(simulation_results_array)
    # Create dictionary mapping zeta values to average pinfty values
    ps = sorted(zeta_result.keys())
    pinfty_avg_dict = {p:pinfty for p, pinfty in zip(ps, pinfty_avg)}
    
    return pinfty_avg_dict

def plotpinfty(ax, resultpath, kcore, zetas, colors, titles, PhaseTypes, ylabel, zeta_simu, simulation):
    """
    Plot pinfty values against p for different zeta values and a given kcore.

    Parameters:
      ax : Axes object
          The Axes object to draw the plot on.
      resultpath : str
          Path to the directory containing simulation results.
      kcore : int
          Core value of the network.
      zetas : list of int
          List of zeta parameter values.
      colors : list of str
          List of colors for plotting.
      titles : str
          Title for the plot.
      PhaseTypes : list of str
          List of phase types (e.g., 'CT' or 'DT').
      ylabel : str
          Label for the y-axis.
      zeta_simu : dict
          Dictionary mapping zeta values to simulation numbers.
      simulation : int
          Number of simulations to be averaged.

    Returns:
      None
    """
  
    #setting of parameter 
    ms = 6
    mew  = 0.85
    mfcc = 'None'
    
    # Plot pinfty for different zeta with a given kcore
    for j, zeta in enumerate(zetas):
        # Get simulation number for the current zeta
        simu = zeta_simu[zeta]
        
        # Load datasets for different zeta
        if PhaseTypes[j] == 'CT':   
            # Compute average pinfty values
            zeta_result = makeaverage(kcore, zeta, simulation)
        elif PhaseTypes[j] == 'DT':
            zeta_result = kp.load(resultpath+f'/kcore{kcore}/zeta{zeta}/kcore{kcore}_zeta{zeta}_simu{simu}_pinfty_dict.pkl')
        
        # Plot pinfty values
        ax.plot(zeta_result.keys(), zeta_result.values(), 'o-', ms= ms, color=colors[j], mfc=mfcc, mew = mew, label= r'$\zeta=$'+str(zeta))
    
    # Customize plot axes and legend
    n_legend = 12
    kp.PlotAxes(ax, r'$p$', ylabel, titles, mode=False) 
    ax.legend(loc='upper left', framealpha=0, fontsize=n_legend)

def PlotAxes1(ax,xlabel,ylabel, title, mode=False):
    '''
    Decorate the axes
    
    Parameters
    ----------
    ax : axes
        axes.
    xlabel : str
        set the xlabel.
    ylabel : str
        set the ylabel.
    title : str
        set the title.
    mode : bool, optional
        whether to show the legend. The default is False.
    
    Returns
    -------
    None.
    
    '''
    fontsize = 8
    font_label = {'family': "Arial", 'size':fontsize}
    
    n_legend = 5
    ax.set_xlabel(xlabel,  fontdict = font_label)
    ax.set_ylabel(ylabel, fontdict = font_label)
    ax.xaxis.labelpad = -1
    ax.yaxis.labelpad = -1

    ax.tick_params(direction='out', which='both',length =2, width=1, pad=0.45,labelsize=n_legend)
    
    #ax.minorticks_on()
    if mode == True:
        ax.legend(loc='best', framealpha=0, fontsize=4)
     
def transform2dict(resultpath, kcore, zeta):
    """
    Transform simulation results into dictionaries and compute delta p values.

    This function loads simulation results for a given kcore and zeta parameter,
    identifies the percolation threshold (pc), computes the average pc, and 
    calculates delta p values for each simulation.

    Parameters:
        resultpath : str
            Path to the directory containing simulation results.
        kcore : int
            Core value of the network.
        zeta : int
            Zeta parameter value.
            
    Returns:
        None
    """
    pc_pinfty_list = []
    pc_list = []

    for i in np.arange(3, 20):
        zeta_lcc_dict = kp.load(resultpath + f'/kcore{kcore}/zeta{zeta}/kcore{kcore}_zeta{zeta}_simu{i}_pinfty_dict.pkl')

        # Identify pc
        non_zero_items = {key: value for key, value in zeta_lcc_dict.items() if value != 0}
        min_value = min(non_zero_items.values())
        pc = [key for key, value in non_zero_items.items() if value == min_value][0]
        pc_list.append(pc)
        pc_pinfty = zeta_lcc_dict[pc]
        pc_pinfty_list.append(pc_pinfty)

    pc_avg = np.average(pc_list)
    delta_p_dict_simu = []

    for i in np.arange(3, 20):
        zeta_lcc_dict = kp.load(resultpath + f'/kcore{kcore}/zeta{zeta}/kcore{kcore}_zeta{zeta}_simu{i}_pinfty_dict.pkl')

        # Identify pc
        non_zero_items = {key: value for key, value in zeta_lcc_dict.items() if value != 0}
        min_value = min(non_zero_items.values())
        pc = [key for key, value in non_zero_items.items() if value == min_value][0]
        pc_pinfty = zeta_lcc_dict[pc]

        delta_p_dict = {}
        for p in non_zero_items.keys():
            if p > pc_avg:
                delta_p = p - pc_avg
                delta_pinfty = non_zero_items[p] - pc_pinfty
                delta_p_dict[delta_p] = delta_pinfty

        delta_p_dict_simu.append(delta_p_dict)

    kp.save(resultpath + f'/kcore{kcore}/zeta{zeta}/kcore{kcore}_zeta{zeta}_deltap_dict.pkl', delta_p_dict_simu)

def PlotBeta(ax, kcore, resultpath, zeta):
    """
    Plot the relationship between delta p and delta pinfty.

    This function calculates delta p values and loads the corresponding simulation
    results to plot the relationship between delta p and delta pinfty for a given kcore and zeta.

    Parameters:
        ax : Axes object
            The Axes object to draw the plot on.
        kcore : int
            Core value of the network.
        resultpath : str
            Path to the directory containing simulation results.
        zeta : int
            Zeta parameter value.
            
    Returns:
        None
    """
    # Calculate the delta p
    transform2dict(resultpath, kcore, zeta)
    
    # Load the datasets
    delta_p_dict_simu = kp.load(resultpath + f'/kcore{kcore}/zeta{zeta}/kcore{kcore}_zeta{zeta}_deltap_dict.pkl')
    i = 7
    
    # Plot settings
    color = plt.get_cmap('Paired')
    ms = 3.5
    mew = 0.55
    lw = 1
    inx = 8  # 12 8 2
    index = -10
    
    # Plot delta p vs. delta pinfty
    x = list(delta_p_dict_simu[inx].keys())
    y = list(delta_p_dict_simu[inx].values())
    ax.loglog(x[index:], y[index:], 'o', ms=ms, mfc='None', mew=mew, color=color(i))
    ax.loglog(x[index:], 0.72 * np.power(x, 0.5)[index:], '-', lw=lw, color=color(i), label=r'$\beta=0.5$')
    
    # Customize plot axes
    PlotAxes1(ax, r'$p-p_c$', r'$P_{\infty}(p)-P_{\infty}(p_c)$', '', mode=True)
    ax.set_yticks([0.02, 0.03, 0.04, 0.06])
    ax.set_yticklabels(['0.02', '0.03', '0.04', '0.06'])

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
if __name__ == "__main__":

    resultpath = '../kcorePercolation/figure1/result'
    figurepath = '../kcorePercolation/figure1/figure'
    networkpath = '../kcorePercolation/figure1/network'
    
    #set the basic information of spatial networks
    kcore = 5
    color = plt.get_cmap('Paired')
    color1 = plt.get_cmap('tab20c')

    #update the result for zeta = 500
    fig = plt.figure(figsize=(6, 8), constrained_layout = True)
    gs1 = fig.add_gridspec(nrows=3, ncols=2, hspace=-0.1, wspace=-0.1) #width_ratios=[1, 1, 1], height_ratios=[2, 2, 2], hspace=0.0, wspace=0.0
    ax1 = fig.add_subplot(gs1[0, 0:1])
    ax2 = fig.add_subplot(gs1[0, 1:2])
    ax3 = fig.add_subplot(gs1[1, 0:1])
    ax4 = fig.add_subplot(gs1[1, 1:2])
    ax5 = fig.add_subplot(gs1[2, 0:1])
    ax6 = fig.add_subplot(gs1[2, 1:2])
    ax5.sharey(ax6)
    
    #plot the illustration of spatial percolation in Sfig1 
    SpatialPercoIllust(networkpath, ax1, ax2, ax3, ax4)
    
    #set the left parameter
    zetas_left = [3,4,5,500]
    colors_left = [color1(2), color1(1), color1(0), color(7)]
    zeta_simu_left = {3:1, 4:1, 5:2, 500:3}
    PhaseTypes_left = ['CT', 'CT', 'CT','DT']
    ylabel_left = r'$P_{\infty}(p)$'
    
    #set the right parameter    
    zetas_right = [10,8,6]
    colors_right = [color(3),color(5), color(9)]
    zeta_simu_right = {10:1, 8:1, 6:1}
    PhaseTypes_right = ['DT','DT', 'DT', 'DT', 'DT']
    ylabel_right = ''
    ax6.set_yticklabels([''])

    #plot the left figure: Sfig2
    plotpinfty(ax5, resultpath, kcore, zetas_left, colors_left, '(b)', PhaseTypes_left, ylabel_left, zeta_simu_left, simulation=10)
    ax5.set_xlim(0.59, 0.765)

    #plot the right figure: Sfig3    
    plotpinfty(ax6, resultpath, kcore, zetas_right, colors_right, '(c)', PhaseTypes_right, ylabel_right, zeta_simu_right, simulation=10)
    ax6.set_xlim(0.66, 0.745)
    
    ax5_subax = ax5.inset_axes([0.16,0.18, 0.29, 0.27])
    PlotBeta(ax5_subax, kcore, resultpath, zeta=500)
    
    plt.savefig(figurepath + '/Fig1.png', dpi=500)
    plt.savefig(figurepath + '/Fig1.pdf')
    plt.savefig(figurepath + '/Fig1.eps')
    