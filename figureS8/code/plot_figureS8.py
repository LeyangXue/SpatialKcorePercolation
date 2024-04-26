# -*- coding: utf-8 -*-
"""
Created on Sat Sep  2 01:30:28 2023

@author: Leyang XUe
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx 

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def GetfailedGCC(G, removed_node_set):
    """
    Get the largest connected component for node removal.

    Parameters:
        G: networkx.Graph 
            The original graph.
        removed_node_set : set 
            Set of removed nodes.

    Returns:
        fgcc_nodes : set 
            Nodes in the largest connected component after node removal.
    """
    
    G_copy = G.copy()
    H = G_copy.subgraph(removed_node_set)
    fgcc_nodes = max(nx.connected_components(H), key=len)
    
    return fgcc_nodes

def transform2M(lcc,L):
    """
    Transform the node coordinates to a matrix representation.
    
    Parameters:
        lcc: set
            Set of node coordinates.
        L: int
            Size of the matrix.
    
    Returns:
        net_matrix: numpy.ndarray 
            Matrix representation of the connected component.
    """  
    
    # Draw the connected component 
    net_matrix = np.zeros((L, L))
    if len(lcc) > 0:
        for index in lcc:
            x = int(index[0])
            y = int(index[1])
            net_matrix[x,y] = 1
    
    return net_matrix

def LoadPicture(networkpath, resultpath, zeta, simulation, ps, zeta_steps):
    '''
    Load network data and process failure nodes.

    Parameters
    ----------
    networkpath : str
        Path to the network folder.
    resultpath : str
        Path to the result folder.
    zeta : int
        Zeta value.
    simulation : dict
        Dictionary of zeta values and corresponding simulation numbers.
    ps : dict
        Dictionary of zeta values and corresponding percolation thresholds.
    zeta_steps : dict
        Dictionary of zeta values and corresponding steps.

    Returns
    -------
    None.

    '''
    p = ps[zeta] 
    simu = simulation[zeta]
    steps = zeta_steps[zeta]
    
    # Load the network 
    G = kp.load(networkpath + f'/NetID0_avgk10_zeta{zeta}_spatialNet.pkl')
    nodeset = set(G.nodes())
    
    # Load the lcc-steps and fgcc nodes 
    lcc_steps_nodes = kp.load(resultpath+f'/zeta{zeta}/kcore5_zeta{zeta}_p{p}_simu{simu}_steps_lcc_nodes_dict.pkl')
        
    for i,step in enumerate(steps):
        print(f'step={step}')
        
        # Load the nodes for fgcc
        pinfty_nodes = lcc_steps_nodes[step]
        removed_node_set = nodeset - set(pinfty_nodes)
        fcc_nodes = GetfailedGCC(G, removed_node_set)
        kp.save(resultpath+f'/fgcc_nodes/zeta{zeta}_failure_nodes_step{step}.pkl', fcc_nodes)
        
def PlotAxes(ax, xlabel, ylabel, title, mode=False, log=True):
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
    fontsize = 11
    nlegend = 9
    
    font_label = {'family': "Arial", 'size': fontsize}
    if log == True:
        ax.set_yscale('log')
    
    if mode == True:
        ax.legend(loc='best', framealpha=0, fontsize=nlegend)

    ax.set_xlabel(xlabel,  fontdict=font_label)
    ax.set_ylabel(ylabel, fontdict=font_label)
    ax.set_title(title, loc='left', fontdict={
                 'family': "Arial", 'size': fontsize}, pad=2)
    ax.tick_params(direction='out', which='both', length=4,
                   width=1, pad=1, labelsize=fontsize-2)
    
def PlotSubAxes(ax, xlabel, ylabel, title, mode=False, log=True):
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
    fontsize = 9
    nlegend = 10
    
    font_label = {'family': "Arial", 'size': fontsize}
    if log == True:
        ax.set_yscale('log')
    
    if mode == True:
        ax.legend(loc='best', framealpha=0, fontsize=nlegend)

    ax.set_xlabel(xlabel,  fontdict=font_label)
    ax.set_ylabel(ylabel, fontdict=font_label)
    ax.text(380,950, title, color='white', size=fontsize)
    ax.tick_params(direction='out', which='both', length=4,
                   width=1, pad=1, labelsize=fontsize-2)
    
def shiftmatrix(zeta_net_matrix, remove_cols):
    """
    Shift the matrix to the left.

    Parameters:
        zeta_net_matrix: numpy.ndarray
            Matrix representation of the network.
        remove_cols: int 
            Number of columns to remove.

    Returns:
        numpy.ndarray: Shifted matrix.
    """
    removedcols = zeta_net_matrix[:,:remove_cols]
    shifted_matrix = np.concatenate((zeta_net_matrix[:,remove_cols:], removedcols), axis=1)
    
    return shifted_matrix

def plot_figure(networkpath, resultpath, figurepath, zeta, ps, simulation, zeta_steps, zeta_seg, titles, index, legend_mode):
    '''
    Plot figures showing percolation dynamics.

    Parameters
    ----------
    networkpath : str
        Path to the network folder.
    resultpath : str
        Path to the result folder.
    figurepath : str
        Path to save the figure.
    zeta : int
        Zeta value.
    ps : dict
        Dictionary of zeta values and corresponding percolation thresholds.
    simulation : dict
        Dictionary of zeta values and corresponding simulation numbers.
    zeta_steps : dict
        Dictionary of zeta values and corresponding steps.
    zeta_seg : dict
        Dictionary of zeta values and corresponding segment information.
    titles : list
        List of subplot titles.
    index: int
        Index of titles
    legend_mode: bool
        Whether to display legend in the plot.
        
    Returns
    -------
    None.

    '''
    
    # Set figure parameters
    cmp = plt.get_cmap('viridis_r')
    color = plt.get_cmap('Paired')
    color1 = plt.get_cmap('viridis_r')(0)
    mfc = 'None'
    ms = 7
    p = ps[zeta]
    simu = simulation[zeta]

    # Plot the picture at each step for a given zeta
    fig = plt.figure(figsize=(4.125, 5.49), tight_layout=True)
    gs1 = fig.add_gridspec(nrows=9, ncols=3, hspace=0.0)
    ax1 = fig.add_subplot(gs1[0:2, 0:3])
    
    # Plot Sfig1
    pinfty_steps = kp.load(resultpath + f'/zeta{zeta}/kcore5_zeta{zeta}_p{p}_simu{simu}_steps_lcc_num_dict.pkl')
    fgcc_steps = kp.load(resultpath + f'/zeta{zeta}/kcore5_zeta{zeta}_p{p}_simu{simu}_steps_fgcc_nodes_dict.pkl')
    
    x_fgcc = np.array(list(fgcc_steps.keys()))
    y_fgcc = np.array(list(fgcc_steps.values()))
    x_pinfty = np.array(list(pinfty_steps.keys()))
    y_pinfty = np.array(list(pinfty_steps.values()))

    inx = np.arange(0, len(x_fgcc), 5)
    ax1.plot(x_fgcc[inx], y_fgcc[inx], 'o-', color=color(9), ms=ms, mfc=mfc, label=r'$P^f_{\infty}(t)$')
    ax1.plot(x_pinfty[inx], y_pinfty[inx], 'o-', color=color1, ms=ms, mfc=mfc, label=r'$P_{\infty}(t)$')

    xtick = np.arange(0, zeta_seg[zeta][0], zeta_seg[zeta][1])
    ax1.set_xticks(xtick)
    ax1.set_xticklabels(xtick)
    PlotAxes(ax1, 'Step t', 'ratio of nodes', f'({titles[index]})  '+r'$\zeta=$'+f'{zeta}', mode=legend_mode, log=False)

    # The remaining subplots will be placed in subsequent rows
    gs2 = fig.add_gridspec(nrows=9, ncols=3, hspace=0.00, wspace=-0.008)#, height_ratios=[2, 2, 2],
    ax2 = fig.add_subplot(gs2[3:5, 0])
    ax3 = fig.add_subplot(gs2[3:5, 1])
    ax4 = fig.add_subplot(gs2[3:5, 2])
    ax5 = fig.add_subplot(gs2[5:7, 0])
    ax6 = fig.add_subplot(gs2[5:7, 1])
    ax7 = fig.add_subplot(gs2[5:7, 2])
    ax8 = fig.add_subplot(gs2[7:9, 0])
    ax9 = fig.add_subplot(gs2[7:9, 1])
    ax10 = fig.add_subplot(gs2[7:9, 2])
    axes = [ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9, ax10]
    
    # Set the null for ticklabels
    for ax in axes:
        ax.set_xticklabels("")
        ax.set_yticklabels("")
        ax.set_xticks([])
        ax.set_yticks([])
        
    # Plot Sfig2-9
    for i, step in enumerate(zeta_steps[zeta]):
        # Load the nodes for fcc and pinfty
        fcc_nodes = kp.load(resultpath + f'/fgcc_nodes/zeta{zeta}_failure_nodes_step{step}.pkl')

        # Transform the nodes into a matrix
        fcc_matrix = transform2M(fcc_nodes, L)
        if zeta == 10:
            remove_cols = 400
            shiftedmatrix = shiftmatrix(fcc_matrix, remove_cols)
            axes[i].imshow(shiftedmatrix, vmin=0, vmax=1, cmap=cmp)
        else:
            axes[i].imshow(fcc_matrix, vmin=0, vmax=1, cmap=cmp)
        PlotSubAxes(axes[i], '', '', f'({titles[index]}{i+1}) t={step}', mode=False, log=False)
        
    plt.savefig(figurepath+f'/zeta{zeta}.png', dpi=1000)
    plt.savefig(figurepath+f'/zeta{zeta}.pdf')
    plt.savefig(figurepath+f'/zeta{zeta}.eps')

if __name__ == "__main__":

    # Set the result, network, and figure paths
    figurepath  = '../kcorePercolation/figureS8/figure'
    resultpath = '../kcorePercolation/figureS8/result'
    networkpath = '../kcorePercolation/figureS8/network'
    
    # Set the parameter
    L = 1000
    zetas = [4, 10, 100]
    ps = {4:0.701501, 10:0.700501, 100:0.679701}
    simulation = {4:1, 10:1, 100:2}
    zeta_steps = {4:[1,10,20,40,50,60,70,80,320], 10:[1,45,90,135,180,225,270,315,360], 100:[1,35,70,105,140,175,210,245,280]}
    zeta_seg = {4:(328,40), 10:(365,45), 100:(285,35)}
    titles = ['a', 'b', 'c']
    
    # load the picture for different zeta
    for zeta in zetas: 
        LoadPicture(networkpath, resultpath, zeta, simulation, ps, zeta_steps)
    
    # Plot the zeta = 4
    zeta = 4 
    index = 0
    plot_figure(networkpath, resultpath, figurepath, zeta, ps, simulation, zeta_steps, zeta_seg, titles, index, legend_mode=False)
    
    #Plot the zeta = 10
    zeta = 10
    index = 1
    plot_figure(networkpath, resultpath, figurepath, zeta, ps, simulation, zeta_steps, zeta_seg, titles, index, legend_mode=True)

    #Plot the zeta = 100
    zeta = 100
    index = 2
    plot_figure(networkpath, resultpath, figurepath, zeta, ps, simulation, zeta_steps, zeta_seg, titles, index, legend_mode=False)
