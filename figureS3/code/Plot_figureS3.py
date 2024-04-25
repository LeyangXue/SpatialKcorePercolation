# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 18:58:53 2024

@author: Leyang Xue
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import math 

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def setting_ps():
    """
    Define the ps (percolation probability) values for different system sizes.
    
    Returns
    -------
    dict
        Dictionary containing system sizes as keys and corresponding ps values as arrays.
    """
    ps_dict = {}
    ps_dict[100] = np.arange(0.90,0.94,0.001)
    ps_dict[300] = np.arange(0.90,0.94,0.001)
    ps_dict[500] = np.arange(0.90,0.94,0.001)
    ps_dict[1000] = np.arange(0.925,0.94,0.001)
    ps_dict[1500] = np.arange(0.925,0.94,0.001) # 0.925, 0.91-0.95
    ps_dict[2000] = np.arange(0.925,0.94,0.001) #0.91-0.95
    
    return ps_dict

def CalculateDist(links, L, resultpath):
    """
    Calculate the distribution of link lengths for a given network.
    
    Parameters
    ----------
    links : array_like
        Array containing the lengths of links in the network.
    L : int
        System size.
    resultpath : str
        Path to save result files.
    
    Returns
    -------
    None
    """
    
    hist, edges = np.histogram(links, bins=50)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    pdf = hist/np.sum(hist)
    pdf_reverse = pdf[::-1]
    cdf_reverse = np.cumsum(pdf_reverse)[::-1]
    
    kp.save(resultpath + f'/bin_centers_L{L}_k7.5_zeta10.pkl', bin_centers)
    kp.save(resultpath + f'/pdf_L{L}_k7.5_zeta10.pkl', pdf)
    kp.save(resultpath + f'/cdf_reverse_L{L}_k7.5_zeta10.pkl', cdf_reverse)

def CalculateDegreeDist(network, L, resultpath):
    """
    Calculate the degree distribution for a given network.
    
    Parameters
    ----------
    network : NetworkX graph
        The network for which the degree distribution is calculated.
    L : int
        System size.
    resultpath : str
        Path to save result files.
    
    Returns
    -------
    None
    """
    
    degree_seq = list(dict(network.degree()).values())
    degree_distribution = Counter(degree_seq)
    total_nodes = sum(degree_distribution.values())
    
    #calculate probability distribution
    probability_distribution = {k: degree_distribution[k] / total_nodes for k in sorted(degree_distribution.keys())}
    kp.save(resultpath + f'/degree_distribution_L{L}_k7.5_zeta10.pkl', probability_distribution)


def LoadNetLink(networkpath, resultpath, Ls):
    """
    Load network links and calculate distributions for different system sizes.
    
    Parameters
    ----------
    networkpath : str
        Path to load network files.
    resultpath : str
        Path to save result files.
    Ls : list of int
        List of system sizes.
    
    Returns
    -------
    None
    """
    for L in Ls:
        print(f'run the system size L={L}')
        
        #load the link 
        links = kp.load(networkpath+f'/Net_L{L}_avgk7.5_links_zeta10.pkl')
        CalculateDist(links, L, resultpath)
        
        #load the network 
        network = kp.load(networkpath + f'/Net_L{L}_avgk7.5_zeta10_spatialNet.pkl')
        CalculateDegreeDist(network, L, resultpath)

def IdentifyPcNOI(step, ps, pinfty):
    """
    Identify pc using the number of iterations (NOI).
    
    Parameters
    ----------
    step : array-like
        Array containing the number of occupied iterations.
    ps : array-like
        Array containing the values of p.
    pinfty : array-like
        Array containing the values of pinfty.
    
    Returns
    -------
    pc : float
        Value of pc.
    pc_pinfty : float
        Value of pinfty at pc.
    first_point : int
        Index of the point where step is maximum.
    """

    first_point = np.argmax(step)
    pc = ps[first_point]
    pc_pinfty = pinfty[first_point]
    
    return pc, pc_pinfty, first_point

def IdentifyPcThreshold(pinfty, ps):
    """
    Identify pc using a pinfty threshold.
    
    Parameters
    ----------
    pinfty : array-like
        Array containing the values of pinfty.
    ps : array-like
        Array containing the values of p.
    
    Returns
    -------
    pc : float
        Value of pc.
    pc_pinfty : float
        Value of pinfty at pc.
    first_point : int
        Index of the first point where pinfty exceeds the threshold.
    """
    # Define the threshold
    threshold = 0.05
    
    # Find the index of the first point where pinfty exceeds the threshold
    first_point = min(np.where(np.array(pinfty)>threshold)[0])
    
    # Get the corresponding values of pc and pinfty
    pc = ps[first_point]
    pc_pinfty = pinfty[first_point]
        
    return pc, pc_pinfty, first_point

def ParseResult(results, ps_dict, L):
    """
    Parse percolation results to calculate averages and standard deviations of pc and pinfty.
    
    Parameters
    ----------
    results : list
        List of percolation results, each containing kcore, pinfty, NOI, and pinfty_step.
    ps_dict : dict
        Dictionary containing values of p for each system size.
    L : int
        System size.
    
    Returns
    -------
    pc_avg : float
        Average value of pc.
    pc_std : float
        Standard deviation of pc.
    pinfty_pc_avg : float
        Average value of pinfty at pc.
    pinfty_pc_std : float
        Standard deviation of pinfty at pc.
    """
    
    # Get the values of p for the given system size
    ps = ps_dict[L]
    
    # Initialize lists to store pc, pc_pos, and pinfty
    pc_list = []
    pc_pos = []
    pinfty_list = []
        
    # Iterate over each simulation result
    for each in results:
        
        # Extract kcore, pinfty, NOI, and pinfty_step
        [kcore, pinfty, NOI, pinfty_step] = each
        pinfty_list.append(pinfty)
        
        # Identify pc using pinfty
        #[pc, pc_pinfty, first_point] = IdentifyPcNOI(NOI, ps, pinfty)
        [pc, pc_pinfty, first_point] = IdentifyPcThreshold(pinfty, ps)
        
        # Store pc and its position
        pc_list.append(pc)
        pc_pos.append(first_point)
   
    # Calculate the average position and pc
    pc_pos_avg = math.ceil(np.average(pc_pos))
    pc_avg = np.average(pc_list)
    pc_std = np.std(pc_list)
    
    # Convert pinfty_list to a numpy array
    pinfty_array = np.array(pinfty_list)
    
    # Get pinfty at the criticality
    pinfty_pc = pinfty_array[:,pc_pos_avg]
    pinfty_pc_nonzero = pinfty_pc[pinfty_pc!=0]
    pinfty_pc_avg = np.average(pinfty_pc_nonzero)
    pinfty_pc_std = np.std(pinfty_pc_nonzero)
        
    return pc_avg, pc_std, pinfty_pc_avg, pinfty_pc_std

def LoadPcPinty(resultpath, ps_dict, Ls):
    """
    Load percolation results and calculate averages and standard deviations of pc and pinfty for different system sizes.
    
    Parameters
    ----------
    resultpath : str
        Path to load percolation result files and save calculated results.
    ps_dict : dict
        Dictionary containing values of p for each system size.
    Ls : list of int
        List of system sizes.
    
    Returns
    -------
    None
    """
    pc_avg_dict = {}
    pc_std_dict = {}
    pinfty_std_dict = {}
    pinfty_avg_dict = {}
    
    for L in Ls:
        print(f'run the system size L={L}')
        #load the results 
        results = kp.load(resultpath+f'/Results_L{L}_avgk7.5_zeta10_kcore5_Percolation.pkl')
        #parse the result and return average of pc and pinfty  
        [pc_avg, pc_std, pinfty_avg, pinfty_std] = ParseResult(results, ps_dict, L)
        
        #save the result 
        pc_avg_dict[L] = pc_avg
        pc_std_dict[L] = pc_std
        pinfty_avg_dict[L] = pinfty_avg
        pinfty_std_dict[L] = pinfty_std
        
        kp.save(resultpath + '/pc_avg_k7.5_zeta10.pkl', pc_avg_dict)
        kp.save(resultpath + '/pc_std_k7.5_zeta10.pkl', pc_std_dict)
        kp.save(resultpath + '/pinfty_avg_k7.5_zeta10.pkl', pinfty_avg_dict)
        kp.save(resultpath + '/pinfty_std_k7.5_zeta10.pkl', pinfty_std_dict)

def PlotAxes(ax,xlabel,ylabel, title, mode=False):
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
    fontsize = 12
    font_label = {'family': "Arial", 'size':fontsize}
    
    n_legend = 10
    ax.set_xlabel(xlabel,  fontdict = font_label)
    ax.set_ylabel(ylabel, fontdict = font_label)
    ax.set_title(title, loc='left',fontdict = {'family': "Arial", 'size':fontsize})
    ax.tick_params(direction='out', which='both',length =4, width=1, pad=1,labelsize=n_legend)

    #ax.minorticks_on()
    if mode == True:
        ax.legend(loc='best', framealpha=0, fontsize=n_legend)

def PlotSubAxes(ax,xlabel,ylabel, title, mode=False):
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
    fontsize = 10
    font_label = {'family': "Arial", 'size':fontsize}
    
    n_legend = 8
    ax.set_xlabel(xlabel,  fontdict = font_label)
    ax.set_ylabel(ylabel, fontdict = font_label)
    ax.set_title(title, loc='left',fontdict = {'family': "Arial", 'size':fontsize})
    ax.tick_params(direction='out', which='both',length =2, width=0.5, pad=0.5,labelsize=n_legend)

    #ax.minorticks_on()
    if mode == True:
        ax.legend(loc='best', framealpha=0, fontsize=n_legend)

def PoissonDistribution(mean_value = 7.5):
    """
    Generate a Poisson distribution.
    
    Parameters
    ----------
    mean_value : float, optional
        Mean value of the Poisson distribution. The default is 7.5.
    
    Returns
    -------
    tuple
        Tuple containing x values and the corresponding Poisson distribution probabilities.
    """
    
    poisson_sequence = np.random.poisson(mean_value, size=1000000)
    poisson_distribution = np.bincount(poisson_sequence) / len(poisson_sequence)
    x = range(len(poisson_distribution))
    
    return x, poisson_distribution

def TheoreticalLine(x_data, slope):
    """
    Generate a theoretical line based on a given slope and x data.
  
    Parameters
    ----------
    x_data : array_like
        Array containing x values.
    slope : float
        Slope of the line.
  
    Returns
    -------
    array_like
        Array containing the corresponding y values of the theoretical line.
    """

    y_data  = np.exp(slope * x_data)    
    return y_data

def plot_system_size_dependency(resultpath, figurepath, Ls, slope):
    """
    Plot the system size dependency.

    Parameters
    ----------
    resultpath : str
        Path to load the result files.
    figurepath : str
        Path to save the generated figures.
    Ls : list of int
        List of system sizes.
    slope : float
        Slope value for theoretical line.

    Returns
    -------
    None
    """
    
    # Load the result and plot the figures
    colors = plt.get_cmap("Paired")
    line_color = plt.get_cmap("tab20c")(17)
    slope_line_color = plt.get_cmap("tab20c")(17)

    lw = 1.5
    mfc = 'None'
    ms = 7
    mks = ['o', 's', '^', '*','h','x', 'p'] 
    
    fig, ax = plt.subplots(2, 2, figsize=(7.5, 6), constrained_layout=True)
    
    # Plot the Poisson distribution
    [x, poisson_dist] = PoissonDistribution(mean_value=7.5)
    ax[0, 1].plot(x, poisson_dist, '-', color=line_color, lw=lw, label=r'Poisson($\lambda=7.5$)') 

    for j, L in enumerate(Ls):
        
        # Degree distribution
        degree_distribution = kp.load(resultpath + f'/degree_distribution_L{L}_k7.5_zeta10.pkl')
        
        # Link distribution
        bin_centers = kp.load(resultpath + f'/bin_centers_L{L}_k7.5_zeta10.pkl')
        pdf = kp.load(resultpath + f'/pdf_L{L}_k7.5_zeta10.pkl')
        
        # Plot the distribution of link length
        ax[0, 0].plot(bin_centers, pdf, mks[j], ms=ms, color=colors(j), mfc=mfc) 
        ax[0, 1].plot(degree_distribution.keys(), degree_distribution.values(), mks[j], color=colors(j), mfc=mfc, label=f'L={L}') 

    # Fitted lines
    y_data = TheoreticalLine(bin_centers, slope)
    x_mid = np.convolve(bin_centers[15:17], [0.5, 0.5], mode='valid')
    x_tri = np.vstack((x_mid, x_mid + 40))
    y_tri = np.interp(x_tri, bin_centers, y_data)
    
    ax[0, 0].plot(bin_centers, y_data, '-', color=line_color, lw=lw, label=f'slope = {slope}')
    ax[0, 0].plot(x_tri+5, np.tile(y_tri[0, :], [2, 1]), color=slope_line_color)  # red horizontal line
    ax[0, 0].plot(np.tile(x_tri[1, :], [2, 1])+5, y_tri, color=slope_line_color)  # red vertical line
    ax[0, 0].text(x_tri[1, 0] + 10, np.exp(np.mean(np.log(y_tri[:, 0]))), r'slope=$-0.1$')

    ax[0, 0].set_yscale('log')
    PlotAxes(ax[0, 0], 'link length, '+'$l$', r'$p(l)$', r'(a) $\langle k \rangle$ = 7.5, $\zeta=10$', mode=False)
    PlotAxes(ax[0, 1], 'degree, '+'$k$', r'$p(k)$', r'(b) $\langle k \rangle$ = 7.5, $\zeta=10$', mode=True)

    # Error bars for pc_avg and pinfty_avg
    pc_avg_dict = kp.load(resultpath + '/pc_avg_k7.5_zeta10.pkl')
    pc_std_dict = kp.load(resultpath + '/pc_std_k7.5_zeta10.pkl')
    pinfty_avg_dict = kp.load(resultpath + '/pinfty_avg_k7.5_zeta10.pkl')  
    pinfty_std_dict = kp.load(resultpath + '/pinfty_std_k7.5_zeta10.pkl')  
    
    ax[1, 0].errorbar(pc_avg_dict.keys(), pc_avg_dict.values(), yerr=list(pc_std_dict.values()), fmt='o-', capsize=3, ms=ms, mfc=mfc, color=line_color)
    ax[1, 1].errorbar(pinfty_avg_dict.keys(), pinfty_avg_dict.values(), yerr=list(pinfty_std_dict.values()), fmt='o-', capsize=3, ms=ms, mfc=mfc, color=line_color)

    PlotAxes(ax[1, 0], '$L$', r'$p_c$', r'(c) $\langle k \rangle$ = 7.5, $\zeta=10$', mode=False)
    PlotAxes(ax[1, 1], '$L$', r'$P_{\infty}(p_c)$', r'(d) $\langle k \rangle$ = 7.5, $\zeta=10$', mode=True)
    
    plt.savefig(figurepath+'/FigS3.png', dpi=500)
    plt.savefig(figurepath+'/FigS3.pdf')
    plt.savefig(figurepath+'/FigS3.eps')

if __name__  == "__main__":
    
    #set the path
    networkpath = '../kcorePercolation/figureS3/network'
    resultpath = '../kcorePercolation/figureS3/result'
    figurepath = '../kcorePercolation/figureS3/figure'
    
    #set the parameter 
    Ls = [100, 300, 500, 1000, 1500, 2000]
    zeta = 10
    slope = -1/zeta
    
    #load the link seq and degree distribution 
    LoadNetLink(networkpath, resultpath, Ls)

    #load the pc and pinfty
    ps_dict = setting_ps()
    LoadPcPinty(resultpath, ps_dict, Ls)
    
    # Plot the system size dependence as FigS3
    plot_system_size_dependency(resultpath, figurepath, Ls, slope)