# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 13:47:12 2024

@author: Leyang Xue
"""

import sys
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def logbins(data, bins_num):
    """
    Create logarithmic bins for histogram.

    Parameters
    ----------
    data : array-like
       Input data.
    bins_num : int
       Number of bins.

    Returns
    -------
    hist : array
       Histogram values.
    bin_edges : array
       Bin edges.
    """
    bins = np.logspace(np.log10(min(data)), np.log10(max(data)), num=bins_num)
    hist, bin_edges = np.histogram(data, bins)
        
    return hist, bin_edges

def CalculateDist(links, k, zeta, resultpath):
    """
    Calculate the distribution of link lengths.

    Parameters
    ----------
    links : array-like
        Lengths of links.
    k : int
        Average degree.
    zeta : int
        Zeta value.
    resultpath : str
        Path to save the results.

    Returns
    -------
    bin_centers : array
        Bin centers.
    pdf : array
        Probability density function.
    cdf_reverse : array
        Reverse cumulative distribution function.
    """    
    hist, edges = np.histogram(links, bins=50)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])
    pdf = hist/np.sum(hist)
    pdf_reverse = pdf[::-1]
    cdf_reverse = np.cumsum(pdf_reverse)[::-1]
    
    kp.save(resultpath + f'/bin_centers_k{k}_zeta{zeta}.pkl', bin_centers)
    kp.save(resultpath + f'/pdf_k{k}_zeta{zeta}.pkl', pdf)
    kp.save(resultpath + f'/cdf_reverse_k{k}_zeta{zeta}.pkl', cdf_reverse)

    return bin_centers, pdf, cdf_reverse

def CalculateDegreeDist(network, k, zeta, resultpath):
    """
    Calculate the distribution of node degrees.

    Parameters
    ----------
    network : NetworkX graph
        Input network.
    k : int
        Average degree.
    zeta : int
        Zeta value.
    resultpath : str
        Path to save the results.

    Returns
    -------
    probability_distribution : dict
        Degree distribution.
    """
    degree_seq = list(dict(network.degree()).values())
    degree_distribution = Counter(degree_seq)
    total_nodes = sum(degree_distribution.values())
    
    # Calculate probability distribution
    probability_distribution = {k: degree_distribution[k] / total_nodes for k in sorted(degree_distribution.keys())}
    kp.save(resultpath + f'/degree_distribution_k{k}_zeta{zeta}.pkl', probability_distribution)

    return probability_distribution
    
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
     
def TheoreticalLine(x_data, slope):
    """
    Generate theoretical line data.
  
    Parameters
    ----------
    x_data : array-like
        x-axis data.
    slope : float
        Slope of the line.
  
    Returns
    -------
    y_data : array-like
        y-axis data.
    """
    y_data  = np.exp(slope * x_data)    
    return y_data

def PoissonDistribution(mean_value = 7.5):
    """
    Generate Poisson distribution data.

    Parameters
    ----------
    mean_value : float, optional
        Mean value for the Poisson distribution. The default is 7.5.

    Returns
    -------
    x : array-like
        x-axis data.
    poisson_distribution : array-like
        Poisson distribution data.
    """
    poisson_sequence = np.random.poisson(mean_value, size=1000000)
    poisson_distribution = np.bincount(poisson_sequence) / len(poisson_sequence)
    x = range(len(poisson_distribution))
    
    return x, poisson_distribution
 
def plot_figureS2(networkpath, resultpath, figurepath, zetas, avgks):

    # Set the plot parameters 
    colors = plt.get_cmap("Paired")
    line_color = plt.get_cmap("tab20c")(17)
    #slope_line_color = plt.get_cmap("tab20c")(17)

    lw = 1.5
    mfc = 'None'
    n_legend = 12
    ms = 6.5
    mks = ['o', 's', '^', '*'] 
    
    fig, ax = plt.subplots(2, 2, figsize=(7.5, 6), constrained_layout = True)

    ax[0,0].sharex(ax[1,0])
    ax[0,1].sharex(ax[1,1])
    
    for i, k in enumerate(avgks):
        
        # Plot the Poisson distribution
        [x, poisson_dist] = PoissonDistribution(mean_value = k)
        ax[i,1].plot(x, poisson_dist, '-', color=line_color, lw=lw, label=r'Poisson($\lambda$'+f'={k})') 
            
        for j, zeta in enumerate(zetas):  
       
             # # Load the network result 
             # links = kp.load(networkpath + f'/network_L_1000_avg_k_{k}_zeta_{zeta}_links.pkl')
             # network = kp.load(networkpath + f'/network_L_1000_avg_k_{k}_zeta_{zeta}_spatialNet.pkl')
             
             # # Calculate the distribution of links 
             # [bin_centers, pdf, cdf_reverse] = CalculateDist(links, k, zeta, resultpath)
             # # Calculate the degree distribution
             # degree_distribution = CalculateDegreeDist(network, k, zeta, resultpath)
             
             # Load the figure result 
             bin_centers = kp.load(resultpath + f'/bin_centers_k{k}_zeta{zeta}.pkl')
             pdf = kp.load(resultpath + f'/pdf_k{k}_zeta{zeta}.pkl')
             #cdf_reverse = kp.load(resultpath + f'/cdf_reverse_k{k}_zeta{zeta}.pkl')
             degree_distribution = kp.load(resultpath + f'/degree_distribution_k{k}_zeta{zeta}.pkl')
             
             # Plot the distribution of link length 
             ax[i,0].plot(bin_centers, pdf, mks[j], ms=ms, color=colors(2*j+1), mfc = mfc) 
             y_data = TheoreticalLine(bin_centers, slope=-1/zeta)
             ax[i,0].plot(bin_centers, y_data/sum(y_data), '-', color=colors(2*j+1), lw=lw, label=r'slope=-{:.2f},'.format(1/zeta)+r' $\zeta$'+f'={zeta}')
             
             # Plot the distribution of degree
             if i == 0:
                 ax[i,1].plot(degree_distribution.keys(), degree_distribution.values(), mks[j], color=colors(2*j+1), mfc = mfc, label=r'$\zeta$='+ f'{zeta}') 
             else:
                 ax[i,1].plot(degree_distribution.keys(), degree_distribution.values(), mks[j], color=colors(2*j+1), mfc = mfc) 

        ax[i,0].set_yscale('log')
        
    PlotAxes(ax[0,0],'',r'$p(l)$', r'(a) $\langle k \rangle=5$', mode=True)
    PlotAxes(ax[0,1],'',r'$p(k)$', r'(b) $\langle k \rangle=5$', mode=False)
    PlotAxes(ax[1,0],'link length, '+r'$l$',r'$p(l)$', r'(c) $\langle k \rangle=10$', mode=False)
    PlotAxes(ax[1,1],'degree, '+r'$k$',r'$p(k)$', r'(d) $\langle k \rangle=10$', mode=False)

    ax[0,1].legend(loc='upper center', bbox_to_anchor=(0.7, 0.93), framealpha=0, fontsize=n_legend)
    ax[1,1].legend(loc='upper center', bbox_to_anchor=(0.7, 0.93), framealpha=0, fontsize=n_legend)

    plt.savefig(figurepath+'/FigS2.png', dpi=500)
    plt.savefig(figurepath+'/FigS2.pdf')
    plt.savefig(figurepath+'/FigS2.eps')
       
if __name__  == "__main__":
    
    # Set the path
    networkpath = '../kcorePercolation/figureS2/network'
    resultpath = '../kcorePercolation/figureS2/result'
    figurepath = '../kcorePercolation/figureS2/figure'
    
    # Set the zeta
    zetas = [4, 5, 10, 100]
    avgks= [5, 10]
    
    # Plot the figureS2
    plot_figureS2(networkpath, resultpath, figurepath, zetas, avgks)
    