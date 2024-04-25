# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 12:57:17 2023

@author: Leyang Xue
"""
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def AvgPcZetaKcore(resultpath, kcores, zetas):
    """
    Calculate the average critical probability for each k-core and zeta.
    
    Parameters
    ----------
    resultpath : str
        Path to the directory containing result files.
    kcores : list
        List of k-core values.
    zetas : list
        List of zeta values.
    
    Returns
    -------
    None
    """
    
    pc_dict = {}
    for kcore in kcores:
        kcore_pc = {}
        for zeta in zetas:
            # Load the result
            pc_result = kp.load(resultpath + f'/kcore{kcore}/kcore{kcore}_zeta{zeta}_binarysearchpc.pkl')
            kcore_pc[zeta] = np.average(pc_result)
            
        # Save the result
        pc_dict[kcore] = kcore_pc
    kp.save(resultpath + '/kcore_zeta_binary_avgpc.pkl', pc_dict)

def PlotPc2Zeta(resultpath, ax, kcore, i):
    """
    Plot the critical probability versus zeta for a specific k-core.
 
    Parameters
    ----------
    resultpath : str
        Path to the directory containing result files.
    ax : matplotlib.axes.Axes
        Axes object for plotting.
    kcore : int
        Value of k-core.
    i : int
        Index of the subplot.
 
    Returns
    -------
    None
    """
    kcore_zeta_binary_avgpc = kp.load(resultpath + '/kcore_zeta_binary_avgpc.pkl')
    
    titles = ['(a)', '(b)', '(c)', '(d)', '(e)']
    color = plt.get_cmap('Paired')
    mew = 1.2
    lw = 1.2
    ms = 8
    big_size = 12  
    
    zeta_threshold_pc = kcore_zeta_binary_avgpc[kcore]
    x = np.array(list(zeta_threshold_pc.keys()))
    y = np.array(list(zeta_threshold_pc.values()))
    max_index = np.argmax(y)
    
    ax[i].axvline(x[max_index], color = color(2*i+1), lw =lw, ls='--')
    ax[i].text(x[max_index]+0.5, min(y)+(max(y)-min(y))/10, r'$\zeta_c$='+str(x[max_index]),color= color(2*i+1), size=big_size)
    ax[i].text(50, max(y)-(max(y)-min(y))/8, str(titles[i]) + f' {kcore}-core',color='black', size=big_size)
    ax[i].plot(x, y,'o', color=color(2*i+1), ms=ms, mfc='None', mew =mew, ls='-', lw =lw)
    ax[i].set_xscale('log')
    kp.PlotAxes(ax[i], '', r'$p_c$','', mode=False)

    
if __name__ == "__main__":
    
    resultpath = '../kcorePercolation/figureS5/result'
    figurepath = '../kcorePercolation/figureS5/figure'
    
    #set the parameter
    kcores = [1,2,3,4,6]
    zetas = [3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,50,100,500,1000]
    # Calculate average critical probability for each k-core and zeta
    AvgPcZetaKcore(resultpath, kcores, zetas)
    
    # Create subplots
    fig, ax = plt.subplots(5,1, figsize=(6,8), constrained_layout=True, sharex=True) 
    
    # Plot Pc vs. Zeta for each k-core
    for i, kcore in enumerate(kcores):
        PlotPc2Zeta(resultpath, ax, kcore, i)
        
    # Set labels and ticks for the last subplot
    kp.PlotAxes(ax[4], r'$\zeta$', r'$P_c$', '')
    ax[4].set_xticks([3, 10, 100])
    ax[4].set_xticklabels([r'$3$', r'$10^1$', r'$10^2$'])
    
    # Save the figure in multiple formats
    plt.savefig(figurepath + '/FigS5.png', dpi=500)
    plt.savefig(figurepath + '/FigS5.pdf')
    plt.savefig(figurepath + '/FigS5.eps')
    plt.close()
    