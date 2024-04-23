# -*- coding: utf-8 -*-
"""
Created on Sun Feb 18 14:26:11 2024

@author: Leyang Xue 

"""

import sys
import matplotlib.pyplot as plt
import numpy as np

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def AveragePinfty(simulation_results_array):
    '''
    Calculate the average pinfty

    Parameters
    ----------
    simulation_results_array : numpy.ndarray
        Array containing simulation results.

    Returns
    -------
    pinfty_avg : numpy.ndarray
        Array containing the average pinfty.
    ''' 
    
    #set the array to store the results
    pinfty_avg = np.zeros(simulation_results_array.shape[1])
    threshold = 0.0001
    
    #identify the index of avg. pc
    avg_pc_indexs = [] 
    for each in simulation_results_array:
       pc_index =  np.argwhere(each>threshold)[0, 0]
       avg_pc_indexs.append(pc_index)
    avg_pc_index = int(round(np.average(avg_pc_indexs),0))
    
    #it will be zero below the index of avg.pc, make a average over cases that is larger than the threshold e.g. above the index of avg.pc
    pinfty = simulation_results_array[:, avg_pc_index:]
    count_nonzero = pinfty.copy()
    count_nonzero[count_nonzero>threshold] = 1
    count_nonzero[count_nonzero!=1] = 0
    nonzero_sum = np.sum(count_nonzero, axis=0)
    pinfty_avg[avg_pc_index:] =  np.sum(pinfty, axis=0)/nonzero_sum
        
    return pinfty_avg

def makeaverage(resultpath, kcore, zeta, simulation):
    '''
    Calculate the average pinfty for a given kcore and zeta.

    Parameters
    ----------
    resultpath : str
        The path to load the kcore simulation result.
    kcore : int
        The kcore value.
    zeta : int
        The zeta value.
    simulation : int
        Simulation times. 
        
    Returns
    -------
    pinfty_avg_dict : dict
        A dictionary containing pinfty values for each parameter value.
    '''
        
    simulation_results =  []
    
    for simu in np.arange(simulation):
        zeta_result = kp.load(resultpath+f'/kcore{kcore}/zeta{zeta}/kcore{kcore}_zeta{zeta}_simu{simu}_pinfty_dict.pkl')
        pinfty = [zeta_result[p] for p in sorted(zeta_result.keys())] 
        simulation_results.append(pinfty)
    
    simulation_results_array = np.array(simulation_results)
    pinfty_avg = AveragePinfty(simulation_results_array)
    ps = sorted(zeta_result.keys())
    pinfty_avg_dict = {p:pinfty for p, pinfty in zip(ps, pinfty_avg)}
    
    return pinfty_avg_dict
    
def plot_pinfty2p(ax, resultpath, kcore, zetas, simu, PhaseTypes, xlim_min, xlim_max, ylim_min, ylim_max, title, xlabel, ylabel, simualtion):
    '''
    Plot pinfty for different zeta values.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes object to plot on.
    resultpath : str
        The path to the result files.
    kcore : int
        The kcore value.
    zetas : list
        List of zeta values.
    simu : int
        The simulation number.
    PhaseTypes : list
        List of phase types.
    xlim_min : float
        Minimum value for x-axis limit.
    xlim_max : float
        Maximum value for x-axis limit.
    ylim_min : float
        Minimum value for y-axis limit.
    ylim_max : float
        Maximum value for y-axis limit.
    title : str
        Title of the plot.
    xlabel : str
        Label for x-axis.
    ylabel : str
        Label for y-axis.
    simulation : int
        Number of the simulation used to make the average.

    Returns
    -------
    None
    '''
    
    #setting of parameter 
    color = plt.get_cmap('Paired')
    color1 = plt.get_cmap('tab20')
    ms = 6
    mew  = 0.85
    mfcc = 'None'
    
    #plot the lcc for different zeta with a given kcore
    for j, zeta in enumerate(zetas):
        
        if PhaseTypes[j] == 'CT':
            zeta_result = makeaverage(resultpath, kcore, zeta, simualtion)
        elif PhaseTypes[j] == 'DT':
            #load the results for different kcore and zeta
            zeta_result = kp.load(resultpath+f'/kcore{kcore}/zeta{zeta}/kcore{kcore}_zeta{zeta}_simu{simu}_pinfty_dict.pkl')
        
        if j == 0:
            ax.plot(zeta_result.keys(), zeta_result.values(), 'o-', ms= ms, color=color1(14), mfc=mfcc, mew = mew, label= r'$\zeta=$'+str(zeta))
        else:    
            ax.plot(zeta_result.keys(), zeta_result.values(), 'o-', ms= ms, color=color(10-2*j+1), mfc=mfcc, mew = mew, label= r'$\zeta=$'+str(zeta))
    
    # Set labels and title
    kp.PlotAxes(ax, xlabel, ylabel, title, mode=False) 
    ax.text(xlim_min+(xlim_max-xlim_min)*1/10, ylim_min+(ylim_max-ylim_min)*8/9,f'{kcore}-core', size=12)
    ax.set_xlim(xlim_min, xlim_max)
    ax.set_ylim(ylim_min, ylim_max)
    
if __name__ == "__main__":

    #set the path 
    resultpath = '../kcorePercolation/figureS4/result'
    figurepath = '../kcorePercolation/figureS4/figure'
    networkpath = '../kcorePercolation/figureS4/network'
    
    #set the basic information of spatial networks
    kcores = [2,3,4,6]
    zetas = [500, 14, 10, 8, 6, 3]
    simu = 1
    PhaseTypes_kcore = {2:['CT','CT','CT','CT','CT','CT'], 3:['DT','DT','CT','CT','CT','CT'],4:['DT','DT','DT','DT','CT','CT'], 6:['DT','DT','DT','DT','DT','CT']}    
    titles = ['(a)', '(b)', '(c)', '(d)']
    xlabels = ['','',r'$p$',r'$p$']
    ylabels = [r'$P_{\infty}(p)$', '',  r'$P_{\infty}(p)$', '']
    
    xlim_mins = [0.095,0.305, 0.485, 0.78]
    xlim_maxs = [0.205,0.355, 0.565, 0.91]
    ylim_mins = [0.0, 0.0, 0.0, 0.0]
    ylim_maxs = [0.10,0.20, 0.45, 0.8]
    simulation = 10
    
    color = plt.get_cmap('Paired')
    color1 = plt.get_cmap('tab20c')
    
    # Create subplots and plot data
    fig, ax = plt.subplots(2,2, figsize=(6.5, 5.8), constrained_layout=True)
    for s, kcore in enumerate(kcores):
        i = int(s/2)
        j = int(s%2)
        PhaseTypes = PhaseTypes_kcore[kcore]
        plot_pinfty2p(ax[i,j], resultpath, kcore, zetas, simu, PhaseTypes, xlim_mins[s], xlim_maxs[s], ylim_mins[s], ylim_maxs[s], titles[s], xlabels[s], ylabels[s], simulation)
        if s == 3:
            ax[1,1].legend(loc='lower left', bbox_to_anchor=(0.05,0.05,0.4,0.4),ncol=1,framealpha=0, fontsize=10,borderpad=0.0,handlelength=1.0)
    
    # Save the plot
    plt.savefig(figurepath+'/FigS4.png', dpi = 500)
    plt.savefig(figurepath+'/FigS4.pdf') 
    plt.savefig(figurepath+'/FigS4.eps') 
    
