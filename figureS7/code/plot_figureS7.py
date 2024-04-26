# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 03:56:24 2024

@author: Leyang Xue

"""
import sys
import numpy as np
import matplotlib.pyplot as plt

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def SetAxShare(ax):
    '''
    Set up shared x and y axes for subplots and remove tick labels as specified.

    Parameters
    ----------
    ax : array_like
        Array of axes objects representing subplots.

    Returns
    -------
    None.

    '''
    # Share y axis
    ax[0,1].sharey(ax[0,0])
    ax[1,1].sharey(ax[1,0])
    ax[2,1].sharey(ax[2,0])
    
    # Share x axis for one column
    ax[2,0].sharex(ax[1,0])
    ax[1,0].sharex(ax[0,0])
    
    # Share x axis for two columns
    ax[2,1].sharex(ax[1,1])
    ax[1,1].sharex(ax[0,1])
    
    # Set xticks as null
    plt.setp(ax[0,0].get_xticklabels(), visible=False)
    plt.setp(ax[0,1].get_xticklabels(), visible=False)
    plt.setp(ax[1,0].get_xticklabels(), visible=False)
    plt.setp(ax[1,1].get_xticklabels(), visible=False)
    
    # Set yticks as null
    plt.setp(ax[0,1].get_yticklabels(), visible=False)
    plt.setp(ax[1,1].get_yticklabels(), visible=False)
    plt.setp(ax[2,1].get_yticklabels(), visible=False)

def calculateSt(result):
    '''
    Calculate st from the result.

    Parameters
    ----------
    result : dict
        Dictionary containing the result.

    Returns
    -------
    x : array_like
        Array containing x values.
    st : array_like
        Array containing st values.

    '''
    x = np.array(list(result.keys()))[1:]
    y = np.array(list(result.values()))
    st = y[:-1] - y[1:]
    
    return x, st

def calculateetat(t_st, st):
    '''
    Calculate etat from t_st and st.

    Parameters
    ----------
    t_st : array_like
        Array containing t_st values.
    st : array_like
        Array containing st values.

    Returns
    -------
    t_etat : array_like
        Array containing t_etat values.
    etat : array_like
        Array containing etat values.

    '''
    etat = st[1:]/st[:-1]
    t_etat = t_st[1:]
    
    return t_etat, etat

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
    fontsize = 14
    font_label = {'family': "Arial", 'size':fontsize}
    
    n_legend = 10
    ax.set_xlabel(xlabel,  fontdict = font_label)
    ax.set_ylabel(ylabel, fontdict = font_label)
    ax.set_title(title, loc='left',fontdict = {'family': "Arial", 'size':fontsize})
    ax.tick_params(direction='out', which='both',length =4, width=1, pad=1,labelsize=n_legend)
    
    #ax.minorticks_on()
    if mode == True:
        ax.legend(loc='best', framealpha=0, fontsize=n_legend)
        
def PlotDynamics(resultpath, figurepath, kcore, zetas, zeta_pc, deltap, simu_dict, zeta_labels, texts, figname):
    
    #plot the plateau
    titles = ['(a)', '(b)', '(c)', '(d)', '(e)', '(f)','(g)', '(h)']
    color = plt.get_cmap("Paired")
    ms = 6
    sms = 6
    mew = 0.4
    lw  = 1
    mfc = 'None'
    size = 14
    
    fig, ax = plt.subplots(3,2, figsize=(5.5,6.25), constrained_layout=True)
    SetAxShare(ax)
    
    #set the filenames with different ps
    markers = ['o', '^', 'X','*', 's']

    for i, zeta in enumerate(zetas):
        pc = zeta_pc[zeta] #get the pc
        deltap_zeta = deltap[zeta]
        labels = zeta_labels[zeta]
        simu = simu_dict[zeta]

        #for a give zeta 10 to deal with 
        for j, delta_p in enumerate(deltap_zeta):
            #update the p 
            p = round(pc - delta_p,6)
            
            #load the dynamic process at the criticality 
            if p == pc:
                zeta_pc_pinfty_dynamic = kp.load(resultpath+f'/zeta{zeta}/kcore{kcore}_zeta{zeta}_p{pc}_simu{simu}_steps_lcc_num_dict.pkl')
                zeta_pc_fgcc_dynamic = kp.load(resultpath+f'/zeta{zeta}/kcore{kcore}_zeta{zeta}_p{pc}_simu{simu}_steps_fgcc_nodes_dict.pkl')
            else:    
                #load the dynamic process below the criticality 
                zeta_pc_pinfty_dynamic = kp.load(resultpath+f'/zeta{zeta}/kcore{kcore}_zeta{zeta}_dp{delta_p}_p{p}_simu{simu}_steps_lcc_num_dict.pkl')
                zeta_pc_fgcc_dynamic = kp.load(resultpath+f'/zeta{zeta}/kcore{kcore}_zeta{zeta}_dp{delta_p}_p{p}_simu{simu}_steps_fgcc_nodes_dict.pkl')
                
            #calculate the st and etat of zeta100
            [t_st, st] = calculateSt(zeta_pc_pinfty_dynamic)
            [t_etat, etat]= calculateetat(t_st, st)
            
            #plot the sfig for pinfty, fgcc, etat
            ax[0,i].plot(list(zeta_pc_pinfty_dynamic.keys())[1:], list(zeta_pc_pinfty_dynamic.values())[1:], marker = markers[j], ls = '-', color=color(2*j+1), ms = ms, mfc=mfc, mew = mew)
            ax[1,i].plot(list(zeta_pc_fgcc_dynamic.keys())[1:], list(zeta_pc_fgcc_dynamic.values())[1:], marker = markers[j], ls = '-', color=color(2*j+1), ms = ms,  mfc=mfc, mew = mew,  label=labels[j])
            ax[2,i].plot(t_etat, etat, marker = markers[j], ls = '-', color=color(2*j+1), ms = sms, mfc=mfc, mew = mew, lw=lw)

    ax[0,0].text(70, 0.62, texts[0], color='black', size=size)
    ax[0,1].text(50, 0.62, texts[1], color='black', size=size)
    
    #set the ticks and titles
    PlotAxes(ax[0,0], '', r'$P_{\infty}(t)$', titles[0])
    PlotAxes(ax[0,1], '', '', titles[1])
    PlotAxes(ax[1,0], '', r'$P^{f}_{\infty}(t)$', titles[2])
    PlotAxes(ax[1,1], '', '', titles[3])
    PlotAxes(ax[2,0], 'Step t', r'$\eta_t$', titles[4])
    PlotAxes(ax[2,1], 'Step t', '', titles[5])
    
    #set the legend
    ax[1,1].legend(loc='lower center', bbox_to_anchor=(0.55, -0.05), framealpha=0, ncol=2, fontsize=8,  labelspacing=0.15, columnspacing=0.2)# bbox_to_anchor=(0.05,-0.05)
    ax[1,0].legend(loc='lower center', bbox_to_anchor=(0.55, -0.05),framealpha=0, ncol=2, fontsize=8,  labelspacing=0.15, columnspacing=0.2)# bbox_to_anchor=(0.05,-0.05)
    ax[0,1].set_xlim(-13,305)
    ax[2,0].set_ylim(-0.1,2.1)
    
    plt.savefig(figurepath+f'/{figname}.png', dpi = 500)
    plt.savefig(figurepath+f'/{figname}.pdf')
    plt.savefig(figurepath+f'/{figname}.eps')
    

if __name__ == "__main__":
    
    # Set the result, network, and figure paths
    figurepath  = '../kcorePercolation/figureS7/figure'
    resultpath = '../kcorePercolation/figureS7/result'
    networkpath = '../kcorePercolation/figureS7/network'
    
    # Set the value of the k-core
    kcore = 5  

    # Set the parameter for zeta = 10 and 100
    zetas = [10, 100]
    deltap = {10:np.array([0,0.002,0.005,0.01,0.1]), 100:np.array([0,0.0001,0.0005, 0.001, 0.005])}
    zeta_pc = {10:0.700501,100:0.679701}
    simu_dict = {10:1, 100:2}
    texts = [r'$\zeta$=10', r'$\zeta=100$']
    zeta_labels =  {10:['$p_c$', '$p_c-0.002$', '$p_c-0.005$', '$p_c-0.01$', '$p_c-0.1$'],100:['$p_c$', '$p_c-0.0001$', '$p_c-0.0005$', '$p_c-0.001$', '$p_c-0.005$']}
    figname  = 'FigS7'
    
    # Plot the figures
    PlotDynamics(resultpath, figurepath, kcore, zetas, zeta_pc, deltap, simu_dict, zeta_labels, texts, figname)
        
    