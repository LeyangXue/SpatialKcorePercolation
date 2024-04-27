# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 00:46:20 2023

@author: Leyang Xue
"""
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def loadAvgkKcoreZetaRhc(resultpath,kcore,avgks,zetas):
    """
    Load average k, kcore, zeta, and Rhc data from files.
    
    Parameters
    ----------
    resultpath : str
        Path to the result directory.
    kcore : int
        The kcore value.
    avgks : list
        List of average k values.
    zetas : list
        List of zeta values.
    
    Returns
    -------
    dict
        A dictionary containing the RHC data for each combination of average k and zeta.
    """
    # Initialize an empty dictionary to store RHC data
    rhc_avgk_dict = {}
    
    # Iterate over each average k value
    for avgk in avgks:
        # Initialize an empty dictionary to store RHC data for the current average k value
        rhc_zeta_dict = {}
        
        # Iterate over each zeta value
        for zeta in zetas:
            # Construct the file name
            file = 'avgk'+str(avgk)+'_kcore'+str(kcore)
            filename = file+'_zeta'+str(zeta)+'.pkl'
            
            # Load RHC data from the file
            rhc = kp.load(resultpath+'/'+file+'/'+filename)
            
            # Check if RHC data is valid
            if rhc[1] != -1:    
                rhc_zeta_dict[rhc[0]] = rhc[1]
            else:
                print(f"exception with avgk: {avgk} and zeta: {zeta}")
                rhc_zeta_dict[rhc[0]] = rhc_zeta_dict[zeta-1]
        
        # Store RHC data for the current average k value in the main dictionary      
        rhc_avgk_dict[avgk] = rhc_zeta_dict
    
    return rhc_avgk_dict

def process_datasets(Lfiles, resultpath, avgks, zetas, kcores):
    """
    Process datasets for different L values and kcores.

    Parameters
    ----------
    Lfiles : list
        List of L values.
    resultpath : str
        Path to the result directory.
    avgks : list
        List of average k values.
    zetas : list
        List of zeta values.
    kcores: list
        List of kcore values.
    Returns
    -------
    None
    """
    
    # Iterate over each L value
    for L in Lfiles:
        
        # Define the path to load data for the current L value
        loadpath = resultpath + '/L' + str(L)
        
        # Iterate over each kcore value
        for kcore in kcores:
            print(f'Loading files with L: {L} and kcore: {kcore}')
            # Load data and process
            kcores_rhc = loadAvgkKcoreZetaRhc(loadpath, kcore, avgks, zetas)
            # Convert processed data to a DataFrame
            kcores_rhc_pd = pd.DataFrame.from_dict(kcores_rhc, orient='index')
            # Define the path to save the processed data
            output_file = resultpath + '/total/L' + str(L) + '_kcore' + str(kcore) + '.csv'
            # Save processed data to a CSV file
            kcores_rhc_pd.to_csv(output_file)
            
            # Print a message confirming the saved file
            print(f'Saved processed data to: {output_file}')
            
def Transdict2Array(rootpath, L, kcore, avgks, zetas):
    """
    Transform pandas DataFrame into a matrix and create mappings for average k and zeta values.
    
    Parameters
    ----------
    rootpath : str
        The root path where the CSV file is located.
    L : int
        The L value.
    kcore : int
        The kcore value.
    avgks : list
        List of average k values.
    zetas : list
        List of zeta values.
    
    Returns
    -------
    tuple
        A tuple containing the matrix representation of RHC data, and dictionaries mapping average k and zeta values to indices.
    """
    
    # Load the results from the CSV file
    data = pd.read_csv(rootpath+'/total/L'+str(L)+'_kcore'+str(kcore)+'.csv', index_col=0)
    
    # Initialize dictionaries for mappings
    avgk2index = {}
    zeta2index = {}
    
    # Initialize the matrix to store RHC data
    kcore_array = np.zeros((len(avgks), len(zetas)))
    
    # Iterate over average k values
    for i, avgk in enumerate(avgks):
        
        # Build the mapping between avgk and index
        avgk2index[avgk] = i
        zetas_rhc = data.loc[avgk]
        
        # Iterate over zeta values
        for j, zeta in enumerate(zetas):
            # Build the mapping between zeta and index
            if i == 0:
                zeta2index[zeta] = j
            
            # Fill the matrix with RHC data
            kcore_array[i,j] = zetas_rhc.loc[str(zeta)]

    return kcore_array, avgk2index, zeta2index

def PlotHeatMap(rootpath, ax, L, kcore, xlabel, ylabel, title, cbars=False):
    """
     Plot a heatmap of RHC data.
    
     Parameters
     ----------
     rootpath : str
         The root path where the rhc array and mappings files are located.
     ax : plt.Axes
         The matplotlib Axes object where the heatmap will be plotted.
     L : int
         The L value.
     kcore : int
         The kcore value.
     xlabel : str
         The label for the x-axis.
     ylabel : str
         The label for the y-axis.
     title : str
         The title of the plot.
     cbars : bool, optional
         Whether to show the colorbar, by default False.
    
     Returns
     -------
     None
     
    """
    # Load rhc array and mappings
    kcores_rhc_array = kp.load(rootpath+'/total/L'+str(L)+'_kcore'+str(kcore)+'_rhc_array.pkl')
    avgk2index = kp.load(rootpath+'/total/L'+str(L)+'_kcore'+str(kcore)+'_avgk2index.pkl')
    
    # Choose colormap
    cmap = "turbo" #cividis,viridis,plasma,coolwarm,YlGnBu,magma,mako, #turbo,crest_r viridis,magma,jet
    
    # Plot heatmap
    sns.heatmap(kcores_rhc_array, ax=ax, vmin=0, vmax=250, cbar=cbars, cbar_kws={'label': '$r^c_h$'},cmap=cmap)
    
    # Set y-axis ticks and labels
    yticks = np.arange(0,len(avgk2index.keys()),10)
    yticklabels = [int(each) for each in np.array(list(avgk2index.keys()))[yticks]]
    ax.set_yticks(yticks+0.5)
    ax.set_yticklabels(yticklabels,rotation='horizontal')
        
    # Set x-axis ticks and labels
    xticks = np.arange(0,len(zetas),15)
    xticklabels = zetas[xticks]
    ax.set_xticks(xticks+0.5)
    ax.set_xticklabels(xticklabels)
    
    # Set axis labels and title
    kp.PlotAxes(ax, xlabel, ylabel, title)

def plot_figure(resultpath, figurepath, Lfiles, kcores):
    """
    Plot the results and save the figure.

    Parameters
    ----------
    resultpath : str
        Path to the result directory.
    figurepath : str
        Path to save the figure.
    Lfiles: list
        Size of different system
    kcores: list
        Different kcore values
        
    Returns
    -------
    None
    """
    
    # Titles for subplots
    titles = [['(a1)', '(a2)', '(a3)'], ['(b1)', '(b2)', '(b3)'], ['(c1)', '(c2)', '(c3)'], ['(d1)', '(d2)', '(d3)']]
   
    # Create subplots
    fig, ax = plt.subplots(3, 3, figsize=(9, 9), sharex=True, sharey=True, constrained_layout=True)
    
    # Iterate over L values and kcore values
    for i, L in enumerate(Lfiles):
        for j, kcore in enumerate(kcores):
            if L != 1000 or (L == 1000 and (kcore != 5 or kcore != 6)):
                # Plot heatmap for each combination of L and kcore
                if i != 2 and j == 0:
                    PlotHeatMap(resultpath, ax[i, j], L, kcore, '',  ' \n', titles[i][j] + f' L={L},  {kcore}-core', cbars=False)
                elif i != 2 and j == 2:
                    PlotHeatMap(resultpath, ax[i, j], L, kcore, '',  '', titles[i][j] + f' L={L},  {kcore}-core', cbars=True)
                elif i == 2 and j == 0:
                    PlotHeatMap(resultpath, ax[i, j], L, kcore, ' \n',  ' \n', titles[i][j] + f' L={L},  {kcore}-core', cbars=False)
                elif i == 2 and j == 2:
                    PlotHeatMap(resultpath, ax[i, j], L, kcore, ' \n',  '', titles[i][j] + f' L={L},  {kcore}-core', cbars=True)
                elif i == 2 and (j == 1 or j == 1):
                    PlotHeatMap(resultpath, ax[i, j], L, kcore, ' \n',  '', titles[i][j] + f' L={L},  {kcore}-core', cbars=False)
                else:
                    PlotHeatMap(resultpath, ax[i, j], L, kcore, '',  '', titles[i][j] + f' L={L},  {kcore}-core', cbars=False)
    
    # Set labels and save the figure
    fontsize = 14
    font_label = {'family': "Arial", 'size': fontsize}
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel(r"$\zeta$", labelpad=10, fontdict=font_label)
    plt.ylabel(r"average degree, $\langle k \rangle$", labelpad=10, fontdict=font_label)
    
    plt.savefig(figurepath+'/FigS9.png', dpi=500)
    plt.savefig(figurepath+'/FigS9.pdf')
    plt.savefig(figurepath+'/FigS9.eps')
    
            
if __name__ == "__main__":
    
    # Set the result, network, and figure paths
    resultpath = '../kcorePercolation/figureS9/result' 
    figurepath = '../kcorePercolation/figureS9/figure'
    
    #one part: load the results
    Lfiles = [300, 500, 700, 1000]#system size
    zetas = np.arange(3,101,1)
    kcores = [3, 4, 5, 6] # Define the list of kcores
    avgks = [round(each,1) for each in np.arange(5.0, 10.01, 0.1)]
    
    # Process the data and transform into the pandas
    process_datasets(Lfiles, resultpath, avgks, zetas, kcores)

    # Transform pd into the array
    avgks = [round(each,1) for each in np.arange(10.0, 4.95, -0.1)]
    for L in Lfiles:
        for kcore in kcores:
            [kcores_rhc_array, avgk2index, zeta2index] = Transdict2Array(resultpath, L, kcore, avgks, zetas)
            kp.save(resultpath+'/total/L'+str(L)+'_kcore'+str(kcore)+'_rhc_array.pkl', kcores_rhc_array)
            kp.save(resultpath+'/total/L'+str(L)+'_kcore'+str(kcore)+'_avgk2index.pkl', avgk2index)
            kp.save(resultpath+'/total/L'+str(L)+'_kcore'+str(kcore)+'_zeta2index.pkl', zeta2index)

    # Plot the figure 
    kcores = [4, 5, 6]
    Lfiles = [500, 700, 1000]#system size
    plot_figure(resultpath, figurepath, Lfiles, kcores)