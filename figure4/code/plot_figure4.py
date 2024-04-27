# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 00:46:20 2023

@author: Leyang Xue

"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.image as mpimg

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def loadAvgkKcoreZetaRhc(resultpath, kcore, avgks, zetas):
    '''
    Load average k-core percolation results for different zeta values.

    Parameters
    ----------
    resultpath : str
        Path to the folder containing result files.
    kcore : int
        Value of the k-core.
    avgks : list
        List of average k values.
    zetas : list
        List of zeta values.

    Returns
    -------
    rhc_zeta_avgk_core3_dict : dict
        A dictionary containing average k-core percolation results for different zeta values.
        Keys are average k values, and values are dictionaries where keys are zeta values
        and values are the corresponding to the critical radius of hole under localized attack.
    '''
    # Initialize an empty dictionary to store the results
    rhc_zeta_avgk_core_dict = {}
    
    for avgk in avgks:
        # Initialize an empty dictionary for each average k value
        rhc_zeta_dict = {}
        for zeta in zetas:
            # Construct the filename based on average k, k-core value, and zeta
            file = 'avgk'+str(avgk)+'_kcore'+str(kcore)
            filename = file+'_zeta'+str(zeta)+'.pkl'
            # Load the percolation results from the file
            rhc = kp.load(resultpath+'/L1000/'+file+'/'+filename)
            # Store the percolation result in the dictionary with zeta as key
            rhc_zeta_dict[rhc[0]] = rhc[1]
            
        # Store the dictionary for the current average k value in the main dictionary
        rhc_zeta_avgk_core_dict[avgk] = rhc_zeta_dict
    
    return rhc_zeta_avgk_core_dict

def Transdict2Array(kcores_rhc, avgks, zetas):
    '''
    Transform the dictionary of k-core percolation results into a numpy array.
    
    Parameters
    ----------
    kcores_rhc : dict
        A dictionary containing k-core percolation results for different average k values.
        Keys are average k values, and values are dictionaries where keys are zeta values
        and values are the corresponding to the critical radius of hole under localized attack.
    avgks : list
        List of average k values.
    zetas : list
        List of zeta values.
    
    Returns
    -------
    kcore_array : numpy.ndarray
        2D array containing k-core percolation results.
    avgk2index : dict
        A mapping from average k values to array indices.
    zeta2index : dict
        A mapping from zeta values to array indices.
    '''

    # Initialize dictionaries to store mappings from values to indices
    avgk2index = {}
    zeta2index = {}
    
    # Initialize a 2D numpy array to store k-core percolation results
    kcore_array = np.zeros((len(avgks), len(zetas)))
    
    # Iterate over average k values
    for i, avgk in enumerate(avgks):
        # Build the map between avgk and index 
        reversed_x = int(len(avgks)-1-i)
        avgk2index[avgk] = reversed_x
        
        # Get k-core percolation results for the current average k value
        zetas_rhc = kcores_rhc[avgk]
        
        # Iterate over zeta values
        for j, zeta in enumerate(zetas):
            # Build the map between zeta and index 
            if i == 0:
                zeta2index[zeta] = j
            # Store the k-core percolation result in the array
            kcore_array[reversed_x,j] = zetas_rhc[zeta]

    return kcore_array, avgk2index, zeta2index

def Transdict2ArrayKcores(kcores_rhc, avgks, zetas):
    
    #transform the dict into the matrix
    kcores_rhc_array_dict = {}
    avgk2index = {}
    zeta2index = {}
    for its, kcore in enumerate(kcores_rhc.keys()):
        kcore_array = np.zeros((len(avgks), len(zetas)))
        avgk_zetas_rhc = kcores_rhc[kcore]
        for i, avgk in enumerate(avgks):
            #build the map between avgk and index 
            reversed_x = int(len(avgks)-1-i)
            if its == 0:
               avgk2index[avgk] = reversed_x
            zetas_rhc = avgk_zetas_rhc[avgk]
            for j, zeta in enumerate(zetas):
                #build the map between avgk and index 
                if its == 0 and i == 0:
                   zeta2index[zeta] = j
                kcore_array[reversed_x,j] = zetas_rhc[zeta]
        #save the result in a array 
        kcores_rhc_array_dict[kcore] = kcore_array
    
    return kcores_rhc_array_dict, avgk2index, zeta2index

def PlotSchematic1(ax, figurepath):
    '''
    Plot a schematic image on the given axes.
    
    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to plot the schematic image on.
    figurepath : str
        Path to the folder containing the schematic image file.
    
    Returns
    -------
    None
    '''

    # Load the schematic image
    figure = mpimg.imread(figurepath + '/schematic.png')
    
    # Plot the image on the axes
    ax.imshow(figure,interpolation='bicubic')
    
    # Set aspect ratio, turn off axis, and remove ticks
    ax.set_aspect('auto')
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])   
    
    # Add annotation to the axes
    kp.PlotAxes(ax,'','','(b)')

def Identify_curve(matrix):
    '''
     Identify the top and bottom curves in a matrix.
    
     Parameters
     ----------
     matrix : numpy.ndarray
         2D array containing data.
    
     Returns
     -------
     rows_bottom : numpy.ndarray
         Array containing row indices of the bottom curve.
     rows_top : numpy.ndarray
         Array containing row indices of the top curve.
    '''
    
    # Define initial values for max and min row indices
    max_row_value = 1000
    min_row_value = -1
    
    # Initialize lists to store row indices for bottom and top curves
    rows_bottom = []
    rows_top = []
    
    # Iterate over columns of the matrix
    for i in np.arange(matrix.shape[1]):
       columns = matrix[:,i]
       
       # Identify the bottom curve
       zero_index = np.where(columns==0)[0]
       if len(zero_index) != 0:
           rows_bottom.append(zero_index[0])
       else:
           rows_bottom.append(max_row_value)
       
       # Identify the top curve
       zero_index = np.where(columns==500)[0]
       if len(zero_index) != 0:
          if zero_index[-1] != (matrix.shape[0]-1):
              rows_top.append(zero_index[-1]) # Add the last zero index
          else:
              rows_top.append(max_row_value) # Add a maximum row value if the last element is 500
       else:
          rows_top.append(min_row_value)  # Add a minimum row value if no 500 is found
          
    # Convert lists to numpy arrays
    return np.array(rows_bottom), np.array(rows_top)

def Tran2avgK(kcore_top, kcore_bottom, index2avgk, zetas):
    '''
    Translate k-core indices to corresponding average degrees for top and bottom curves.
    
    Parameters
    ----------
    kcore_top : list
        List of k-core indices for the top curve.
    kcore_bottom : list
        List of k-core indices for the bottom curve.
    index2avgk : dict
        Dictionary mapping k-core indices to average degrees.
    zetas : list
        List of zeta values.
    
    Returns
    -------
    top_value : dict
        Dictionary containing zeta values as keys and corresponding average degrees for the top curve.
    bottom_value : dict
        Dictionary containing zeta values as keys and corresponding average degrees for the bottom curve.
    '''
    
    # Initialize dictionaries to store zeta-average degree mappings for top and bottom curves
    top_value = {}
    bottom_value = {}
    
    # Iterate over zeta values and corresponding k-core indices for top curve
    for zeta, each in zip(zetas, kcore_top):
        # Check if the k-core index exists in the mapping
        if each in index2avgk.keys():
           top_value[zeta]= index2avgk[each] # Map k-core index to average degree for the top curve
           
    # Iterate over zeta values and corresponding k-core indices for bottom curve
    for zeta, each in zip(zetas, kcore_bottom):
        # Check if the k-core index exists in the mapping
         if each in index2avgk.keys():
            bottom_value[zeta] = index2avgk[each]  # Map k-core index to average degree for the bottom curve
            
    return  top_value,  bottom_value         

def MapPc2HeatMap(pc, avgk2index, kcore):
    '''
    Map percolation threshold values to heatmap coordinates.
    
    Parameters
    ----------
    pc : dict
        Dictionary containing critical thresholds for different k-core values and zeta values.
    avgk2index : dict
        Dictionary mapping average degrees to corresponding indices.
    kcore : int
        The k-core value.
    
    Returns
    -------
    x : list
        List of x-coordinates for the heatmap.
    y : list
        List of y-coordinates for the heatmap.
    '''
    
    # Get percolation thresholds for k-core 5
    kcore5_pc = pc[kcore]
    
    # Initialize dictionary to store heatmap indices for percolation thresholds
    PcshowInx = {}
    # Iterate over zeta values and percolation thresholds for k-core 5
    for zeta in kcore5_pc.keys():
        # Calculate the index of the percolation threshold in the heatmap
        PcshowInx[int(zeta-3)] =avgk2index[np.floor(kcore5_pc[zeta]*100)/10]
        
    # Convert dictionary keys to list of x-coordinates
    x = list(PcshowInx.keys())
    x.append(197) # Add an extra point for plotting purposes
    
    # Convert dictionary values to list of y-coordinates
    y = list(PcshowInx.values())
    y.append(y[-1]) # Add an extra point for plotting purposes
    
    return x, y

def PlotHeatMap(resultpath2, ax, kcores_rhc_array, avgk2index, zetas, kcore):
    '''
    Plot a heatmap of k-core percolation results.
    
    Parameters
    ----------
    resultpath2 : str
        Path to the folder containing the critical threshold of kcore percolation.
    ax : matplotlib.axes.Axes
        The axes to plot the heatmap on.
    kcores_rhc_array : numpy.ndarray
        2D array containing k-core percolation results.
    avgk2index : dict
        A mapping from average k values to array indices.
    zetas : list
        List of zeta values.
    kcore : int
        The k-core value.
    
    Returns
    -------
    None
    '''
    
    # Build a mapping from index to average k value
    index2avgk= {each[1]:each[0] for each in avgk2index.items()}
    
    # Load the critical threshold of k-core percolation
    pc = kp.load(resultpath2 + '/kcore_zeta_threshold_pc')

    # Plotting parameters
    cbars = False
    color2 ='white'
    cmap = "turbo"
    
    # Identify the curve
    [kcore_bottom, kcore_top]= Identify_curve(kcores_rhc_array)
    
    # Plot the heatmap
    h = sns.heatmap(kcores_rhc_array, ax=ax, vmin=0, vmax=500, cbar=cbars, cmap=cmap)
    cbar = plt.colorbar(h.collections[0], ax=ax, location='right', pad=0.02, orientation='vertical')
    cbar.set_label('$r^c_h$')  # Set color bar label
    
    # Set yticks and yticklabels
    yticks = np.arange(0,len(avgk2index.keys()),5)
    yticklabels = [8.5, 8.0, 7.5, 7.0, 6.5, 6.0]
    ax.set_yticks(yticks+0.5)
    ax.set_yticklabels(yticklabels)
        
    # Set xticks and xticklabels
    xticks = np.arange(2,100,15)
    xticklabels = zetas[xticks]
    ax.set_xticks(xticks+0.5)
    ax.set_xticklabels(xticklabels)
    
    # Plot kcore = 5 for the top
    top5_value, bottom5_value = Tran2avgK(kcore_top, kcore_bottom, index2avgk, zetas)        
    top5_x = np.array(list(top5_value.keys()))
    top5_y = np.array(list(top5_value.values()))
    top5y_index = np.array([avgk2index[each]+1 for each in top5_y]+[0.5])
    top5x_index = np.array([each-2.5 for each in top5_x]+[49.5])
    inx5_top = [0,1,2,3,7,11,16,21,32,40,49]  # Indices for plotting
    ax.plot(top5x_index[inx5_top], top5y_index[inx5_top], '--', color=color2)
    
    # Plot kcore = 5 for the bottom
    [x,y]= MapPc2HeatMap(pc, avgk2index, kcore=kcore)
    x[0] = 0.5
    index  = [0,1,2,3,4,10,11,12,13,14,15,16,17]
    ax.plot(np.array(x)[index], np.array(y)[index], '--', color=color2)
    
    # Add text annotations
    ax.text(30, 23, 'unstable', color ='white', size=12)
    ax.text(30, 12, 'metastable', color ='white', size=12)
    ax.text(5.8, 6, 'stable', color ='white', size=12, rotation=90)
    
    # Plot axes and add annotation
    kp.PlotAxes(ax, r"$\zeta$",  r"average degree, $\langle k \rangle$", '(a)')

def PLotFiniteSize(resultpath4, ax3):
    '''
    Plot finite size scaling curves.
    
    Parameters
    ----------
    resultpath4 : str
        Path to the directory containing result files.
    ax3 : AxesSubplot
        Matplotlib subplot to plot the finite size scaling curves.
    
    Returns
    -------
    None
    '''
    # Define average degrees
    avgks = [7.0, 7.1, 7.2, 7.3, 7.4]
    
    # Load data dictionary containing average degree and corresponding rhc values
    avgks_dict = kp.load(resultpath4+'/avgk_rhc.pkl')
    
    # Set the parameter of figure
    mfc = "None"
    colors = plt.get_cmap("Paired")
    tiltes = ['(c)']
    zeta = 18
    fontsize = 10
    
    # Plot the finite size scaling curves for each average degree
    for j, avgk in enumerate(avgks):
        # Extract rhc values for the current average degree
        rhc = avgks_dict[avgk]
        # Plot rhc values
        ax3.plot(rhc.keys(), rhc.values(), 'o-', mfc=mfc, color = colors(2*j+1), label = r'$\langle k \rangle$='+f'{avgk}')
        # Set subplot axes labels and title
        kp.PlotAxes(ax3, '$L$', r'$r_h^{c}$', tiltes[0], mode=False)
        
    # Add text indicating zeta value
    ax3.text(1500, 150, r"$\zeta$="+str(zeta), size = fontsize)    
    ax3.legend(bbox_to_anchor=(0.0,0.72), loc='center left', framealpha=0)
    # Set y-axis limits and ticks
    ax3.set_ylim(0,180)
    ax3.set_yticks([0,50,100,150])
    
def PlotZetaSystemSzie(rootpath, ax5):
    '''
    Plot system size scaling curves for different zeta values.
    
    Parameters
    ----------
    rootpath : str
        Path to the directory containing result files.
    ax5 : AxesSubplot
        Matplotlib subplot to plot the system size scaling curves.
    
    Returns
    -------
    None
    '''

    # Load the dataset for different zeta values
    zetas = [10, 20, 30, 50] 
    zetas_dict = kp.load(rootpath+'/zeta_rhc.pkl') # Load data dictionary containing zeta and corresponding rhc values
    
    # Set the parameter of figure
    fontsize = 10
    mfcs = 'None'
    colors = plt.get_cmap('Paired')
    
    # Plot system size scaling curves for each zeta value
    for i, zeta in enumerate(zetas):
        ax5.plot(zetas_dict[zeta].keys(), zetas_dict[zeta].values(), 'o-', mfc=mfcs, color=colors(2*i+1), label =r'$\zeta$='+f'{zeta}')
        
    # Set y-axis limits and add text annotation
    ax5.set_ylim(0,580)
    ax5.text(1500, 490, r"$\langle k \rangle$="+str(7.5), size = fontsize)    
    ax5.legend(bbox_to_anchor=(0.0,0.70), loc='center left', framealpha=0)
    
    # Set subplot axes labels and title
    kp.PlotAxes(ax5, r'$L$', r'$r_h^{c}$', '(d)', mode=False)
    
if __name__ == "__main__":
    
    # Set the result, network, and figure paths
    resultpath2 = '../kcorePercolation/figure4/result/result2'
    resultpath3 = '../kcorePercolation/figure4/result/result3'
    resultpath4 = '../kcorePercolation/figure4/result/result4'
    figurepath = '../kcorePercolation/figure4/figure'
    
    # Set the parameter
    zetas = np.arange(3,101,1)
    avgks = [6.0, 6.1, 6.2, 6.3, 6.4, 6.5, 6.6, 6.7, 6.8, 6.9, 7.0, 7.1, 7.2, 7.3, 7.4, 7.5, 7.6, 7.7, 7.8, 7.9, 
             8.0, 8.1, 8.2, 8.3, 8.4, 8.5] #5.0, 5.1, 5.2, 5.3, 5.4, 5.5, 5.5, 5.6, 5.7, 5.8, 5.9, 8.6, 8.7, 8.8, 8.9, 9.0, 9.1, 9.2, 9.3, 9.4, 9.5, 9.6, 9.7, 9.8, 9.9, 10.0
    
    # Load the datasets
    kcore = 5
    kcores_rhc = loadAvgkKcoreZetaRhc(resultpath3, kcore,avgks,zetas)
      
    # Transform dicts into the array
    [kcores_rhc_array, avgk2index, zeta2index] = Transdict2Array(kcores_rhc, avgks, zetas)
    
    #plot the result
    fig, ax = plt.subplots(2,2, figsize=(6.5, 6.5), tight_layout=True) 
    
    PlotSchematic1(ax[0,1], figurepath)
    PlotHeatMap(resultpath2, ax[0,0], kcores_rhc_array, avgk2index, zetas, kcore)
    PLotFiniteSize(resultpath4, ax[1,0])
    PlotZetaSystemSzie(resultpath4,ax[1,1])
    
    plt.savefig(figurepath+'/Fig4.png', dpi = 500)
    plt.savefig(figurepath+'/Fig4.pdf')
    plt.savefig(figurepath+'/Fig4.eps')