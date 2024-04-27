# -*- coding: utf-8 -*-
"""
Created on Mon Jul 10 14:49:57 2023

@author: Leyang Xue
"""

import sys
import numpy as np
import os

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def CalculateAvgkRhc(resultpath, avgks, Ls, savepath):
    """
    Calculate the average value of Rhc for different avgk and lattice sizes (Ls).

    Parameters:
    -----------
    resultpath : str
        Path to the folder containing the result files.
    avgks : list
        List of average degree values.
    Ls : list
        List of lattice sizes.
    savepath : str
        Path to save the calculated average Rhc values.

    Returns:
    --------
    dict
        A nested dictionary containing the average Rhc values for different avgk and lattice sizes.
    """
    # Dictionary to store the average Rhc values for different avgk and lattice sizes
    avgk_rhc = {} 
    
    # Iterate over each avgk value
    for avgk in avgks:
        L_rhc = {}
        
        # Iterate over each lattice size (L)
        for L in Ls:
            simulation = []
            file_folder = resultpath + f'/L{L}/avgk{avgk}_kcore5'
            files = os.listdir(file_folder)
            
            # Load the Rhc values from result files
            for file in files:
                rhc = kp.load(file_folder + '/' + file)
                simulation.append(rhc[1])
            
            # Calculate the average Rhc value for the current avgk and lattice size
            rhc_avg = np.average(np.array(simulation))
            L_rhc[L] = rhc_avg
        
        # Store the average Rhc values for the current avgk
        avgk_rhc[avgk] = L_rhc
    
    # Save the calculated average Rhc values
    kp.save(savepath + '/avgk_rhc.pkl', avgk_rhc)

    return avgk_rhc

def CalculateZetaRhc(resultpath, zetas, Ls, savepath):
    """
    Calculate the average value of Rhc for different zeta values with a given k.

    Parameters:
    -----------
    resultpath : str
        Path to the folder containing the result files.
    zetas : list
        List of zeta values.
    Ls : list
        List of lattice sizes.
    savepath : str
        Path to save the calculated average Rhc values.

    Returns:
    --------
    dict
        A nested dictionary containing the average Rhc values for different zeta values and lattice sizes.
    """
    # Dictionary to store the average Rhc values for different zeta values and lattice sizes
    zetas_rhc = {}
    
    # Iterate over each zeta value
    for zeta in zetas:
        L_rhc = {}
        
        # Iterate over each lattice size (L)
        for L in Ls:
            simulation = []
            file_folder = resultpath + f'/L{L}/avgk7.5_kcore5/zeta{zeta}'
            files = os.listdir(file_folder)
            
            # Load the Rhc values from result files
            for file in files:
                rhc = kp.load(file_folder + '/' + file)
                simulation.append(rhc[1])
            
            # Calculate the average Rhc value for the current zeta and lattice size
            rhc_avg = np.average(np.array(simulation))
            L_rhc[L] = rhc_avg
        
        # Store the average Rhc values for the current zeta
        zetas_rhc[zeta] = L_rhc
    
    # Save the calculated average Rhc values
    kp.save(savepath + '/zeta_rhc.pkl', zetas_rhc)

    return zetas_rhc

if __name__ == "__main__":
    
    # Set the network and result paths
    resultpath = '../kcorePercolation/figure4/result/result'
    savepath = '../kcorePercolation/figure4/result/result4'
    
    Ls = [300,500,700,1000,1500,2000]
    avgks = [7.0, 7.1, 7.2, 7.3, 7.4]
    zetas = [10,20,30,50,80]

    # Calculate the result for different avgk with zeta = 18 
    avgk_rhc = CalculateAvgkRhc(resultpath, avgks, Ls, savepath)
    
    # Calculate the result for different zeta with avgk = 7.5 
    zetas_rhc = CalculateZetaRhc(resultpath, zetas, Ls, savepath)