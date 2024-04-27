# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 20:15:15 2023

@author: Leyang XUe
"""

import sys 
import numpy as np
from multiprocessing import Pool
 
# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def GenerateZetaNet(args):
    """
    Generate a network with a specific zeta value.
    
    Parameters
    ----------
    args : list
        List containing parameters for network generation:
        i : int
            Network ID.
        avg_k : float
            Average degree.
        zeta : int
            Zeta value.
        M : int
            Number of edges.
        L : int
            Lattice length.
        networkpath : str
            Path to save the network.
    
    Returns
    -------
    None
    
    """
    
    i, avg_k, zeta, M, L, networkpath= args
    range_coordinate = np.arange(0,L,1) #range of coordinate 
    
    print('netId:%d, zeta:%d'%(i, zeta))    
    
    # Generate the sequence of links with length l
    links = kp.expvariate_sequence(zeta, M, L)
    networkpath_link = networkpath + '/' + 'NetID'+str(i)+'_avgk' + str(avg_k) + '_links_zeta'+str(zeta)+'.pkl'
    kp.save(networkpath_link, links)
 
    # Create the network 
    filename = 'NetID'+str(i)+'_avgk' + str(avg_k) + '_zeta'+str(zeta)
    G = kp.NetworkCreate(range_coordinate, links, networkpath, filename)
    
def GNParallel(networkpath, i, zetas, avg_k, M, L):
    """
    Generate networks with different zeta values in parallel.
    
    Parameters:
    -----------
    networkpath : str
        Path to save the generated networks.
    i : int
        Network ID.
    zetas : list
        List of zeta values for which networks need to be generated.
    avg_k : float
        Average degree.
    M : int
        Number of edges.
    L : int
        Lattice length.
    """
    
    args = []
    # Create arguments for GenerateZetaNet function
    for zeta in zetas:
        args.append([i, avg_k, zeta, M, L, networkpath])  
        
    # Create a pool of worker processes
    pool = Pool(processes=int(len(zetas)))    
    # Distribute the tasks to the worker processes
    pool.map(GenerateZetaNet,args)
    # Close the pool to prevent further tasks submission
    pool.close()

def CalculateRhC(args):
    """
    Calculate the critical radius of the hole (Rh_c) for a given zeta and kcore.
    
    Parameters:
    -----------
    args : list
        List containing the following elements in order:
            i : int
                Simulation ID.
            ResultPath : str
                Path to save the results.
            networkpath : str
                Path where the network files are stored.
            avg_k : float
                Average degree.
            zeta : int
                Zeta value.
            kcore : int
                K-core value.
            Rh_low : float
                Lower bound of the critical radius.
            Rh_high : float
                Upper bound of the critical radius.
            threshold : float
                Convergence threshold for binary search.
    
    Returns:
    --------
    tuple
        A tuple containing the zeta value and the calculated critical radius (Rh_c).
    """
    
    # Unpack the arguments
    i, ResultPath, networkpath, avg_k, zeta, kcore, Rh_low, Rh_high, threshold = args 
    
    # Load the network
    G = kp.load(networkpath+'/NetID0_avgk'+str(avg_k)+'_zeta'+str(zeta)+'_spatialNet.pkl') 
    # Calculate the critical radius of the hole (Rh_c)
    Rh_c = kp.BinarySearch(G, kcore, Rh_low, Rh_high, threshold)
    
    # Save the result
    result = (zeta,Rh_c)        
    kp.save(ResultPath+'/avgk'+str(avg_k)+'_kcore'+str(kcore)+'_zeta'+str(zeta)+'_simulation'+str(i)+'.pkl', result)
    
    return zeta, Rh_c

def RhcParallel(ResultPath, networkpath, zetas, avg_k, kcore, Rh_low, Rh_high, threshold, simulation):
    """
    Calculate the critical hole radius in parallel for multiple zeta values.
    
    Parameters
    ----------
    ResultPath : str
        Path to save the results.
    networkpath : str
        Path to load the network.
    zetas : list
        List of zeta values.
    avg_k : float
        Average degree of the network.
    kcore : int
        K-core value.
    Rh_low : float
        Lower bound for critical radius.
    Rh_high : float
        Upper bound for critical radius.
    threshold : float
        Threshold value for binary search.
    simulation: int
        Time of simution
        
    Returns
    -------
    None
    
    """    
    
    args = []
    
    # Calculate the Rhc
    for zeta in zetas:
        print(f'Began to run for zeta:{zeta}')
        for i in np.arange(simulation):
            args.append([i, ResultPath, networkpath, avg_k, zeta, kcore, Rh_low, Rh_high, threshold])        

        pool = Pool(processes=int(simulation))    
        results = pool.map(CalculateRhC,args)
        pool.close()
        kp.save(ResultPath+'/avgk'+str(avg_k)+'_kcore'+str(kcore)+'_zeta'+str(zeta)+'.pkl', results)  
    
if __name__ == "__main__":
    
    # Set the network and result paths
    networkpath = '../kcorePercolation/figure4/network'
    resultpath = '../kcorePercolation/figure4/result/result'
 
    # Set parameters for creating the network
    L = 300 # [500,700,1000,1500,2000]
    N = L*L  # The number of nodes 
    i = 0 # Simulation times
    kcores = [5]
    
    # Define average degrees and zetas to iterate over
    avg_ks = [round(each,2) for each in np.arange(7.5, 7.52, 0.1)]
    simulation = 10
    zetas = [10,20,30,50] #characteristic length #80 
    
    # Loop over each average degree
    for avg_k in avg_ks:
         
         print('Begin to run for average degree:', avg_k)
         # Update network path for current average degree
         networkpath = networkpath+f'/size_dependency/L{L}/avgk'+str(avg_k)
         kp.mkdirectory(networkpath) #create the folder
         
         # Calculate the number of edges (M) based on average degree
         M = int((N * avg_k)/2) #the number of edges
         GNParallel(networkpath, i, zetas, avg_k, M, L)
         print('Finish creating network')
         
         for kcore in kcores:
             
             print('Begin calculating for k-core:', kcore)
             resultpath = resultpath+f'/L{L}/avgk'+str(avg_k)+'_kcore'+str(kcore)
             kp.mkdirectory(resultpath)

             # Identify the critical radius of hole    
             Rh_low = 0
             Rh_high = L/2
             threshold = 0.00001
             RhcParallel(resultpath, networkpath, zetas, avg_k, kcore, Rh_low, Rh_high, threshold, simulation)
    
