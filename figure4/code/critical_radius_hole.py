# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 16:18:59 2023

@author: Leyang Xue

"""
import sys
import numpy as np
from multiprocessing import Pool 
import datetime

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def GenerateZetaNet(args):
    """
    Generate network for a specific zeta value.

    Parameters
    ----------
    args : list
        List containing the parameters:
            i : int
                Simulation index.
            avg_k : float
                Average degree of the network.
            zeta : int
                Zeta value for characteristic length.
            M : int
                Number of edges in the network.
            L : int
                Lattice length.
            networkpath : str
                Path to save the generated networks.

    Returns
    -------
    None

    """
    
    i, avg_k, zeta, M, L, networkpath= args
    range_coordinate = np.arange(0,L,1) # Range of coordinates 
    
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
    Generate networks in parallel for different zeta values.
    
    Parameters
    ----------
    networkpath : str
        Path to save the generated networks.
    i : int
        Simulation times.
    zetas : list
        List of zeta values for characteristic length.
    avg_k : float
        Average degree of the network.
    M : int
        Number of edges in the network.
    L : int
        Lattice length.
    
    Returns
    -------
    None
    
    """
    
    args = []
    for zeta in zetas:
        args.append([i, avg_k, zeta, M, L, networkpath])   
       
    pool = Pool(processes=int(len(zetas)/5))    
    pool.map(GenerateZetaNet,args)
    pool.close()
       
def CalculateRhC(args):
    """
    Calculate the critical hole radius for a given zeta and kcore.
    
    Parameters
    ----------
    args : list
        List containing the parameters:
            ResultPath : str
                Path to save the results.
            networkpath : str
                Path to load the network.
            avg_k : float
                Average degree of the network.
            zeta : int
                Zeta value for characteristic length.
            kcore : int
                K-core value.
            Rh_low : float
                Lower bound for critical radius.
            Rh_high : float
                Upper bound for critical radius.
            threshold : float
                Threshold value for binary search.
    
    Returns
    -------
    tuple
        Tuple containing the zeta value and the calculated critical hole radius.
    
    """
    
    # For a given zeta and kcore
    ResultPath, networkpath, avg_k, zeta, kcore, Rh_low, Rh_high, threshold = args 
    
    # Load the network
    G = kp.load(networkpath+'/NetID0_avgk'+str(avg_k)+'_zeta'+str(zeta)+'_spatialNet.pkl')   
    # Calculate the critical hole radius
    Rh_c = kp.BinarySearch(G, kcore, Rh_low, Rh_high, threshold)
    
    result = (zeta,Rh_c)        
    kp.save(ResultPath+'/avgk'+str(avg_k)+'_kcore'+str(kcore)+'_zeta'+str(zeta)+'.pkl', result)
    
    return zeta, Rh_c

def RhcParallel(resultPath, networkpath, zetas, avg_k, kcore, Rh_low, Rh_high, threshold):
    """
    Calculate the critical hole radius in parallel for multiple zeta values.
    
    Parameters
    ----------
    resultPath : str
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
    
    Returns
    -------
    None
    
    """    
    
    args = []
    
    # Calculate the Rhc
    for zeta in zetas:    
        args.append([resultPath, networkpath, avg_k, zeta, kcore, Rh_low, Rh_high, threshold])        

    pool = Pool(processes=int(len(zetas)/5))    
    results = pool.map(CalculateRhC,args)
    pool.close()
    kp.save(resultPath+'/avgk'+str(avg_k)+'_kcore'+str(kcore)+'.pkl', results)  
    
    # Transform into a dictionary
    Rhc_dict = {each[0]:each[1] for each in results}
    kp.save(resultPath+'/avgk'+str(avg_k)+'_kcore'+str(kcore)+'_dict.pkl', Rhc_dict)  

if __name__ == "__main__":
    
    # Set the network and result paths
    networkpath = '../kcorePercolation/figure4/network'
    resultpath = '../kcorePercolation/figure4/result/result3'

    # Set parameters for creating the network
    L = 1000 # Lattice length
    N = L*L   # Number of nodes
    i = 0  # Simulation times
    kcores = [3,4,5,6] # K-core values to consider
    
    # Define average degrees to iterate over
    avg_ks = [round(each,2) for each in np.arange(5.0, 10.01, 0.1)]
    # Define zeta values for characteristic length
    zetas = np.arange(3,101,1)
    
    # Loop over each average degree
    for avg_k in avg_ks:
        
         print('Begin to run for average degree:', avg_k)
         start_time=datetime.datetime.now()
         
         # Update network path for current average degree
         networkpath_avgk = networkpath+'/avgk'+str(avg_k)
         kp.mkdirectory(networkpath) # Create the folder if it doesn't exist
         
         # Calculate the number of edges (M) based on average degree
         M = int((N * avg_k)/2) 
         # Run GNParallel to generate networks for various zeta values
         GNParallel(networkpath_avgk, i, zetas, avg_k, M, L)
         
         end_time=datetime.datetime.now()
         print('Finish creating network, time taken:', end_time - start_time)
         
         # Loop over each k-core value
         for kcore in kcores:
             
             print('Begin calculating for k-core:', kcore)
             start_time=datetime.datetime.now()
             
             # Define the path to save results for the current k-core and average degree
             savepath = resultpath+f'/L{L}/avgk'+str(avg_k)+'_kcore'+str(kcore)
             kp.mkdirectory(savepath)

             # Identify the critical radius of the hole using RhcParallel
             Rh_low = 0
             Rh_high = L/2
             threshold = 0.0001 #depend on the system size
             RhcParallel(savepath, networkpath_avgk, zetas, avg_k, kcore, Rh_low, Rh_high, threshold)
             
             end_time=datetime.datetime.now()
             print('Finish calculating k-core', kcore, 'time taken:', end_time - start_time)
    
