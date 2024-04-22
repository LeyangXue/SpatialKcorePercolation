# -*- coding: utf-8 -*-
"""
Created on Sun Mar 5 16:18:59 2023

@author: LeyangXue

"""
import sys
import os
import numpy as np
from multiprocessing import Pool 

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def GenerateZetaNet(args):
    """
    Generate network for specific parameters

    Args:
        args (list): List of parameters including network ID (i), average degree (avg_k),
                     zeta value (zeta), number of edges (M), lattice length (L), and network save path (networkpath)

    Returns:
        None
    """
    i, avg_k, zeta, M, L, networkpath = args
    range_coordinate = np.arange(0, L, 1)  # Coordinate range
    
    print('Generating network for netId:%d, zeta:%d' % (i, zeta))    
    # Generate sequence of links with specified length and parameters             
    links = kp.expvariate_sequence(zeta, M, L)
    networkpath_link = os.path.join(networkpath, f'NetID{i}_avgk{avg_k}_links_zeta{zeta}.pkl')
    kp.save(networkpath_link, links)
 
    # Create network 
    filename = f'NetID{i}_avgk{avg_k}_zeta{zeta}'
    kp.NetworkCreate(range_coordinate, links, networkpath, filename)

def GNParallel(networkpath, i, zetas, avg_k, M, L):
    '''
        Generate networks in parallel.

    Parameters
    ----------
    networkpath : str
             Path to save the networks.
    i : int
        Network ID.
    zetas : list
        List of zeta values.
    avg_k : int
        Average degree.
    M : int
        Number of edges.
    L : int
        Lattice length.

    Returns
    -------
    None.
    '''
   
    args = []
    for zeta in zetas:
        args.append([i, avg_k, zeta, M, L, networkpath])   
    
    # Create a pool of worker processes for parallel execution
    # The number of processes is set to one-fifth of the length of the zetas list
    # This ensures efficient utilization of system resources while parallelizing the network generation process
    pool = Pool(processes=int(len(zetas)/5))    
    pool.map(GenerateZetaNet, args)
    pool.close()
       
if __name__ == "__main__":
    
    # Set root path
    rootpath = '../kcorePercolation/figure1' # Manually specify the path to your script
    
    # Create network 
    L = 1000  # Lattice length
    N = L * L  # Number of nodes
    
    # Define range of zeta values
    zetas = list(np.arange(3, 21)) + [30, 40, 60, 70, 80, 90, 100, 500, 1000]
    avg_k = 10
    i = 0

    print('Begin to generate the network...') 
    networkpath = os.path.join(rootpath, 'network')
    kp.mkdirectory(networkpath)  # Create folder to save networks
    M = int((N * avg_k) / 2)  # Calculate number of edges
    GNParallel(networkpath, i, zetas, avg_k, M, L)
    print('Finish generating the network.')    
