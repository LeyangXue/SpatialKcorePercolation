# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 19:43:01 2024

@author: Leyang Xue
"""

import sys
import numpy as np
from multiprocessing import Pool 
import random 
import networkx as nx

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def kcorepercolation(G,p,k):
    '''
    Perform k-core percolation on a graph G with node ramining probability p and k-core value k.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    p : float
        The node ramining probability.
    k : int
        The k-core value.

    Returns
    -------
    lcc : float
        The size of the largest connected component normalized by the total number of nodes.

    '''
    
    #obtain the removed nodes 
    rm_nodes_num = int((1-p)*G.order())
    nodes = list(G.nodes())
    rm_nodes = random.sample(nodes, rm_nodes_num)
    
    #remove the nodes 
    Gr = G.copy()
    N = G.order()
    Gr.remove_nodes_from(rm_nodes)
    degree_dict = {node: degree for node, degree in Gr.degree()}

    while len(degree_dict)!=0 and min(degree_dict.values()) < k:
                
        #remove nodes with degree less than k
        nodes = list(Gr.nodes())
        for node in nodes:
            #print('node',node)
            if degree_dict[node] < k:
                Gr.remove_node(node)
        
        #update the degree dict
        degree_dict = {node: degree for node, degree in Gr.degree()}
    
    #calculate the connected components
    connect_components = sorted(nx.connected_components(Gr), key=len, reverse=True)
    if len(connect_components) > 0:
        lcc = len(connect_components[0])/N 
    else:
        lcc = 0
        
    return lcc

def find_critical_point(G, p_min, p_max, kcore, threshold):
    '''
    Find the critical point for k-core percolation on a graph G.

    Parameters
    ----------
    G : networkx.Graph
        The input graph.
    p_min : float
        The minimum remaining probability of nodes.
    p_max : float
        The maximum remaining probability of ndoes.
    kcore : int
        The k-core value.
    threshold : float
        The threshold value.

    Returns
    -------
    p_max : float
        The critical point.

    '''
    
    precision= 1 / G.order() 
    while abs(p_max - p_min) > precision:
        p_mid = (p_min + p_max) / 2
        pinfty = kcorepercolation(G, p_mid, kcore)
        print(f'current pmin={p_min}, pmax={p_max}, p_mid = {p_mid} and pinfty={pinfty}')
        
        if pinfty < threshold:
            p_min = p_mid
        else:
            p_max = p_mid
            
    return p_max

def BinarySearchPc(args):
    '''
    Perform a binary search to find the critical point for k-core percolation.

    Parameters
    ----------
    args : list
        List of arguments including avgk, zeta, p_min, p_max, kcore, threshold, itr, networkpath.

    Returns
    -------
    pc : float
        The critical point.
    
    '''
    
    avgk, zeta, p_min, p_max, kcore, threshold, itr, networkpath = args
    G = kp.load(networkpath+f'/NetID0_avgk{avgk}_zeta{zeta}_spatialNet.pkl')
    print(f'zeta:{zeta} with simulation:{itr}')
    
    #find the critical point
    pc = find_critical_point(G, p_min, p_max, kcore, threshold)

    return pc 

def PercolationPcParalell(resultpath, networkpath, kcore, avgk, zeta, p_min, p_max, threshold, simulation_num):
    '''
    Perform parallel binary search to find the critical point for k-core percolation.

    Parameters
    ----------
    resultpath : str
        Path to save the result files.
    networkpath : str
        Path to the network files.
    kcore : int
        The k-core value.
    avgk : float
        Average degree of the network.
    zeta : float
        Zeta parameter.
    p_min : float
        The minimum remaining probability of nodes.
    p_max : float
        The maximum remaining probability of nodes.
    threshold : float
        The threshold value.
    simulation_num : int
        Number of simulations.

    Returns
    -------
    None.

    '''
    args = []    
    #perform the percolation process
    for itr in np.arange(simulation_num):
        args.append([avgk, zeta, p_min, p_max, kcore, threshold, itr, networkpath])
    
    # The number of processes should be chosen based on the number of CPU cores available
    # Please adjust the number of processes according to the actual machine configuration   
    # You can change the divisor in the following line to adjust the number of processes
    pool = Pool(processes=int(simulation_num/15))
    results = pool.map(BinarySearchPc,args)
    pool.close()

    #save the result for a given zeta        
    kcore_path = resultpath+f'/kcore{kcore}'
    kp.mkdirectory(kcore_path)
    savepath= kcore_path+f'/kcore{kcore}_zeta{zeta}_binarysearchpc.pkl'
    kp.save(savepath, results)

if __name__ == "__main__":
    
    # Set the paths
    resultpath = '../kcorePercolation/figure2/result' 
    networkpath = '../kcorePercolation/figure2/network' 
    # Please notice: download and put the network structure into the current path or run the code to generate the network structure before running this code

    # Set the parameters for k-core percolation analysis
    kcore = 5
    zetas = [5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,50,100,500,1000]#
    avg_k = 10  # Average degree of the network
    simu = 50   # Number of simulations
    p_window = 0.1  # Window size for defining the range of node reamining probabilities
    threshold = 0.001  # Threshold value for k-core percolation

    # Load pre-computed critical points
    pc = kp.load(resultpath+'/kcore_zeta_threshold_pc')
    kcore_pc = pc[kcore]
    print('kcore pc', kcore_pc)

    # Perform k-core percolation analysis for each zeta value
    for zeta in zetas:

        # Set the probability range based on pre-computed critical points
        if zeta in [500,1000]:
            zeta_pc = pc[kcore][200] 
        elif zeta in [16,17,18,19]:          
            zeta_pc = pc[kcore][15] 
        else:
            zeta_pc = pc[kcore][zeta]
            
        # Define the range of node reamining probabilities
        p_min = zeta_pc - p_window
        p_max = zeta_pc + p_window

        # Run k-core percolation analysis in parallel
        print(f'Starting k-core percolation analysis for zeta={zeta} with p_min={p_min} and p_max={p_max}')
        PercolationPcParalell(resultpath, networkpath, kcore, avg_k, zeta, p_min, p_max, threshold, simu)
        print(f'Completed k-core percolation analysis for zeta={zeta}')
