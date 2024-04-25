# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 19:44:26 2024

@author:  Leyang Xue
"""

import sys
import networkx as nx

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def kcorepercolationPlateau(G, rm_nodes, k, p, z, simu, resultpath):
    '''
    Perform k-core percolation and save the largest connected component (LCC) at each step.

    Parameters
    ----------
    G : networkx.Graph
        The spatial network.
    rm_nodes : list
        List of nodes to remove.
    k : int
        The k-core value.
    p : float
        Probability of node reserved.
    z : int
        Zeta value for spatial network generation.
    simu : int
        Simulation number.
    resultpath : str
        Path to save result files.

    Returns
    -------
    lcc : set
        Nodes in the largest connected component after stabilization.
    NOI : int
        Number of iterations until k-core percolation stabilizes..
    '''
    # Make a copy of the graph
    Gr = G.copy()
    N = G.order() # Total number of nodes in the network
    
    # Initialize variables
    NOI = 0  # Number of iterations until stabilization
    lccs_step_dict = {NOI:1} # Dictionary to track LCC fraction over steps
    
    # Save initial LCC
    initial_lcc  = sorted(nx.connected_components(Gr), key=len, reverse=True)[0]
    kp.save(resultpath+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', initial_lcc)

    # Remove nodes and update the degree dictionary
    Gr.remove_nodes_from(rm_nodes)
    degree_dict = {each[0]:each[1] for each in Gr.degree()}
    
    # Perform k-core percolation until stabilization
    while len(degree_dict)!=0 and min(degree_dict.values()) < k:
        # Increment the step count
        NOI += 1
        
        # Save the LCC at each step
        lcc = sorted(nx.connected_components(Gr), key=len, reverse=True)[0]
        lccs_step_dict[NOI] = len(lcc)/N
        kp.save(resultpath+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', lcc)

        # Remove nodes with degree less than k
        nodes_to_remove = [node for node, degree in degree_dict.items() if degree < k]
        Gr.remove_nodes_from(nodes_to_remove)
        
        # Update the degree dictionary
        degree_dict = {node: degree for node, degree in Gr.degree()}
    
    # Calculate the connected components after stabilization
    connect_components = sorted(nx.connected_components(Gr), key=len, reverse=True)
    if connect_components:
        lcc = connect_components[0]
        lcc_p = len(lcc) / N 
        if lcc_p != lccs_step_dict[NOI]:
            NOI  += 1  
            lccs_step_dict[NOI] = lcc_p 
            kp.save(resultpath+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', lcc)
    else:
        lcc = 0
        NOI  += 1  
        lccs_step_dict[NOI] = lcc
        kp.save(resultpath+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', connect_components)

    # Save the results
    savepath_lcc_num = resultpath+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_steps_lcc_num_dict.pkl'
    kp.save(savepath_lcc_num, lccs_step_dict)

    return lcc, NOI

if __name__ == "__main__":

    # Set the paths for result, network, and gcc files
    resultpath = '../kcorePercolation/figure2/result' 
    networkpath = '../kcorePercolation/figure2/network'
    gccpath =  '../kcorePercolation/figure2/gcc'

    # Set basic information about the spatial networks
    kcore = 5
    zetas = [4, 6, 7, 8, 9, 10, 12, 14, 16, 18, 20, 50, 60, 70, 80, 100, 500] 
    avg_k = 10 #average degree
    simu = 1
    
    # Define the probability of reserved nodes at criticality for different zeta values
    zeta_pc_removed_nodes = {4:(0.701501, 50), 6:(0.721001,23), 7:(0.718501,9), 8:(0.710351,18), 9:(0.709751,30), 10:(0.700501,40), 12:(0.69765,27), 14:(0.692701,34), 16:(0.691051,36), 18:(0.687351,73), 20:(0.6876,28),
    50:(0.681601,20), 60:(0.679701,39), 70:(0.680201,34), 80:(0.681001,26), 100:(0.680601,25), 500:(0.679151,37)}
    
    # Iterate over different zeta values
    for zeta in zetas:
        print(f'run the task with zeta={zeta}')    
        
        # Get the probability to reserve the nodes at the criticality 
        (zeta_pc,step) = zeta_pc_removed_nodes[zeta]

        # Load the removed nodes at the criticality for the current zeta
        rm_nodes = kp.load(gccpath+f'/kcore{kcore}/rmnodes/kcore{kcore}_zeta{zeta}_simu{simu}_p{zeta_pc}_step{step}_rmnodes.pkl')
            
        # Load the network
        G = kp.load(networkpath+f'/NetID0_avgk{avg_k}_zeta{zeta}_spatialNet.pkl')
        
        # Set the path to save connected components at each step
        savegccpath = gccpath+f'/kcore{kcore}/zeta{zeta}_lcc_nodes/'
        kp.mkdirectory(savegccpath)
        
        # Run k-core percolation for the current zeta
        [lcc, NOI] =  kcorepercolationPlateau(G, rm_nodes, kcore, zeta_pc, zeta, simu, savegccpath)
        print('zeta:',zeta,'lcc:',lcc,'NOI:',NOI)