# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 19:44:26 2024

@author:  Leyang Xue
"""

import sys
import networkx as nx

# Add package path
packagepath = 'F:/work/work12/kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def kcorepercolationPlateau(G, rm_nodes, k, p, z, simu, savegccpath, resultpath):
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
    savegccpath : str
        Path to save the lcc nodes at each step
    resultpath : str
        Path to save the number of lcc nodes at each step

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
    kp.save(savegccpath+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', initial_lcc)

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
        kp.save(savegccpath+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', lcc)

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
            kp.save(savegccpath+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', lcc)
    else:
        lcc = 0
        NOI  += 1  
        lccs_step_dict[NOI] = lcc
        kp.save(savegccpath+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', connect_components)

    # Save the results
    savepath_lcc_num = resultpath+f'/zeta{z}/kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_steps_lcc_num_dict.pkl'
    kp.save(savepath_lcc_num, lccs_step_dict)

    return lcc, NOI

if __name__ == "__main__":

    # Set the paths for result, network, and gcc files
    resultpath = 'F:/work/work12/kcorePercolation/figure3/result' 
    networkpath = 'F:/work/work12/kcorePercolation/figure3/network'

    # Set basic information about the spatial networks
    kcore = 5
    zetas = [7, 500] 
    avg_k = 10 #average degree
    simu_dict = {7:1, 500:2}
    
    # Define the probability of reserved nodes at criticality for different zeta values
    zeta_pc_removed_nodes = {7:(0.718501,9), 500:(0.679851,50)}
    
    # Iterate over different zeta values
    for zeta in zetas:
        print(f'run the task with zeta={zeta}')    
           
        # Get the probability to reserve the nodes at the criticality 
        (zeta_pc,step) = zeta_pc_removed_nodes[zeta]

        # Load the removed nodes at the criticality for the current zeta
        rm_nodes = kp.load(networkpath+f'/rmnodes/kcore{kcore}_zeta{zeta}_simu{simu_dict[zeta]}_p{zeta_pc}_step{step}_rmnodes.pkl')
            
        # Load the network
        G = kp.load(networkpath+f'/NetID0_avgk{avg_k}_zeta{zeta}_spatialNet.pkl')
        
        # Set the path to save connected components at each step
        savegccpath = networkpath+f'/zeta{zeta}_lcc_nodes/'
        kp.mkdirectory(savegccpath)
        
        # Run k-core percolation for the current zeta
        [lcc, NOI] =  kcorepercolationPlateau(G, rm_nodes, kcore, zeta_pc, zeta, simu_dict[zeta], savegccpath, resultpath)
        print('zeta:',zeta,'lcc:',lcc,'NOI:',NOI)