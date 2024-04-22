# -*- coding: utf-8 -*-
"""
Created on Sun Mar 5 16:18:59 2023

@author: LeyangXue

"""
import sys
import networkx as nx
import random 
from multiprocessing import Pool 

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def kcore_percolation_plateau(G,rm_nodes,k, p, z, simu, resultpath):
    '''
    Perform k-core percolation simulation

    Parameters
    ----------
    G : Graph
        Input graph.
    rm_nodes : list
        Nodes to remove initially.
    k : int
        k-core value.
    p : float
        Probability threshold for node removal.
    z : int
        Zeta value.
    simu : int
        Simulation number.
    resultpath : str
        Path to save results.

    Returns
    -------
    lcc : float
        Largest connected component size (normalized).
    NOI : int
        Number of iterations until plateau.
    '''
    
    #remove the nodes 
    Gr = G.copy()
    N = G.order()
    Gr.remove_nodes_from(rm_nodes)
    degree_dict = {each[0]:each[1] for each in Gr.degree()}

    #count the number of stages, second_lcc in the iterative process
    NOI = 0
    lccs_step_dict = {NOI:1}
    
    while len(degree_dict)!=0 and min(degree_dict.values()) < k:
        
        #update the steps and lcc of each steps
        NOI += 1
        lccs_num = len(sorted(nx.connected_components(Gr), key=len, reverse=True)[0])
        lccs_step_dict[NOI] = lccs_num/N
        
        #remove nodes with degree less than k
        nodes = list(Gr.nodes())
        for node in nodes:
            if degree_dict[node] < k:
                Gr.remove_node(node)
        
        #update the degree dict
        degree_dict = {each[0]:each[1] for each in Gr.degree()}
    
    #calculate the connected components
    connect_components = sorted(nx.connected_components(Gr), key=len, reverse=True)
    if len(connect_components) > 0:
        lcc = len(connect_components[0])/N 
        if lcc != lccs_step_dict[NOI]:
            NOI  += 1  
            lccs_step_dict[NOI] = lcc 
    else:
        lcc = 0
        NOI  += 1  
        lccs_step_dict[NOI] = lcc
    
    #create the path to save the result
    zeta_path = resultpath+f'/kcore{k}/zeta{z}/'
    kp.mkdirectory(zeta_path)
    
    savepath = zeta_path+f'kcore{k}_zeta{z}_p{round(p,7)}_simu{simu}_lcc_dict.pkl'
    kp.save(savepath, lccs_step_dict)
    
    return lcc, NOI

def run_kcore_percolation(avg_k, p_start, p_end, p_iterate, kcore, z, simu, resultpath):
    '''
    Perform k-core percolation simulation over a range of probabilities.

    Parameters
    ----------
    avg_k : int
        Average degree.
    p_start : float
        Starting probability for node removal.
    p_end : float
        Ending probability for node removal.
    p_iterate : float
        Probability step size.
    kcore : int
        k-core value.
    z : int
        Zeta value.
    simu : int
        Simulation number.
    resultpath : str
        Path to save results.

    Returns
    -------
    None.

    '''
    
    G = kp.load(networkpath + f'/NetID0_avgk{avg_k}_zeta{z}_spatialNet.pkl')

    N = G.order()
    initial_rm_nodes_num = int((1 - p_start) * N)
    end_rm_nodes_num = int((1 - p_end) * N)
    update_node_num = int(p_iterate * N)

    nodes = set(list(G.nodes()))
    rm_nodes = random.sample(nodes, initial_rm_nodes_num)
    rm_nodes_set = set()

    pinfty_dict = {}
    NOI_dict = {}

    while len(rm_nodes) <= end_rm_nodes_num:
        p = 1 - len(rm_nodes) / N
        lcc, NOI = kcore_percolation_plateau(G, rm_nodes, kcore, p, z, simu, resultpath)
        pinfty_dict[p] = lcc
        NOI_dict[p] = NOI

        print('run the kcore percolation for kcore={kcore}, zeta={z}, p={p}, simu={simu}, lcc={lcc}, NOI={NOI}'.format(kcore=kcore, z=z, p=round(p,7), simu=simu, lcc=lcc, NOI=NOI))

        rm_nodes_set.update(rm_nodes)
        rest_nodes_set = nodes - rm_nodes_set
        update_nodes = random.sample(rest_nodes_set, update_node_num)
        rm_nodes += update_nodes

    savepath = resultpath + f'/kcore{kcore}/zeta{z}/'
    kp.save(savepath + f'kcore{kcore}_zeta{z}_simu{simu}_pinfty_dict.pkl', pinfty_dict)
    kp.save(savepath + f'kcore{kcore}_zeta{z}_simu{simu}_NOI_dict.pkl', NOI_dict)

def parallel_computation(args):
    """
    Perform parallel computation of k-core percolation simulation.

    Args:
        args (list): List of arguments for simulation.

    Returns:
        None
    """
    [avg_k, p_start, p_end, p_iterate, kcore, z, simu, resultpath] = args

    print('zeta:', z)
    run_kcore_percolation(avg_k, p_start, p_end, p_iterate, kcore, z, simu, resultpath)
    print(f'zeta:{z} has done')


if __name__ == "__main__":

    #set the path 
    resultpath = '../kcorePercolation/figure1/result'
    networkpath = '../kcorePercolation/figure1/network'

    #set the basic information of spatial networks
    kcore = 5  #2,3,4,6
    zetas = [3,4,5,6,7,8,9,10,11,12,13,14,15,500]
    avg_k = 10 
    simu = 10 
    p_iterate = 0.001
    p_window = 0.02
    
    # Load the percolation threshold values for different zeta values
    # kcore_zeta_threshold_pc records the percolation threshold for different zeta values for L=500
    pc = kp.load(resultpath + '/kcore_zeta_threshold_pc')
    print(f'{kcore}-core pc', pc[kcore])

    args = []
    for z in zetas:
        if z == 500:
            zeta_pc = pc[kcore][200] 
        else:          
            zeta_pc = pc[kcore][z] 
            
        p_start = zeta_pc + p_window
        p_end = zeta_pc - p_window
        args.append([avg_k, p_start, p_end, p_iterate, kcore, z, simu, resultpath])
        

    pool = Pool(processes=int(len(zetas)))    
    results = pool.map(parallel_computation, args)
    pool.close()
