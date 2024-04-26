import sys
import networkx as nx

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def GetfailedGCC(G, removed_node_set):
    '''
    Calculate the size of the largest connected component (GCC) 
    formed by the removed nodes.

    Parameters
    ----------
    G : networkx.Graph
        The original spatial network.
    removed_node_set : set
        Set of nodes that have been removed.

    Returns
    -------
    fgcc : float
        Fraction of nodes in the largest connected component formed by the removed nodes.
    '''
    
    G_copy = G.copy()
    H = G_copy.subgraph(removed_node_set)
    fgcc_size = len(max(nx.connected_components(H), key=len))
    fgcc = fgcc_size/G_copy.order()
    
    return fgcc

def Calculategccfailnodes(G, lccs_step_num_dict, kcore, zeta, pc, simu, networkpath, resultpath):
    '''
    Calculate the size of the largest connected component (GCC) 
    formed by the removed nodes at each step of k-core percolation.

    Parameters
    ----------
    G : networkx.Graph
        The original spatial network.
    lccs_step_num_dict : dict
        Dictionary containing the number of nodes in the largest connected component at each step.
    kcore : int
        The k-core value.
    zeta : int
        Zeta value for spatial network generation.
    pc : float
        Critical probability of node reservation.
    simu : int
        Simulation number.
    networkpath : str
        Path to the network files.
    resultpath : str
        Path to save the results.

    Returns
    -------
    None
    '''
   
    # Get the set of all nodes in the network
    nodeset = set(G.nodes())
    
    # Dictionary to store the GCC size at each step
    fgcc_step = {}
    
    # Iterate over each step
    for step in lccs_step_num_dict.keys():
        
        # Load the nodes in the largest connected component at the current step
        lcc_nodes_name = networkpath+f'/zeta{zeta}_lcc_nodes/kcore{kcore}_zeta{zeta}_p{round(pc,7)}_simu{simu_dict[zeta]}_step{step}_lcc_nodes.pkl'
        lcc_nodes = kp.load(lcc_nodes_name)
        
        # Calculate the set of removed nodes
        removed_node_set = nodeset - set(lcc_nodes)
        
        # Calculate the GCC size formed by the removed nodes
        fgcc = GetfailedGCC(G, removed_node_set)
        
        # Store the result in the dictionary
        fgcc_step[step] = fgcc

    # Save the results
    kp.save(resultpath+f'/zeta{zeta}/kcore{kcore}_zeta{zeta}_p{round(pc,7)}_simu{simu}_steps_fgcc_nodes_dict.pkl', fgcc_step)

if __name__ == "__main__":
    
    # Set the result, and network paths
    resultpath = '../kcorePercolation/figureS7/result'
    networkpath = '../kcorePercolation/figureS7/network'
    
    # Set the parameters 
    avg_k = 10 
    kcore = 5 
    zetas = [10, 100]  # Zeta values 10, 100
    zetapc = {10:0.700501,100:0.679701}
    simu_dict = {10:1, 100:2}

    # Deal with each zeta
    for zeta in zetas:
        
        pc = zetapc[zeta]
        print(f'run the task with zeta={zeta}')

        # Load the network 
        G = kp.load(networkpath+f'/NetID0_avgk{avg_k}_zeta{zeta}_spatialNet.pkl')
        
        # Load the number of nodes in the largest connected component at criticality 
        lcc_num_name = resultpath+f'/zeta{zeta}/kcore{kcore}_zeta{zeta}_p{round(pc,7)}_simu{simu_dict[zeta]}_steps_lcc_num_dict.pkl'
        lccs_step_num_dict = kp.load(lcc_num_name)

        # Calculate the GCC size formed by the removed nodes at each step
        Calculategccfailnodes(G, lccs_step_num_dict, kcore, zeta, pc, simu_dict[zeta], networkpath, resultpath)
        