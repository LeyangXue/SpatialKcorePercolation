import sys
import random 
import networkx as nx 

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def kcorepercolationPlateau(G, rm_nodes, kcore, dp, p, zeta, simu,networkpath, resultpath):
    """
    Perform k-core percolation and save the largest connected component (LCC) at each step.

    Parameters
    ----------
    G : networkx.Graph
        The spatial network.
    rm_nodes : list
        List of nodes to remove.
    kcore : int
        The k-core value.
    dp : float
        Density parameter for node addition.
    p : float
        Probability of node reserved.
    zeta : int
        Zeta value for spatial network generation.
    simu : int
        Simulation number.
    networkpath : str
        Path to save the lcc nodes at each step.
    resultpath : str
        Path to save the number of lcc nodes at each step.

    Returns
    -------
    None
    """
    savepath = networkpath + f'/zeta{zeta}_dp{dp}/'
    kp.mkdirectory(savepath)

    # Remove the nodes 
    Gr = G.copy()
    N = G.order()
    
    # Count the number of lcc in the iterative process
    NOI = 0
    lccs_step_dict = {NOI:1}
    lcc = sorted(nx.connected_components(Gr), key=len, reverse=True)[0]
    lccs_step_nodes_dict = {NOI:lcc}#initial connected component
    kp.save(savepath+f'kcore{kcore}_zeta{zeta}_dp{dp}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', lcc)

    # Remove the nodes and update the degree dict
    Gr.remove_nodes_from(rm_nodes)
    degree_dict = {each[0]:each[1] for each in Gr.degree()}

    while len(degree_dict)!=0 and min(degree_dict.values()) < kcore:
        
        # Ipdate the steps and lcc of each steps
        NOI += 1
        lcc = sorted(nx.connected_components(Gr), key=len, reverse=True)[0]
        lccs_step_nodes_dict[NOI] = lcc
        lccs_step_dict[NOI] = len(lcc)/N
        kp.save(savepath+f'kcore{kcore}_zeta{zeta}_dp{dp}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', lcc)

        # Remove nodes with degree less than k
        nodes = list(Gr.nodes())
        for node in nodes:
            #print('node',node)
            if degree_dict[node] < kcore:
                Gr.remove_node(node)
        
        #update the degree dict
        degree_dict = {each[0]:each[1] for each in Gr.degree()}
    
    # Calculate the connected components
    connect_components = sorted(nx.connected_components(Gr), key=len, reverse=True)
    if len(connect_components) > 0:
        lcc = connect_components[0]
        lcc_p = len(lcc)/N 
        NOI  += 1  
        lccs_step_nodes_dict[NOI] = lcc
        lccs_step_dict[NOI] = lcc_p 
        kp.save(savepath+f'kcore{kcore}_zeta{zeta}_dp{dp}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', lcc)
    else:
        lcc = 0
        NOI  += 1  
        lccs_step_nodes_dict[NOI] = connect_components
        lccs_step_dict[NOI] = lcc
        kp.save(savepath+f'kcore{kcore}_zeta{zeta}_dp{dp}_p{round(p,7)}_simu{simu}_step{NOI}_lcc_nodes.pkl', connect_components)

    # Save the results
    savepath_lcc_num = resultpath+f'/zeta{zeta}/kcore{kcore}_zeta{zeta}_dp{dp}_p{round(p,7)}_simu{simu}_steps_lcc_num_dict.pkl'
    kp.save(savepath_lcc_num, lccs_step_dict)
    
def Runkcorepercolation(rm_nodes, avg_k, dp, kcore, zeta, simu, networkpath, resultpath):
    """
    Run k-core percolation with density parameter variation.

    Parameters
    ----------
    rm_nodes : list
        List of nodes to remove.
    avg_k : int
        Average degree.
    dp : float
        Density parameter.
    kcore : int
        The k-core value.
    zeta : int
        Zeta value for spatial network generation.
    simu : int
        Simulation number.
    networkpath : str
        Path to the network files.
    resultpath : str
        Path to save the results.

    Returns
    -------
    None
    """    
    # Load the network
    G = kp.load(networkpath+f'/NetID0_avgk{avg_k}_zeta{zeta}_spatialNet.pkl')
    
    # Add nodes to the removal list based on density parameter
    rm_nodes_copy = rm_nodes.copy()
    N = G.order()
    update_node_num  = int(dp*N)
    nodes = set(list(G.nodes()))
    rm_nodes_copy += random.sample(nodes - set(rm_nodes_copy), update_node_num)
    
    # Recalculate the critical probability
    p = 1 - len(rm_nodes_copy)/N
    kcorepercolationPlateau(G, rm_nodes_copy, kcore, dp, p, zeta, simu, networkpath, resultpath)

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

def Calculategccfailnodes(G, lccs_step_num_dict, kcore, zeta, pc, dp, simu, networkpath, resultpath):
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
    dp : float
        Change in probability for node reservation.
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
        lcc_nodes_name = networkpath+f'/zeta{zeta}_dp{dp}/kcore{kcore}_zeta{zeta}_dp{dp}_p{round(pc-dp,7)}_simu{simu_dict[zeta]}_step{step}_lcc_nodes.pkl'
        lcc_nodes = kp.load(lcc_nodes_name)
        
        # Calculate the set of removed nodes
        removed_node_set = nodeset - set(lcc_nodes)
        
        # Calculate the GCC size formed by the removed nodes
        fgcc = GetfailedGCC(G, removed_node_set)
        
        # Store the result in the dictionary
        fgcc_step[step] = fgcc

    # Save the results
    kp.save(resultpath+f'/zeta{zeta}/kcore{kcore}_zeta{zeta}_dp{dp}_p{round(pc-dp,7)}_simu{simu}_steps_fgcc_nodes_dict.pkl', fgcc_step)
    
if __name__ == "__main__":  

    # Set the result, and network paths
    resultpath = '../kcorePercolation/figure3/result'
    networkpath = '../kcorePercolation/figure3/network'
    
    #set the parameters
    avg_k = 10 #average degree
    kcore = 5 
    zetas = [7, 500] #10,100,500
    zetapc = {7:(0.718501,9), 500:(0.679851,50)} #10:(0.700501,40), 100:(0.680501,26) 
    simu_dict = {7:1, 500:2}
    zeta_updatedp = {7:[0.002, 0.005, 0.01, 0.1], 500:[0.0001, 0.0005, 0.001, 0.005]}

    #deal with each zeta
    for zeta in zetas:
        
        (pc, step) = zetapc[zeta]
        print(f'run the task with zeta={zeta}')
        
        #load the removed nodes at the criticality
        rm_nodes = kp.load(networkpath+f'/rmnodes/kcore{kcore}_zeta{zeta}_simu{simu_dict[zeta]}_p{pc}_step{step}_rmnodes.pkl')
        G = kp.load(networkpath+f'/NetID0_avgk{avg_k}_zeta{zeta}_spatialNet.pkl')
       
        #update the p as -0.002 -0.005 -0.01 -0.1
        for dp in zeta_updatedp[zeta]:
            
            print(f'run the task with zeta={zeta}, dp={round(dp,7)}')            
            Runkcorepercolation(rm_nodes, avg_k, dp, kcore, zeta, simu_dict[zeta], networkpath, resultpath)
            
            # Load the number of nodes in the largest connected component at criticality 
            lcc_num_name = resultpath+f'/zeta{zeta}/kcore{kcore}_zeta{zeta}_dp{dp}_p{round(pc-dp,7)}_simu{simu_dict[zeta]}_steps_lcc_num_dict.pkl'
            lccs_step_num_dict = kp.load(lcc_num_name)
            
            # Calculate the GCC size formed by the removed nodes at each step 
            Calculategccfailnodes(G, lccs_step_num_dict, kcore, zeta, pc, dp, simu_dict[zeta], networkpath, resultpath)
