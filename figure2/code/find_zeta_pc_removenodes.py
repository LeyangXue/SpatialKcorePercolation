import sys
import random 
import networkx as nx

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def kcorepercolationPlateau(G, rm_nodes, k):
    '''
    Perform k-core percolation and track the largest connected component (LCC) over steps.

    Parameters
    ----------
    G : networkx.Graph
        The spatial network.
    rm_nodes : list
        List of nodes to remove.
    k : int
        The k-core value.

    Returns
    -------
    lcc : float
        Fraction of nodes in the largest connected component.
    NOI : int
        Number of iterations until k-core percolation stabilizes.

    '''
    
    Gr = G.copy()
    N = G.order()
    
    # Initialize variables
    NOI = 0
    lccs_step_dict = {NOI:1}

    # Remove the nodes and update the degree dictionary
    Gr.remove_nodes_from(rm_nodes)
    degree_dict = {each[0]:each[1] for each in Gr.degree()}
    
    # Perform k-core percolation until stabilization
    while len(degree_dict)!=0 and min(degree_dict.values()) < k:
        
        # Update the step count and LCC fraction
        NOI += 1
        lccs_num = len(sorted(nx.connected_components(Gr), key=len, reverse=True)[0])
        lccs_step_dict[NOI] = lccs_num/N
        
        # Remove nodes with degree less than k
        nodes = list(Gr.nodes())
        for node in nodes:
            #print('node',node)
            if degree_dict[node] < k:
                Gr.remove_node(node)
        
        # Update the degree dictionary
        degree_dict = {each[0]:each[1] for each in Gr.degree()}
    
    # Calculate the connected components after stabilization
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
    
    return lcc, NOI

def Runkcorepercolation(avg_k, p_start, p_end, p_iterate, kcore, z, simu, resultpath, networkpath, gccpath):
    '''
    Run k-core percolation on a spatial network with varying probabilities of node removal.

    Parameters
    ----------
    avg_k : int
        Average degree of the spatial network.
    p_start : float
        Starting probability of node removal.
    p_end : float
        Ending probability of node removal.
    p_iterate : float
        Probability increment for each iteration.
    kcore : int
        The k-core value.
    z : int
        Zeta value for spatial network generation.
    simu : int
        Simulation number.
    resultpath : str
        Path to save result files.
    networkpath : str
        Path to load spatial network files.
    gccpath : str
        Path to save largest connected component files.

    Returns
    -------
    None.

    '''
    
    # Load the spatial network
    G = kp.load(networkpath+f'/NetID0_avgk{avg_k}_zeta{z}_spatialNet.pkl')

    # Obtain the initial removed nodes
    N = G.order()
    initial_rm_nodes_num = int((1-p_start)*N)
    end_rm_nodes_num = int((1-p_end)*N)
    update_node_num  = int(p_iterate*N)

    # Initialize removed nodes set and step count
    nodes = set(list(G.nodes()))
    rm_nodes = random.sample(nodes, initial_rm_nodes_num)
    rm_nodes_set = set()
    step = 0
    
    # Main loop for k-core percolation with varying probabilities
    while len(rm_nodes) <= end_rm_nodes_num: 

        # Calculate probability to reserve nodes
        p = 1 - len(rm_nodes)/N
        # Perform k-core percolation
        [lcc, NOI] =  kcorepercolationPlateau(G, rm_nodes, kcore, z, simu, step, resultpath)
        
        # Handle different conditions based on lcc and current number of removed nodes
        if lcc == 0 and (len(rm_nodes) == initial_rm_nodes_num):
            # Reset removed nodes to initial state if lcc is zero and at the beginning
            rm_nodes = random.sample(nodes, initial_rm_nodes_num)
            
        elif lcc != 0 and (len(rm_nodes) == end_rm_nodes_num):
            # Reset removed nodes to initial state if lcc is non-zero and at the end
            rm_nodes = random.sample(nodes, initial_rm_nodes_num) 
            
        elif (len(rm_nodes) != initial_rm_nodes_num) and lcc == 0:
            # Save removed nodes if lcc is zero and not at the beginning or end
            savegccpath = gccpath+f'/kcore{kcore}/rmnodes'
            kp.mkdirectory(savegccpath)
            kp.save(savegccpath+f'kcore{kcore}_zeta{z}_simu{simu}_p{p}_step{step}_rmnodes.pkl', rm_nodes)
            break
        else:   
            # Update removed nodes by adding new nodes to remove
            rm_nodes_set.update(rm_nodes)
            rest_nodes_set = nodes - rm_nodes_set
            update_nodes = random.sample(rest_nodes_set, update_node_num)
            rm_nodes += update_nodes
            step += 1

        # Print progress information
        print(f'run the kcore percolation for kcore={kcore}, zeta={z}, step={step}, p={round(p,7)}, simu={simu}, lcc={lcc}, NOI={NOI}')

if __name__ == "__main__":  

    # Set the paths for files of result, network, and removed node at criticality (gcc) 
    resultpath = '../kcorePercolation/figure2/result'
    networkpath = '../kcorePercolation/figure2/network'
    gccpath = '../kcorePercolation/figure2/gcc'

    # Set basic information about spatial networks
    kcore = 5
    zetas = [4, 6, 7, 8, 10, 12, 14, 16, 18, 20, 50, 100, 500] #50,60,70,80,100,500
    avg_k = 10 # Average degree
    simu = 1
    
    # Set the window of probability to reserve the nodes and the probability of updation
    p_iterate = 0.0001 # Incremental probability
    p_upper_window = 0.005 # Upper window for probability range
    p_lower_window = 0.01  # Lower window for probability range

    #Load the pre-computed percolation threshold for k-core for L = 500
    pc = kp.load(resultpath+'/kcore_zeta_threshold_pc')
    print(f'{kcore}-core pc', pc[kcore])
    
    # Iterate over different zeta values
    for z in zetas:
        # Set the probability to reserve the nodes
        if z == 500:
            zeta_pc = pc[kcore][200] 
        elif z in [60,70,80,90]:
            zeta_pc = pc[kcore][50]
        else:          
            zeta_pc = pc[kcore][z] 
            
        # Calculate start and end probabilities for percolation
        p_start = zeta_pc + p_upper_window
        p_end = zeta_pc - p_lower_window
        
        # Perform k-core percolation
        Runkcorepercolation(avg_k, p_start, p_end, p_iterate, kcore, z, simu, resultpath, networkpath, gccpath)


