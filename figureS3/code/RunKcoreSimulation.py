import sys
import numpy as np 
from multiprocessing import Pool 

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def setting_ps():
    """
    Define the ps (percolation probability) values for different system sizes.
    
    Returns
    -------
    dict
        Dictionary containing system sizes as keys and corresponding ps values as arrays.
    """
    ps_dict = {}
    ps_dict[100] = np.arange(0.90,0.94,0.001)
    ps_dict[300] = np.arange(0.90,0.94,0.001)
    ps_dict[500] = np.arange(0.90,0.94,0.001)
    ps_dict[1000] = np.arange(0.925,0.94,0.001)
    ps_dict[1500] = np.arange(0.925,0.94,0.001) # 0.925, 0.91-0.95
    ps_dict[2000] = np.arange(0.925,0.94,0.001) #0.91-0.95
    
    return ps_dict

def SimulationParalellComputation(args):
    """
    Run a single k-core percolation simulation.
    
    Parameters
    ----------
    args : list
        List containing simulation parameters.
        [its, G, ps, k]
        its: int
            Iteration number.
        G : NetworkX graph
            Input graph for k-core percolation.
        ps : numpy array
            Probability to reserve the nodes.
        k : int
            Target k-core number.
    
    Returns
    -------
    results : list
        Results of the k-core percolation simulation.
    """
    its, G, ps, k = args
    results = kp.Runkcorepercolation(its, G, ps, k)
    
    return results
    
def RunSimulations(G,ps,k,simulation):
    """
    Run multiple k-core percolation simulations in parallel.
   
    Parameters
    ----------
    G : NetworkX graph
       Input graph for k-core percolation.
    ps : numpy array
       Probability to reserve the nodes.
    k : int
       Target k-core number.
    simulation : int
       Number of simulations to run.
   
    Returns
    -------
    results : list
        List of results from k-core percolation simulations.
    """
    
    args = []
    for i in np.arange(simulation):
        args.append([i, G, ps, k])
        
    # Set up a pool of worker processes
    # The number of processes is set to simulation/5 for efficiency
    pool = Pool(processes=int(simulation/5))   
    # Map simulation function to input arguments and run in parallel
    results = pool.map(SimulationParalellComputation,args)
    # Close the pool of worker processes
    pool.close()
    
    return results     

if __name__ == "__main__":

    # Set the path
    networkpath = '../kcorePercolation/figureS3/network'
    resultpath = '../kcorePercolation/figureS3/result'
    
    # Set the parameters
    Ls = [300, 500, 1000, 1500, 2000]
    ps = setting_ps() # Probability to reserve the nodes
    simulations = {300:200, 500:200, 1000:100, 1500:25, 2000:25} # Simulation times
    kcore = 5
    
    for L  in Ls:
        
        print(f'run the percolation with L={L}')
        
        # Load the network
        networkname = f'Net_L{L}_avgk7.5_zeta10_spatialNet.pkl'
        G = kp.load(networkpath+'/'+networkname)
        
        # Run the percolation
        print('Begin to run the percolation...')
        results = RunSimulations(G, ps, kcore, simulations[L])
        resultkname = f'Results_L{L}_avgk7.5_zeta10'
        kp.save(resultpath+'/'+ resultkname + '_kcore'+str(kcore)+'_Percolation.pkl', results)
