import sys
import numpy as np

# Add package path
packagepath = '../kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def generate_networks(networkpath, L=1000, ks=[5, 10], zetas=[4, 5, 100]):
    """
    Generate networks for k-core percolation.

    Parameters
    ----------
    networkpath : str
        Path to save the generated networks.
    L : int, optional
        Size of the lattice. Default is 1000.
    ks : list of int, optional
        List of average degrees to generate networks for. Default is [5, 10].
    zetas : list of int, optional
        List of zeta values to generate networks for. Default is [4, 5, 100].

    Returns
    -------
    None
    """
    
    # Define the range of coordinates for the lattice
    range_coordinate = np.arange(0, L, 1)

    # Generate networks for different average degrees and zeta values
    for avg_k in ks:
        print('Generating networks for average degree %d' % avg_k)
        
        for zeta in zetas:
            # Set the filename
            filename = 'network_L_%d_avg_k_%d_zeta_%d' % (L, avg_k, zeta) 
            M = int((L * L * avg_k) / 2)  # Total number of edges
            
            # Generate links
            links = kp.expvariate_sequence(zeta, M, L)

            # Save the links
            linkname = networkpath + '/' + filename + '_links.pkl'
            kp.save(linkname, links)

            # Create and save the network graph
            G = kp.NetworkCreate(range_coordinate, links, networkpath, filename)

if __name__ == "__main__":
    
    #set the parameters
    L = 1000
    N = L*L
    ks = [5, 10]
    zetas = [4, 5, 10, 100]
    range_coordinate = np.arange(0,L,1) #range of coordinate 

    #set the path 
    networkpath = '../kcorePercolation/figureS2/network'
    # Generate networks
    generate_networks(networkpath)
            