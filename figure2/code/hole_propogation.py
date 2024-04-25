# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 19:43:01 2024

@author: Leyang Xue
"""

import sys
# Add package path
packagepath = 'F:/work/work12/kcorePercolation'  # Manually specify the path to your script
sys.path.append(packagepath)
from utils import kcorePercolation as kp # Import custom utility functions

def calculateRh(L, networkpath, gccpath, kcore, zeta, pc, simu):
    '''
    Calculate the radius of holes at different times during percolation.

    Parameters
    ----------
    L : int
        Size of the lattice.
    networkpath : str
        Path to the directory containing network data.
    gccpath : str
        Path to the directory containing giant connected component data.
    kcore : int
        The k-core parameter.
    zeta : int
        The zeta parameter.
    pc : float
        The critical probability.
    simu : int
        The simulation number.

    Returns
    -------
    None
        This function saves the radius of holes at different times as a dictionary.
    '''
    
    gcc_nodes_path = gccpath+f'/kcore{kcore}/zeta{zeta}_lcc_nodes'
    lccs_step_dict = kp.load(gcc_nodes_path+f'/kcore{kcore}_zeta{zeta}_p{pc}_simu{simu}_steps_lcc_num_dict.pkl')
    
    # Load the network
    G = kp.load(networkpath+f'/NetID0_avgk10_zeta{zeta}_spatialNet.pkl')
    
    rh = {}
    for t in lccs_step_dict.keys():
        if t >= 1: 
            lcc_times_t = kp.load(gcc_nodes_path+f'/kcore{kcore}_zeta{zeta}_p{pc}_simu{simu}_step{t}_lcc_nodes.pkl')
            matrix_t = kp.transform2M(lcc_times_t, G)
            [origin, diameter]= kp.LocalizeCircleOrigin(matrix_t, L)  # Calculate the position of origin for the circle
            rh[t] = diameter / 2
            print(f't={t}, matrix_t={rh[t]}')
            
    kp.save(networkpath+f'/hole_propogation/kcore{kcore}_zeta{zeta}_simu{simu}_rh_time_dict.pkl', rh)

if __name__ == '__main__':
    
    # Setting the paths
    networkpath = '../kcorePercolation/figure2/network'
    gccpath = '../kcorePercolation/figure2/gcc'
    
    # Setting parameters
    kcore = 5
    simu = 1
    zetas = [10] 

    # Please download the network and zeta{zeta}_lcc_nodes data for zeta = 10 on the Mendeley dataset if running the code directly.
    # Also, ensure to run the find_zeta_pc_removenodes.py script in advance if running other cases for zeta = 4, 7.... and  zeta = 6(pc=0.721001), zeta = 8(pc=0.710351), zeta = 12(pc=0.69765), zeta = 14(pc=0.692701), zeta = 16(pc=0.691051), zeta = 18(pc= 0.687351), zeta = 20(pc = 0.6876), zeta = 50(pc = 0.681601), zeta = 60(pc = 0.679701), zeta = 100(pc = 0.680501), and zeta = 500(pc = 0.679151).

    pcs_zeta = {10:0.700501}
    L = 1000
    
    # Iterating over different zeta values
    for zeta in zetas:
        
        # Printing current task    
        print(f'run the task with zeta={zeta}')
        
        # Getting the critical probability for the current zeta
        pc = pcs_zeta[zeta]
        
        # Calculating radius of holes over time for the current zeta
        calculateRh(L, networkpath, gccpath, kcore, zeta, pc, simu)
