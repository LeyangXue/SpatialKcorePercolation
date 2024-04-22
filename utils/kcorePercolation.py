#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 02:44:59 2016

@author: leyangx
"""

import numpy as np
import os
import pickle as pk 
import random 
import networkx as nx
import itertools as itert
import matplotlib.pyplot as plt

def load(filename):
   '''
    This function can load some data saved as with a specific format by dump

    Parameters
    ----------
    filename : str 
       path that need to load the file, and the file is saved in dump function.

    Returns
    -------
    TYPE file type, is the same as the previous type 
       file.

    ''' 
   with open(filename, 'rb') as input_file:  
        try:
            return pk.load(input_file)
        
        except EOFError:
            return None

def save(ResultPath, results):
   '''
   Save the results to a file using pickle serialization.

    Parameters:
       ResultPath : str 
           The path to the file where results will be saved.
       results : object
           The results to be saved.
      
    Returns
    -------
    None.

    '''
   pk.dump(results, open(ResultPath, "wb"))


def largestComponent(G):
    '''
    Parameters
    ----------
    dataset : list
       network edgelist.

    Returns
    -------
    G : graph
        simple network.
    '''
    largest_cc = max(nx.connected_components(G), key=len)
    largest_cc_G = G.subgraph(largest_cc)
    G=nx.Graph(largest_cc_G)
    G.remove_edges_from(nx.selfloop_edges(G))
    
    return G

def expvariate_sequence(kesi, n, L):
    '''
    Generate a sequence of link lengths using exponential distribution.
    
    Parameters
    ----------
    kesi : int
        Zeta, denotes the characteristic length of a link.
    n : int
        The number of links to generate.
    L : int
        The system length.

    Returns
    -------
    sample : list
        The sequence of link lengths.
    '''
    
    sample = [] 
    i = 0    
    while i < n: 
        p = random.random()
        x = -kesi*np.log(1-p)
        length = x
        if length <= L/2:
            sample.append(length)
            i += 1
            
    return sample

def kcorepercolation(G,p,k):
    '''
    Perform k-core percolation on a given network.
    
    Parameters
    ----------
    G : graph
        The input network.
    p : float
        The proportion of initial remaining nodes.
    k : int 
        The minimum degree of nodes in the core.

    Returns
    -------
    lcc : int
        The size of the giant connected component.
    NOI : int 
        The number of iteration steps.
    sec_lccs : list
        The size of the second largest cluster at each step.
    '''
    
    #obtain the removed nodes 
    rm_nodes_num = int((1-p)*G.order())
    nodes = list(G.nodes())
    rm_nodes = random.sample(nodes, rm_nodes_num)
    
    #remove the nodes 
    Gr = G.copy()
    Gr.remove_nodes_from(rm_nodes)
    degree_dict = {each[0]:each[1] for each in Gr.degree()}

    #count the number of stages, second_lcc in the iterative process
    NOI = 0
    sec_lccs = []

    while len(degree_dict)!=0 and min(degree_dict.values()) < k:
        
        #update the steps and second lcc
        NOI += 1
        if len(sorted(nx.connected_components(Gr), key=len, reverse=True)) > 1:
            seclcc_num = len(sorted(nx.connected_components(Gr), key=len, reverse=True)[1])
        else:
            seclcc_num = 0
        sec_lccs.append(seclcc_num) 
        
        #remove nodes with degree less than k
        nodes = list(Gr.nodes())
        for node in nodes:
            #print('node',node)
            if degree_dict[node] < k:
                Gr.remove_node(node)
        
        #update the degree dict
        degree_dict = {each[0]:each[1] for each in Gr.degree()}
    
    #calculate the connected components
    connect_components = sorted(nx.connected_components(Gr), key=len, reverse=True)
    if len(connect_components) > 0:
        lcc = len(connect_components[0])
    else:
        lcc = 0
        
    return lcc, NOI, sec_lccs

def kcorepercolationPlateau(G,p,k):
    '''
    Perform k-core percolation and return the giant component size at each step.

    Parameters
    ----------
    G : graph 
        The input network.
    p : float
        The proportion of initial remaining nodes.
    k : int 
        The minimum degree of nodes in the core.
        
    Returns
    -------
    lcc : int
        The number of nodes in the giant connected component.
    NOI : int 
        The number of iteration steps.
    lccs_step : list
        The size of the giant connected component at each step.
    '''
    
    #obtain the removed nodes 
    rm_nodes_num = int((1-p)*G.order())
    nodes = list(G.nodes())
    rm_nodes = random.sample(nodes, rm_nodes_num)
    
    #remove the nodes 
    Gr = G.copy()
    Gr.remove_nodes_from(rm_nodes)
    degree_dict = {each[0]:each[1] for each in Gr.degree()}

    #count the number of stages, second_lcc in the iterative process
    NOI = 0
    lccs_step = []

    while len(degree_dict)!=0 and min(degree_dict.values()) < k:
        
        #update the steps and lcc of each steps
        NOI += 1
        lccs_num = len(sorted(nx.connected_components(Gr), key=len, reverse=True)[0])
        lccs_step.append(lccs_num) 
        
        #remove nodes with degree less than k
        nodes = list(Gr.nodes())
        for node in nodes:
            #print('node',node)
            if degree_dict[node] < k:
                Gr.remove_node(node)
        
        #update the degree dict
        degree_dict = {each[0]:each[1] for each in Gr.degree()}
    
    #calculate the connected components
    connect_components = sorted(nx.connected_components(Gr), key=len, reverse=True)
    if len(connect_components) > 0:
        lcc = len(connect_components[0])
    else:
        lcc = 0
    
    NOI  += 1  
    lccs_step.append(lcc)
    
    return lcc, NOI, lccs_step

def kcorepercolationLccNodes(G,p,k):
    '''
    K-core percolation and return all nodes from the giant component at each step.
    
    Parameters
    ----------
    G : graph 
        The input network.
    p : float
        The proportion of initial remaining nodes.
    k : int 
        The minimum degree of nodes in the core.

    Returns
    -------
    NOI : int
        The number of iteration steps.
    lcc_step_nodes_dict : dict
        A dictionary containing the nodes of the giant component at each step.

    '''
    #obtain the removed nodes 
    rm_nodes_num = int((1-p)*G.order())
    nodes = list(G.nodes())
    rm_nodes = random.sample(nodes, rm_nodes_num)
    
    #remove the nodes 
    Gr = G.copy()
    Gr.remove_nodes_from(rm_nodes)
    degree_dict = {each[0]:each[1] for each in Gr.degree()}

    #count the number of stages, second_lcc in the iterative process
    NOI = 0
    lcc_step_nodes_dict = {}
    lcc_step_nodes_dict[NOI] = sorted(nx.connected_components(Gr), key=len, reverse=True)#initial connected component

    while len(degree_dict)!=0 and min(degree_dict.values()) < k:
        
        #remove nodes with degree less than k
        nodes = list(Gr.nodes())
        for node in nodes:
            #print('node',node)
            if degree_dict[node] < k:
                Gr.remove_node(node)
        
        #update the degree dict
        degree_dict = {each[0]:each[1] for each in Gr.degree()}
    
        #update the steps and lcc nodes of each steps
        NOI += 1
        lcc_step_nodes_dict[NOI] = sorted(nx.connected_components(Gr), key=len, reverse=True)
    
    return NOI, lcc_step_nodes_dict

def Runkcorepercolation(its,G,ps,kcore):
    '''
    Run k-core percolation simulation for multiple values of p.

    Parameters
    ----------
    its : int
        Number of iterations/simulations to run for each p.
    G : graph 
        The input network.
    ps : list
        List of p values (proportion of remaining nodes).
    kcore : int 
        The minimum degree of nodes in the core.

    Returns
    -------
    kcore : int
        The k-core value used in the simulation.
    lcc_list : list
        List containing the number of nodes in the giant connected component for each p.
    NOI_list : list
        List containing the number of iteration steps for each p.
    lccstep_list : list
        List containing the size of the giant connected component at each step for each p.
        
    '''
    print('run the kcore percolation for k=%d'%kcore)

    lcc_list = []
    NOI_list = []
    lccstep_list = []
    for p in ps:
        print('----run current %d simulation times---p:%f'%(its,p))
        [lcc,noi,lccstep] = kcorepercolationPlateau(G,p,kcore)
        lcc_list.append(lcc)
        NOI_list.append(noi)
        lccstep_list.append(lccstep)
    
    return kcore, lcc_list, NOI_list, lccstep_list

def PlotAxes(ax,xlabel,ylabel, title, mode=False):
    '''
    Decorate the axes
    
    Parameters
    ----------
    ax : axes
        The axes object to be decorated.
    xlabel : str
        The label for the x-axis.
    ylabel : str
        The label for the y-axis.
    title : str
        The title of the plot.
    mode : bool, optional
        Whether to show the legend. The default is False.
    
    Returns
    -------
    None.
    
    '''
    fontsize = 14
    font_label = {'family': "Arial", 'size':fontsize}
    
    n_legend = 12
    ax.set_xlabel(xlabel,  fontdict = font_label)
    ax.set_ylabel(ylabel, fontdict = font_label)
    ax.set_title(title, loc='left',fontdict = {'family': "Arial", 'size':fontsize})
    ax.tick_params(direction='out', which='both',length =4, width=1, pad=1,labelsize=n_legend)
    
    #ax.minorticks_on()
    if mode == True:
        ax.legend(loc='best', framealpha=0, fontsize=n_legend)

def PlotSAxes(ax,xlabel,ylabel, title, mode=False):
    '''
    Decorate the axes with smaller font sizes.
    
    Parameters
    ----------
    ax : axes
        The axes object to be decorated.
    xlabel : str
        The label for the x-axis.
    ylabel : str
        The label for the y-axis.
    title : str
        The title of the plot.
    mode : bool, optional
        Whether to show the legend. The default is False.
    
    Returns
    -------
    None.
    
    '''
    fontsize = 10
    font_label = {'family': "Arial", 'size':fontsize}
    
    n_legend = 8
    ax.set_xlabel(xlabel,  fontdict = font_label)
    ax.set_ylabel(ylabel, fontdict = font_label)
    ax.set_title(title, loc='left',fontdict = {'family': "Arial", 'size':fontsize})
    ax.tick_params(direction='out', which='both',length =2, width=0.5, pad=0.5,labelsize=n_legend)
    
    #ax.minorticks_on()
    if mode == True:
        ax.legend(loc='best', framealpha=0, fontsize=n_legend)
        
def findCandidateCoordinate(link_length):
    '''
    Find candidate coordinates for placing nodes in a lattice around a circle.

    Parameters
    ----------
    link_length : float
        The length of the link.

    Returns
    -------
    all_ccoordinates : list
        List of candidate coordinates.
    cc_circle : dict
        Dictionary containing distances from lattice points to the circle.
    cc_lattice : dict
        Dictionary containing distances from lattice points to the other lattice points.
        
    '''
    
    # Determine maximum lattice point distance
    max_l = np.ceil(link_length)
    # Enumerate lattice points
    enumerate_length = np.arange(0, max_l, 1)
    
    # Store distances from lattice points to the circle
    cc_lattice = {} #sotre the distances
    cc_circle = {} #sotre the distances
    
    # Iterate over lattice points
    for deltax in enumerate_length:
        
        # Calculate deltay for the given deltax
        deltay = np.sqrt(np.power(link_length,2)-np.power(deltax, 2))
        
        #find the coordinate (deltax, deltay) on the circle
        lattice_deltax = deltax
        lattice_deltay = round(deltay,0)
        
        # Calculate shortest distance from lattice point to circle
        cc_circle[(deltax,deltay)] = abs(np.sqrt(lattice_deltax**2+lattice_deltay**2)-link_length)
        cc_lattice[(lattice_deltax,lattice_deltay)] = abs((lattice_deltay+deltay)*(lattice_deltay-deltay)) #calculate the shortest distance from the point (lattice_deltax, lattice_deltay) to circle
        
    #Find lattice point closest to the circle
    s_value_k = {cc_lattice[each]:each for each in cc_lattice.keys()} 
    min_s = min(cc_lattice.values())
    coord_sd = s_value_k[min_s] #coordinate that close to circle with shortset distance     
    
    # Generate all candidate coordinates 
    x = coord_sd[0]
    y = coord_sd[1]
    all_ccoordinates = [(x,y),(-x,y),(x,-y),(-x,-y),(y,x),(-y,x),(y,-x),(-y,-x)]
    
    return all_ccoordinates, cc_circle, cc_lattice
    
def CheckPeriodicBoundary(x, y, max_coordinate):
    '''
    Check if a point (x, y) lies within the periodic boundary conditions.

    Parameters
    ----------
    x : int
        x-coordinate of the point.
    y : int
        y-coordinate of the point.
    max_coordinate : int
        Maximum value for the coordinates.

    Returns
    -------
    x_upd : int
        Updated x-coordinate after considering periodic boundary conditions.
    y_upd : int
        Updated y-coordinate after considering periodic boundary conditions.
    '''
    x_upd = x
    y_upd = y
    L = max_coordinate+1
    
    if x > max_coordinate:
        x_upd = x%L
    if y > max_coordinate:
        y_upd = y%L
    if x < 0:
        x_upd = (x+L)%L
    if y < 0:
        y_upd = (y+L)%L
    
    return x_upd, y_upd

def CheckNeighbors(G, source_node, nodex, nodey, cc_shortest, max_coordinate):
    '''
    Check neighboring nodes of a source node and update coordinates if necessary.

    Parameters
    ----------
    G : graph
        The input network.
    source_node : tuple
        Tuple representing the source node coordinates.
    nodex : int
        Current x-coordinate of the node.
    nodey : int
        Current y-coordinate of the node.
    cc_shortest : list
        List of candidate coordinates.
    max_coordinate : int
        Maximum value for the coordinates.

    Returns
    -------
    nodex : int
        Updated x-coordinate of the node.
    nodey : int
        Updated y-coordinate of the node.
    state : bool
        Boolean indicating if there are effective candidate nodes.

    '''
    count = 0
    state = True #judge whether there is a effective candidate 
    
    # Check if the current node is a neighbor of the source node
    while (nodex,nodey) in list(G.neighbors(source_node)):#
        # Reassign select the candidate
        candidate_coord =  random.choice(cc_shortest) #random selects a candidate
        #get a new assigned node
        assigned_node_x = source_node[0]+candidate_coord[0]
        assigned_node_y = source_node[1]+candidate_coord[1]
        # Judge or update the coordinate 
        [nodex, nodey]= CheckPeriodicBoundary(assigned_node_x, assigned_node_y,max_coordinate)
        count += 1
        # If count exceeds the number of candidate coordinates, set state to False and break loop
        if count > len(cc_shortest):
            state = False
            break
        
    return nodex, nodey, state

def NetworkCreate(range_coordinate, links, networkpath, filename): 
    '''
    Create a spatial network with given range of coordinates and links.

    Parameters
    ----------
    range_coordinate : list
        Range of coordinates for the network.
    links : list
        List of link lengths.
    networkpath : str
        Path to save the network file.
    filename : str
        Name of the network file.

    Returns
    -------
    G : graph
        The generated spatial network.

    '''
    # 1. Create the empty graph    
    G = nx.Graph()
    max_coordinate = max(range_coordinate) #calculate the maximum coordinates

    # 2. Generate the nodes with coordinates and add the nodes
    for coordinate in itert.product(range_coordinate,range_coordinate):
        G.add_node(coordinate)
        
    # 3. Add the edges
    # Assign the links to nodes
    state = True
    nodes = list(G.nodes())
    link_num = 0
    while link_num < len(links):
        
        link_length = links[link_num] #obatin the length of links
        
        #1.randomly select a source node and find a candiate coordinate base
        source_node = random.choice(nodes)
        [cc_shortest,cc_circle,cc_lattice]  = findCandidateCoordinate(link_length)#find the candiate coordinates whose distance is closest to l
        cc_shortest_set= list(set(cc_shortest))
        
        # If link_length < 0.5, change (0,0) to [(0,1), (0,-1), (1,0), (-1,0)]
        if cc_shortest_set[0] == (0,0): #change (0,0) as [(0,1),(0,-1),(1,0),(-1,0)] 
            cc_shortest_set = [(0,1),(0,-1),(1,0),(-1,0)]  
        
        # Randomly select a candidate node   
        candidate_coord = random.choice(cc_shortest_set) 
        
        # 2. Get a new assigned node
        assigned_node_x = source_node[0]+candidate_coord[0]
        assigned_node_y = source_node[1]+candidate_coord[1]
            
        # 3. Update the coordinate  
        [nodex, nodey]= CheckPeriodicBoundary(assigned_node_x, assigned_node_y, max_coordinate)

        # 4. Judge whether new node is in the neighbor of source node to avoid multiple edges between two nodes
        if (nodex,nodey) in list(G.neighbors(source_node)):
            [nodex, nodey, state] = CheckNeighbors(G, source_node, nodex, nodey, cc_shortest_set, max_coordinate)
        
        # 5. Successfully add the link into the network if state == True
        if state  == True:
            G.add_edge(source_node,(nodex,nodey))
            link_num += 1
            print('%d links have been added'%link_num)
        else: #change the source node 
            continue 
        
    # Examine whether the network has been created successfully
    if link_num == len(links):
        print('network has been created successfully')
        save(networkpath+'/'+filename + '_spatialNet.pkl', G)
    else:
        print('failed to create')
    
    return G        

def Fn(x, sample, m): 
    '''
    Estimate the empirical distribution function of a sample.
    
    Parameters
    ----------
    x : array_like
        The points at which the distribution function will be estimated.
    sample : array_like
        The sample data.
    m : int
        The number of bins for grouping the sample data.
    
    Returns
    -------
    b : array_like
        The bin edges.
    f : array_like
        The empirical cumulative distribution function values.
    p : array_like
        The empirical probability mass function values.
    y : array_like
        The estimated distribution function values at points `x`.
    
    '''                  
    n = len(sample) #sample size 
    f, b=np.histogram(sample, bins=m) #calculate the frequence number and group 
    f = f/n #calculate the ratio of frequence 
    p = f.copy()
    for k in range(m-1): #calculate the accumulate ratio 
        f[k+1]+=f[k]
        
    y=np.zeros(len(x))#store the x distribution according the experience distribution  
    for i in range(1, m): #calculate the probability within each group
        d=np.where((x>b[i-1])&(x<=b[i]))
        y[d]=f[i-1]
    d=np.where(x>b[m-1]) #calculate the probability at the last group 
    y[d]=1
    
    return b, f, p, y  #return the probability distribution of sample i.e.(p), accumulate distribution of sample and sample within x's interval(f,y) 

def plot_exponent_distribution(sample, kesi, figurepath, filename):
    '''
    Plot the exponent distribution.
    
    Parameters
    ----------
    sample : array_like
        The sample data.
    kesi : float
        The characteristic length.
    figurepath : str
        The path to save the figure.
    filename : str
        The name of the figure file.
    
    Returns
    -------
    None.
    
    '''
    m = int(max(sample))
    x = np.linspace(int(min(sample)), int(max(sample)),m)
    interval_exp, accumprob_exp, prob_exp, yprob_x = Fn(x, sample, m)
    fig, ax = plt.subplots(1,2,figsize=(8,4), tight_layout = True)
    
    ax[0].plot(interval_exp[1:],prob_exp,'o', ms=6, mfc='white')
    ax[0].axvline(kesi, ls='--')
    ax[0].plot(x, np.exp(-x/kesi)/sum(np.exp(-x/kesi)), color='black', ls='-')
    ax[0].set_xscale('log')
    ax[0].set_yscale('log')
    PlotAxes(ax[0], '$l$', r'$p(l)$', '(a)')
    
    ax[1].plot(x, yprob_x,'-')
    ax[1].axvline(kesi, ls='--')
    x_index = np.where(x>kesi)[0][0]
    cum = np.cumsum((np.exp(-x/kesi)/sum(np.exp(-x/kesi))))
    ax[1].axhline(cum[x_index],ls='--')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    PlotAxes(ax[1], '$l$', 'accum. ' + r'$p(l)$', '(b)')
    
    plt.savefig(figurepath + '/' +filename + '_linklength.png', dpi =300)
  
def mkdirectory(path):
    '''
    Create a directory if it doesn't exist.
     
    Parameters
    ----------
    path : str
        The path of the directory to be created.
     
    Returns
    -------
    None.
     
    '''
    
    if not os.path.exists(path):
        os.makedirs(path)

def IdentifyPcNumerical(lcc_p, threshold, ps):
    '''
    Identify percolation threshold numerically.

    Parameters
    ----------
    lcc_p : array_like
        List of giant connected component sizes.
    threshold : float
        Threshold value for identifying percolation.
    ps : array_like
        List of proportions of remaining nodes.

    Returns
    -------
    pc : float
        Percolation threshold.

    '''
    # Find the index of the first point where the LCC size exceeds the threshold
    first_point = np.where(np.array(lcc_p) > threshold)[0][0]
    
    # Obtain the percolation threshold corresponding to the identified index
    pc = ps[first_point]
    
    return pc 


def IdentifyPcNOI(step, ps):
    '''
    Identify percolation threshold based on the number of iteration steps.

    Parameters
    ----------
    step : array_like
        List of iteration steps.
    ps : array_like
        List of proportions of remaining nodes.

    Returns
    -------
    pc : float
        Percolation threshold.

    '''
    # Find the index of the maximum number of iteration steps
    first_point = np.argmax(step)
    
    # Obtain the percolation threshold corresponding to the identified index
    pc = ps[int(first_point + 1)]  # Add 1 to the index to account for zero-based indexing
    
    return pc

def IdentifyPcSeclccSum(seclcc, ps):
    '''
    Identify percolation threshold based on the sum of sizes of second largest clusters.

    Parameters
    ----------
    seclcc : list
        List of sizes of second largest clusters at each step.
    ps : array_like
        List of proportions of remaining nodes.

    Returns
    -------
    pc : float
        Percolation threshold.

    '''
    sec_sum = []  # Initialize list to store sums of second largest cluster sizes
    for each in seclcc:  # Iterate over sizes of second largest clusters
        if len(each) != []:  # Check if the list of sizes is not empty
            sec_sum.append(sum(each))  # Calculate and append the sum of sizes
        else:
            sec_sum.append(0)  # Append 0 if the list of sizes is empty
    
    # Find the index of the maximum sum of sizes
    first_point = np.argmax(np.array(sec_sum))
    
    # Obtain the percolation threshold corresponding to the identified index
    pc = ps[first_point]
    
    return pc

def make_average(pinfty_array, ps, numerical_pc_array):
    '''
    Calculate the average of pinfty_array and set values below the percolation threshold to zero.

    Parameters
    ----------
    pinfty_array : array_like
        2D array of pinfty values.
    ps : array_like
        List of proportions of remaining nodes.
    numerical_pc_array : array_like
        List of numerical percolation thresholds.

    Returns
    -------
    pinfty_avg : array_like
        Averaged pinfty values with values below the percolation threshold set to zero.

    '''
    # Obtain the average of pinfty_array across the rows
    pinfty_avg = np.average(pinfty_array, axis=1)
    
    # Calculate the average numerical percolation threshold
    pc_avg = np.average(numerical_pc_array)
    
    # Find the index of the first ps value greater than the average numerical percolation threshold
    pc_pos = np.argwhere(ps > pc_avg)[0][0]
    
    # Set values of pinfty_avg below the percolation threshold to zero
    pinfty_avg[0:pc_pos] = 0
    
    return pinfty_avg

def transform2MTimes(connected_component, G):
    '''
    Transform the connected component into a matrix.

    Parameters
    ----------
    connected_component : list
        List of connected components.
    G : graph
        The graph.

    Returns
    -------
    G_matrix : array_like
        The matrix representation of the connected component.

    '''
    # Initialize the matrix with zeros
    max_node = max(G.nodes())
    G_matrix = np.zeros((max_node[0] + 1, max_node[0] + 1))
    
    # Fill the matrix with ones based on the connected components
    if len(connected_component) == 1:
        for index in list(connected_component[0]):
            x = int(index[0])
            y = int(index[1])
            G_matrix[x, y] = 1
    else:
        for each_cc in connected_component:
            for index in each_cc:
                x = int(index[0])
                y = int(index[1])
                G_matrix[x, y] = 1
    
    return G_matrix

def LocalizedAttackNodes(center_node, Rh, max_coordinate):
    '''
    Identify candidate nodes within a localized attack radius.

    Parameters
    ----------
    center_node : tuple
        Coordinates of the center node.
    Rh : float
        Radius of the attack.
    max_coordinate : int
        Maximum coordinate value.

    Returns
    -------
    candidatenodes : list
        List of candidate nodes within the attack radius.

    '''
    #create the hole e.g. make a circle, find all nodes within the circle  
    delta_coordinate = np.arange(0, Rh+1,1)
    candidate_coordinate = set()
    for coordinate in itert.product(delta_coordinate,delta_coordinate):
        if np.sqrt(np.power(coordinate[0],2) + np.power(coordinate[1],2)) <= Rh:
            [x,y] = [coordinate[0],coordinate[1]] 
            candidate_coordinate.add((x,y))
            candidate_coordinate.add((-x,y))
            candidate_coordinate.add((x,-y))
            candidate_coordinate.add((-x,-y))
            
    # Convert candidate coordinates to actual nodes
    candidatenodes = []
    for (dx,dy) in candidate_coordinate:
        nodex = center_node[0]+dx
        nodey = center_node[1]+dy
        [nodeux, nodeuy] = CheckPeriodicBoundary(nodex,nodey,max_coordinate)
        candidatenodes.append((nodeux, nodeuy))
    
    return candidatenodes

def LAKcorePercolation(G, Rh, k):
    '''
    Perform K-core percolation with a localized attack.

    Parameters
    ----------
    G : graph
        The network graph.
    Rh : float
        Radius of the localized attack.
    k : int
        The core number.

    Returns
    -------
    NOI : int
        Number of iteration steps.
    lcc_step_nodes_dict : dict
        Dictionary containing connected components at each step.

    '''
    # Randomly select a node as the center of the hole
    nodes = list(G.nodes())
    max_coordinate = max(nodes)[0]
    
    # Find the hole nodes
    center_node = random.choice(nodes)
    holenodes = LocalizedAttackNodes(center_node, Rh, max_coordinate)
    
    # Remove the nodes located in the hole 
    Gr = G.copy()
    Gr.remove_nodes_from(holenodes)
    degree_dict = {each[0]:each[1] for each in Gr.degree()}
    
    # Initialize variables
    NOI = 0
    lcc_step_nodes_dict = {}
    lcc_step_nodes_dict[NOI] = sorted(nx.connected_components(Gr), key=len, reverse=True)#initial connected component
    
    # Perform K-core percolation
    while len(degree_dict)!=0 and min(degree_dict.values()) < k:
        
        # Remove nodes with degree less than k
        nodes = list(Gr.nodes())
        for node in nodes:
            if degree_dict[node] < k:
                Gr.remove_node(node)
        
        # Update the degree dict
        degree_dict = {each[0]:each[1] for each in Gr.degree()}
        
        # Update the steps and connected components at each step
        NOI += 1
        lcc_step_nodes_dict[NOI] = sorted(nx.connected_components(Gr), key=len, reverse=True)
  
    return NOI, lcc_step_nodes_dict

def LocalizedAttackLcc(G, Rh, k):
    '''
    Perform K-core percolation with a localized attack and calculate the size of the largest connected component.

    Parameters
    ----------
    G : graph
        The network graph.
    Rh : float
        Radius of the localized attack.
    k : int
        The core number.

    Returns
    -------
    lcc : float
        Size of the largest connected component as a fraction of the total number of nodes.
    '''
    
    # Perform K-core percolation with a localized attack
    [NOI, lcc_step_nodes_dict]= LAKcorePercolation(G, Rh, k)
    
    # Calculate the size of the largest connected component
    if len(lcc_step_nodes_dict[NOI]) == 0:
        lcc = 0
    else:
        lcc = len(max(lcc_step_nodes_dict[NOI]))/G.order()
    return lcc

def BinarySearch(G, k, Rh_low, Rh_high, threshold):
    '''
    Perform binary search to find the localized attack radius for a given threshold of the largest connected component size.
    
    Parameters
    ----------
    G : graph
        The network graph.
    k : int
        The core number.
    Rh_low : float
        Lower bound of the attack radius.
    Rh_high : float
        Upper bound of the attack radius.
    threshold : float
        Threshold value for the largest connected component size.
    
    Returns
    -------
    Rh_mid : float
        Attack radius that satisfies the threshold.
    
    '''

    #mid_Rhs = {}
    t = 0

    # Calculate the largest connected component size for the lower and upper bounds
    lcc_low = LocalizedAttackLcc(G, Rh_low, k)
    print('lower Rh:%d and lcc:%f at %d step'%(Rh_low, lcc_low, t))
    
    lcc_high = LocalizedAttackLcc(G, Rh_high, k)
    print('higher Rh:%d and lcc:%f at %d step'%(Rh_high, lcc_high, t))
    
    # Perform binary search
    while Rh_low <= Rh_high:
        
        #calculate the middle value 
        Rh_mid = int(round((Rh_low + Rh_high)/2,0))
        lcc_mid = LocalizedAttackLcc(G, Rh_mid, k)
        print('Mid Rh:%d and lcc:%f at %d step'%(Rh_mid, lcc_mid, t))
        t = t+1
        
        # Check conditions to update bounds or return the result
        if (Rh_high-Rh_low) == 0: #return the mid value
            return Rh_mid
        elif lcc_low < threshold:
            return  Rh_low
        elif lcc_high > threshold:
            return  Rh_high
        elif lcc_mid < threshold: #lower rh turn to be larger if lcc_low > threshold
            Rh_high = Rh_mid-1
            lcc_high = LocalizedAttackLcc(G, Rh_high, k)
            print('higher Rh:%d and lcc:%f at %d step'%(Rh_high, lcc_high, t))
        elif lcc_mid > threshold: #higher rh turn to be larger  if lcc_high < threshold
            Rh_low = Rh_mid+1 
            lcc_low = LocalizedAttackLcc(G, Rh_low, k)
            print('lower Rh:%d and lcc:%f at %d step'%(Rh_low, lcc_low, t))
    
    return Rh_mid

def BinarySearchPc(G, kcore, threshold):
    
    '''
    Perform binary search to find the percolation threshold for a given threshold of the largest connected component size.
   
    Parameters
    ----------
    G : graph
        The network graph.
    kcore : int
        The core number.
    threshold : float
        Threshold value for the largest connected component size.
   
    Returns
    -------
    p : float
        Percolation threshold that satisfies the threshold.
   
    '''

    # Set the initial values for the search range
    pleft = 0
    pright = 1
    N = G.order() # Total number of nodes in the network
    
    # Initialize time steps
    t = 0
    
    # Perform binary search
    while (pright-pleft) > 1/N:
        
        # Calculate the middle value
        pmid = (pleft + pright)/2
        lcc_pmid = kcorepercolation(G,pmid,kcore)[0]/N
        
        # Update time steps
        t = t+1

        if lcc_pmid < threshold:
            # Move the left boundary to the middle
            pleft = pmid
        elif lcc_pmid > threshold: 
            # Move the right boundary to the middle
            pright = pmid
    
    p = pright
    
    return p

def RunSimulationBinaryPc(args):
    '''
    Run simulations to find the percolation threshold using binary search.
   
    Parameters
    ----------
    args : tuple
        Tuple containing the following:
            G : graph
                The network graph.
            kcore : int
                The core number.
            threshold : float
                Threshold value for the largest connected component size.
            PT : float
                Percolation threshold obtained from previous iteration.
   
    Returns
    -------
    p : float
        Percolation threshold that satisfies the threshold.
    
    '''
    
    [G, kcore, threshold, PT]= args 
    # Find the percolation threshold using binary search
    p = BinarySearchPc(G, kcore, threshold)
    
    return p
