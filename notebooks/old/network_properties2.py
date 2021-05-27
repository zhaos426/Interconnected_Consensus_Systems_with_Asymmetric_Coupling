from pylab import *
from numpy import *
import networkx as nx

def get_network_properties(network_name):
    """
    Parameters
    ----------
    network_name: the name of a graph model
    
    Returns
    -------
    a list of network properties:
        n: the number of nodes in each layer
        t: the number of layers
        nt: the total number of nodes in a network
        GLs: a list of t individual NetworkX Graphs(directed and/or undirected)
        EI: interlayer edges
    """
    
    ##############################
    ## Initialization: 2 layers ##
    ##############################
    if network_name == 'SimpleNetwork':
        n = 6 # n: number of nodes in one layer
        t = 2 # m: number of layers
        nt = n*t
        
        # Initialization: intralayer connection
        EL1 = array([[1,4],[1,5],[1,6],[2,4],[2,6],[3,4],[3,6],[4,6],[5,6]]) # edges in the first layer
        GL1 = nx.Graph()
        GL1.add_edges_from(EL1)
        
        EL2 = array([[1,4],[1,5],[1,6],[2,3],[2,5],[3,5],[5,6]]) # edges in the second layer
        GL2 = nx.Graph()
        GL2.add_edges_from(EL2)
        GLs = [GL1,GL2]
        
        # Initialization: interlayer connection
        EI = array([[1,2]]) # eg: the element[1,2] means layer 1 is connected with layer 2

    ##################################################################
    ## Initialization: 2 layer networks with Opposite directed ring ##
    ##################################################################
    if network_name == 'RingGraph':
        n = 6 # n: number of nodes in one layer
        t = 2 # t: number of layers
        nt = n*t
        
        #### Initialization: intralayer connection ####
#        GL1 = nx.generators.classic.cycle_graph(6, create_using=nx.DiGraph()) # a directed ring graph
        EL1 = array([[6,5],[5,4],[4,3],[3,2],[2,1],[1,6]]) # edges in the first layer
        GL1 = nx.DiGraph()
        GL1.add_edges_from(EL1)
        
        GL2 = nx.reverse_view(GL1) # reverse direction with G1
        GLs = [GL1,GL2]
        
        #### Initialization: interlayer connection ####
        EI = array([[1,2]]) # eg: the element[1,2] means layer 1 is connected with layer 2


    if network_name == 'RingGraph4':
        n = 6 # n: number of nodes in one layer
        t = 2 # t: number of layers
        nt = n*t
        
        #### Initialization: intralayer connection ####
        EL1 = array([[2,1],[3,2],[6,3],[1,6],[6,5],[5,4],[4,3]]) # edges in the first layer
        GL1 = nx.DiGraph()
        GL1.add_edges_from(EL1)

        EL2 = array([[1,2],[2,3],[3,6],[6,1],[5,6],[4,5],[3,4]]) # edges in the second layer
        GL2 = nx.DiGraph()
        GL2.add_edges_from(EL2)
        GLs = [GL1,GL2]
        
        #### Initialization: interlayer connection ####
        EI = array([[1,2]]) # eg: the element[1,2] means layer 1 is connected with layer 2


    if network_name == 'RingGraph6':
        n = 6 # n: number of nodes in one layer
        t = 2 # t: number of layers
        nt = n*t
        
        #### Initialization: intralayer connection ####
        EL1 = array([[2,1],[3,2],[6,3],[1,6],[3,4],[4,5],[5,6]]) # edges in the first layer
        GL1 = nx.DiGraph()
        GL1.add_edges_from(EL1)
        
        EL2 = array([[1,2],[2,3],[3,6],[6,1],[4,3],[5,4],[6,5]]) # edges in the second layer
        GL2 = nx.DiGraph()
        GL2.add_edges_from(EL2)
        GLs = [GL1,GL2]
        
        #### Initialization: interlayer connection ####
        EI = array([[1,2]]) # eg: the element[1,2] means layer 1 is connected with layer 2


    if network_name == 'LineStar':
        n = 6 # n: number of nodes in one layer
        t = 2 # t: number of layers
        nt = n*t
        
        #### Initialization: intralayer connection ####
        EL1 = array([[1,2],[2,3],[3,4],[4,5],[5,6]]) # edges in the first layer
        GL1 = nx.Graph()
        GL1.add_edges_from(EL1)
        
        EL2 = array([[1,2],[1,3],[1,4],[1,5],[1,6]]) # edges in the second layer
        GL2 = nx.Graph()
        GL2.add_edges_from(EL2)
        GLs = [GL1,GL2]
        
        #### Initialization: interlayer connection ####
        EI = array([[1,2]]) # eg: the element[1,2] means layer 1 is connected with layer 2


    if network_name == 'Krackhardt-High-Tech':
        
        #https://comunelab.fbk.eu/data.php
        n = 21 # n: number of nodes in one layer
        t = 2 # m: number of layers
        nt = n*t
        edges_file_name = 'Krackhardt-High-Tech_Multiplex_Social' + '/Dataset/' + 'Krackhardt-High-Tech_multiplex.edges'
        nodes_file_name = 'Krackhardt-High-Tech_Multiplex_Social' + '/Dataset/' + 'Krackhardt-High-Tech_nodes.txt'
        
        #### Initialization: interlayer connection ####
        EI = array([[1,2]]) # eg: the element[1,2] means layer 1 is connected with layer 2
        
        #### Initialization: intralayer connection ####
        EL = file_to_edges(edges_file_name)
        EL1 = [] # edges in the first layer
        EL2 = [] # edges in the second layer
        for i in range(len(EL)):
            if EL[i][0] == 1:
                EL1.append(EL[i][1:3])
            elif EL[i][0] == 2:
                EL2.append(EL[i][1:3])

        nodes_list = [i+1 for i in range(n)]
        GL1 = nx.DiGraph()
        GL1.add_nodes_from(nodes_list)
        GL1.add_edges_from(EL1)

        GL2 = nx.DiGraph()
        GL2.add_nodes_from(nodes_list)
        GL2.add_edges_from(EL2)

        GL1 = nx.reverse_view(GL1)
        GL2 = nx.reverse_view(GL2)


        # make a layer with random directed edges

        random.seed(1)
        #best so far (d,seed) =  (4,0)
#（5，113）

        d = 4 # node degree
        EL2b = zeros( (n*d,2 ),dtype=int) # edge list
        count = 0
        for i in range(n):
            neighbors = np.random.permutation(n)[:d]# for node i, choose 3 neighbors at random
            for j in neighbors:
                EL2b[count] = [i,j]
                count +=1
        GL2b = nx.DiGraph()
        GL2b.add_edges_from(EL2b)

        #GLs = [GL1,GL2]
        GLs = [GL1,GL2b]

    return [n,t,nt,GLs,EI]



def file_to_edges(file_name):
    """
    from edges file to get the intralayer edge list of the social network
    """
    
    edges_file = open('../RealWorldNetworks/' + file_name,'r')
    edges_file = edges_file.read().split('\n')
    edges = []
    for e in edges_file:
        e = e.split(' ')
        if len(e) == 4:
            e = [int(i) for i in e]
            edges.append(e)

    return edges



def get_networks_properties(network_names):
    """
    Parameters
    ----------
    network_names: a list contains several names of graph models
    
    Returns
    -------
    Gs: graph models dictionary
    """
    
    Gs = {}
    for network_name in network_names:
        Gs[network_name] = get_network_properties(network_name)
    
    return Gs



def get_intra_layer_edges(edges_file_name,n):
    edges = file_to_edges(edges_file_name)
    ELs = [[],[],[]]
    for i in range(len(edges)):
        i_t = edges[i][0]-1
        ELs[i_t].append(tuple(array(edges[i][1:3])+i_t*n))
    return ELs



def get_node_position(G):
    """
    Parameters
    ----------
    G: a list of a network properties

    Returns
    -------
    pos: all nodes positions of a network
    """
    n = G[0]
    t = G[1]
    G1 = G[3][0]
    #     pos = nx.spring_layout(G1) # nodes positions of layer 1
    #     pos = nx.shell_layout(G1) # nodes positions of layer 1
    #     pos = nx.circular_layout(G1) # nodes positions of layer 1
    pos = nx.kamada_kawai_layout(G1) # nodes positions of layer 1
    
    # nodes positions of remaining layers
    for i_t in range(t-1):
        for i in range(1,n+1):
            xy = pos[i+(i_t)*n]
            pos[i+(i_t+1)*n] = array([xy[0],xy[1]-3])
#            pos[i+(i_t+1)*n] = array([xy[0]-.1,xy[1]-3])
#            pos[i+(i_t+1)*n] = array([xy[0],xy[1]-5])

    return pos




