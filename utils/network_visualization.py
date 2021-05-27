# Plotting functions for Multiplex Diffusion

# Zhao Song, April 7, 2020

from pylab import *
from numpy import *
from numpy import matlib

from network_properties import *
from eigenvalues import *


def network_visualization(G,inter_layer_edge_list,intra_layer_edge_colors):
    """
        Parameters
        ----------
        G: a list of a network properties
        inter_layer_edge_list:
        intra_layer_edge_colors: 
        """
    n = G[0] # number of nodes in each layer
    t = G[1] # number of layers
    
    pos,pos_relabel = get_node_position(G)
    Gm = nx.DiGraph()
    Gm.add_nodes_from(pos_relabel)
    
    fig,ax = subplots(1,1,figsize = (7,5))
    
    nx.draw_networkx_nodes(Gm, pos=pos_relabel,
                           node_color = 'k',
                           alpha = .8,
                           node_size = 100,
                           node_label = True,
                           edgecolors = 'k'
                           )

    for i_t in range(t):
        nx.draw_networkx_edges(G[3][i_t],
                              pos=pos['layer '+str(i_t+1)],
                              edgelist = G[3][i_t].edges(),
                              alpha = .6,
                              edge_color = intra_layer_edge_colors[i_t],
                              );
#        nx.draw_networkx_labels(G[3][i_t], pos=pos['layer '+str(i_t+1)], font_size=20, font_family="sans-serif")

#    nx.draw_networkx_edges(Gm,
#                           pos=pos_relabel,
#                           edgelist = inter_layer_edge_list,
#                           alpha = .3,
#                           );

    plt.arrow(-.3, -1.9, 0, .7,
     #         -.2, -3.5, 0, 2,
               head_width = 0.055,head_length = 0.2,
               width = 0.0000005,
               ec ='k',
              facecolor = 'k',linestyle='dashed')
    plt.arrow(.24, -.9, 0, -.7,
     #         .15, -1.5, 0, -2,
               head_width = 0.055,head_length = 0.2,
               width = 0.003,
               ec ='k',
              facecolor = 'k')
    fig.text(.25, .48, '$(1-\delta)\omega$',fontsize=14)
    fig.text(.66, .49, '$(1+\delta)\omega$',fontsize=14)
    
    
#     fig.text(.87, .68, 'Human',fontsize=14)
#     fig.text(.9, .3, 'AI',fontsize=14)
    
    #     x = [-.9,.65,.85,-.7]
    #     y1 = [-1.2,-1.2,1.2,1.2]
    #     ax.add_patch(patches.Polygon(xy=list(zip(x,y1)),facecolor=edge_colors[0],alpha=.2 ))
    #     y2 = [-4.2,-4.2,-1.8,-1.8]
    #     ax.add_patch(patches.Polygon(xy=list(zip(x,y2)),facecolor=edge_colors[1],alpha=.2 ))
    
    plt.box(False)



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
    pos = {}
    #     pos = nx.spring_layout(G1) # nodes positions of layer 1
    #     pos = nx.shell_layout(G1) # nodes positions of layer 1
    #     pos = nx.circular_layout(G1) # nodes positions of layer 1
    pos['layer 1'] = nx.kamada_kawai_layout(G1) # nodes positions of layer 1
    pos_relabel = nx.kamada_kawai_layout(G1)
    
    # nodes positions of remaining layers
    for i_t in range(1,t):
        #        pos['layer '+str(i_t+1)] = pos['layer '+str(i_t)]
        pos['layer '+str(i_t+1)] = {}
        for i in range(1,n+1):
            xy = pos['layer '+str(i_t)][i]
            pos['layer '+str(i_t+1)][i] = array([xy[0],xy[1]-3])
            pos_relabel[i+(i_t)*n] = array([xy[0],xy[1]-3])
    #            pos[i_t+1][i] = array([xy[0]-.1,xy[1]-3])
#            pos[i_t+1][i] = array([xy[0],xy[1]-5])

    return pos,pos_relabel
