from pylab import *
import networkx as nx
from numpy import *


def get_Laplacian(G,delta):
    """
    Parameters
    ----------
    G: a list of a network properties
    delta: asymmetry parameter

    Returns
    -------
    Ls: dictionary
    Laplacian matrix dictionary for different Laplacian matrices
    (intralayer and interlayer (supra)Laplacians)
    """
    
    [n,t,nt,GLs,EI] = G
    Ls ={}
    
#### Intralayer Laplacian: LT ####
#### Intralayer supra-Laplacian: LL ####
    Ls['LTs'] = zeros((t,n,n)) # the intralayer Laplacian t * (n x n)
    Ls['ATs'] = zeros((t,n,n)) # the intralayer Adjacency Mtrices t * (n x n)
    Ls['LL'] = zeros((nt,nt)) # LL: supra-Laplacian of the independent layers
    for i in range(t):
        AT = nx.to_numpy_array(GLs[i]); # AI:intralayer adjacency matrix
        AT = AT.T
        d = sum(AT,1);DT = diag(d) # DI: intralayer degree matrix
        LT = DT - AT
        Ls['LTs'][i,:,:] = LT # LTs contains t intralayer Laplacians LT
        Ls['ATs'][i,:,:] = AT # ATs contains t intralayer adjacency matrices AT
        Ls['LL'][i*n:(i+1)*n,i*n:(i+1)*n] = LT # Intralayer supra-Laplacian

#### Interlayer Laplacian: II ####
    AI = zeros((t,t)) # AI:interlayer adjacency matrix
    for edge in EI:
        i,j = edge-1
        AI[i,j] = 1+delta
        AI[j,i] = 1-delta
    AI = AI.T
    d = sum(AI,1)
    DI = diag(d) # DI: interlayer degree matrix
    Ls['II'] = DI - AI # Interlayer Laplacian

#### Interlayer supra-Laplacian: LI ####
    I = eye(n)
    Ls['LI'] = kron(Ls['II'],I) # Kronecker product. LI: interlayer supra-Laplacian
    
    return Ls



def get_Laplacian_2(G,delta,chi):
    """
        It only works for 2-layer networks. Since the intralayer Laplacians have been rate scaled:
            L1 -> chi * L1
            L2 -> (1-chi) * L2
        
        Parameters
        ----------
        G: a list of a network properties
        delta: asymmetry parameter
        chi: rate-scaling parameter
        
        Returns
        -------
        Ls: dictionary
            Laplacian matrix dictionary for different Laplacian matrices
            (intralayer and interlayer (supra)Laplacians)
        """
    
    [n,t,nt,GLs,EI] = G
    Ls ={}
    
    #### Intralayer Laplacian: LT ####
    #### Intralayer supra-Laplacian: LL ####
    Ls['LTs'] = zeros((t,n,n)) # the intralayer Laplacian t*(n x n)
    Ls['LL'] = zeros((nt,nt)) # LL: supra-Laplacian of the independent layers
    for i in range(t):
        AT = nx.adjacency_matrix(GLs[i]).toarray();
        AT = AT.T;d = sum(AT,1);DT = diag(d)
        LT = DT - AT
        Ls['LTs'][i,:,:] = LT # LTs contains t intralayer Laplacian
        if i == 0:
            Ls['LL'][i*n:(i+1)*n,i*n:(i+1)*n] = chi * LT
        elif i == 1:
            Ls['LL'][i*n:(i+1)*n,i*n:(i+1)*n] = (1-chi) * LT

    ### Interlayer Laplacian: II ####
    AI = zeros((t,t)) # interlayer adjacency matrix
    for edge in EI:
        i,j = edge-1
        AI[i,j] = 1+delta
        AI[j,i] = 1-delta
    AI = AI.T
    d = sum(AI,1)
    DI = diag(d) # interlayer degree matrix
    Ls['II'] = DI - AI

    #### Interlayer supra-Laplacian: LI ####
    I = eye(n)
    Ls['LI'] = kron(Ls['II'],I) # Kronecker product. LI: interlayer supra-Laplacian

    return Ls



def get_supra_Laplacian(G,omega,delta):
    """
    Parameters
    ----------
    G: a list of a network properties
    omega: coupling strength
    delta: asymmetry parameter
    
    Returns
    -------
    supra_L: supra-Laplacian matrix
    """
    
    Ls = get_Laplacian(G,delta)
    supra_L = Ls['LL'] + omega*Ls['LI']
               
    return supra_L











