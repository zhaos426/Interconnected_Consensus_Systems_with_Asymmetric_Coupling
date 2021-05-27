# Functions for Multiplex Diffusion

# Zhao Song, Feb 26, 2021

from pylab import *
#import networkx as nx
from numpy import *

from network_properties import *
from Laplacians import *

def get_eigen(L):
    """
    Parameters
    ----------
    L: a matrix
    
    Returns
    -------
    evals: sorted eigenvalues of matrix L from small to large
    evecs: eigenvectors of matrix L ordered as corresponding eigenvalues
    """
    evals,evecs = eig(L)
    ids = argsort(real(evals))
    evecs = evecs[:,ids]
    evals = evals[ids]
    
    return evals,evecs



def get_lambda2(L):
    """
    Parameters
    ----------
    L: a matrix
    
    Returns
    -------
    the second smallest eigenvalue of matrix L
    """
    
    return get_eigen(L)[0][1]



def get_evals_across_omegas_deltas(G,omegas,deltas):
    """
    Parameters
    ----------
    G: a list of a network's properties
    omegas: a bunch of coupling strengths
    deltas: a bunch of asymmetry parameters
    
    Returns
    -------
    Levals: eigenvalues of supra-Laplacian
    """
    
    [n,t,nt,GLs,GI] = G
    Levals = zeros((nt,len(omegas),len(deltas) ),dtype=complex) # eigenvalues of supra-Laplacian
    Ls = [get_Laplacian(G,delta=0)] * len(deltas) # Dictionary with 4 kinds of Laplacian 'LTs','LL','II','LI'
    L = zeros((nt,nt,len(deltas))) # N X N supra Laplacian
    for i_delta,delta in enumerate(deltas):
        Ls[i_delta] = get_Laplacian(G,delta)
        for i_omega,omega in enumerate(omegas):
            L[:,:,i_delta] = Ls[i_delta]['LL'] + omega*Ls[i_delta]['LI']
            Levals[:,i_omega,i_delta] = get_eigen(L[:,:,i_delta])[0]

    return Levals



def get_lambda2_thm(G,deltas):
    """
        Parameters
        ----------
        G: a list of a network's properties
        deltas: a bunch of asymmetry parameters
        
        Returns
        -------
        lambda2s: predicted smallest nonzero eigenvalue of supra-Laplacian
        """
    n = G[0]
    t = G[1]
    lambda2s = zeros((len(deltas),2),dtype = complex)
    for i_delta,delta in enumerate(deltas):
        Ls = get_Laplacian(G,delta)
        u = get_eigen(Ls['II'].T)[1][:,0] # the left eigenvector associated with zero eigenvalue of interlayer Laplacian
        u_normlized = u/sum(u)
        L_bar = zeros((n,n))
        for i_t in range(t):
            L_bar = L_bar + u_normlized[i_t] * Ls['LTs'][i_t]
        lambda2s[i_delta,:] = get_eigen(L_bar)[0][1:3]
    return lambda2s



def get_Leval2_2layer_eachlayer(network_name,chis):
    """
    Parameters
    ----------
    network_name: the name of the network
    chis: a bunch of rate-scaling parameters
    
    Returns
    -------
    Leval2: the second smallest eigenvalue of individual layer Laplacian
    Ls: dictionary
        Laplacian matrix dictionary for different Laplacian matrices
        (intralayer and interlayer (supra)Laplacians)
    """
    G = get_network_properties(network_name)
    Ls = get_Laplacian(G,0)
    Leval2 = zeros((4,len(chis)),dtype = complex)
    for i,chi in enumerate(chis):
        Leval2[0,i] = get_eigen(chi * Ls['LTs'][0])[0][1] # the second smallest eigenvalue of chi*L1
        Leval2[2,i] = get_eigen((1 - chi) * Ls['LTs'][1])[0][1] # the second smallest eigenvalue of (1-chi)*L2
    Leval2[1,:] = get_eigen(Ls['LTs'][0])[0][1] # the second smallest eigenvalue of L1
    Leval2[3,:] = get_eigen(Ls['LTs'][1])[0][1] # the second smallest eigenvalue of L2
    return Leval2,Ls




def get_lambda2_2layer_thm(network_name,deltas,chis):
    """
    Parameters
    ----------
    network_name: the name of the network
    deltas: a bunch of asymmetry parameter
    chis: a bunch of rate-scaling parameters
    
    Returns
    -------
    lambda2s: the predicted second smallest eigenvalue of supra-Laplacian with rate-scaling parameter chi.
    """
    
    G = get_network_properties(network_name)
    Ls = get_Laplacian(G,0)
    lambda2s = zeros((len(deltas),len(chis)),dtype = complex)
    for i,delta in enumerate(deltas):
        for j,chi in enumerate(chis):
            L_bar = (1+delta)/2 * chi * Ls['LTs'][0] + (1-delta)/2 * (1 - chi) * Ls['LTs'][1]
            lambda2s[i,j] = get_eigen(L_bar)[0][1]
    return lambda2s


def get_lambda2_2layer_exact(network_name,omegas,deltas,chis):
    """
    Parameters
    ----------
    network_name: the name of the network
    deltas: a bunch of asymmetry parameter
    chis: a bunch of rate-scaling parameters
    
    Returns
    -------
    lambda2s: the exact second smallest eigenvalue of supra-Laplacian
    """
    G = get_network_properties(network_name)
    lambda2s = zeros((len(omegas),len(deltas),len(chis)),dtype = complex)
    for i,omega in enumerate(omegas):
        for j,delta in enumerate(deltas):
            for k,chi in enumerate(chis):
                Ls = get_Laplacian_2(G,delta,chi)
                lambda2s[i,j,k] = get_eigen(Ls['LL'] + omega*Ls['LI'])[0][1]
    return lambda2s



def get_dlambda2_ddelta_thm(network_name,chis):
    """
        Parameters
        ----------
        network_name: the name of the network
        chis: a bunch of rate-scaling parameters
        
        Returns
        -------
        -dlambda2_ddelta: the predicted derivative of the second smallest eigenvalue of supra-Laplacian
        """
    G = get_network_properties(network_name)
    Ls = get_Laplacian(G,0)
    dlambda2_ddelta = zeros((len(chis),2),dtype = complex)
    for i,delta in enumerate([-1,1]):
        lambda2s_star = get_eigen(Ls['LTs'][i])[0]
    for i,delta in enumerate([-1,1]):
        u = get_eigen(Ls['LTs'][i])[1][:,1]
        v = get_eigen(Ls['LTs'][i].T)[1][:,1]
        for i_chi,chi in enumerate(chis):
            L = - ( chi * Ls['LTs'][0] - (1-chi) * Ls['LTs'][1] )/2
            if dot(v.T.conjugate(),u) != 0:
                dlambda2_ddelta[i_chi,i] = dot(dot(v.T.conjugate(),L),u) / dot(v.T.conjugate(),u)
    return -dlambda2_ddelta



def get_optima(network_name,omegas,deltas,chis):
    optima_delta = zeros(len(chis))
    lambda2s_thm = get_lambda2_2layer_thm(network_name,deltas,chis) #delta*chi
    for i_chi,chi in enumerate(chis):
        optima = 0
        i_optima = 0
        for i_lam,lam in enumerate(lambda2s_thm[:,i_chi].real):
            if lam > optima:
                optima = lam
                i_optima = i_lam
        optima_delta[i_chi] = deltas[i_optima]
    return optima_delta

