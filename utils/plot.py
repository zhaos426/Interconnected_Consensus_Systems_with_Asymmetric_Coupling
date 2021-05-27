# Plotting functions for Multiplex Diffusion

# Zhao Song, April 7, 2020

from pylab import *
from numpy import *
from numpy import matlib

from network_properties import *
from eigenvalues import *



##########################################################
## Choosing colors
##########################################################
def get_color1(i_delta,deltas):
    """
    color range from purple to green
    """
    return ( .6,.4 + i_delta/len(deltas)/2  ,(len(deltas)-i_delta)/len(deltas)/2 + .5   )

def get_color2(i_omega,omegas):
    """
    color range from orange to blue
    """
    return ( (len(omegas) - i_omega)/len(omegas)  ,.6,i_omega/len(omegas) )



######################################################
## fig 1 ##
######################################################
def plot_timeseries(ax,G,deltas,omegas,tt,dt):
    """
        Parameters
        ----------
        ax: axes
        G: a list of a network's properties
        deltas: a bunch of asymmetry parameters
        omegas: a bunch of coupling strengths
        tt: length of time series
        dt: time step
        
        
        Returns
        -------
        """
    
    supraL = get_supra_Laplacian(G,omegas[0],deltas[0])
    #supraL = supraL.T
    nt = len(supraL)
    T = 2
    n = int(nt/ T)
    
    x0 = np.random.rand(nt)
    x0[:n] = x0[:n] + 2
    #x0[n:] = x0[n:] - .5
    #x0 = np.linspace(-.5,.5,nt)
    x0 = x0 - mean(x0)
    #x0 = x0 - mean(x0) + 1
    #     x0 = x0 / sum(x0)
#    print((x0))

    errors = np.zeros((tt,len(deltas)))
#     colors = ['#ff7f0e','#2ca02c','#8c564b']
#     colors = ['#9467bd','#8c564b','#bcbd22']
    colors = ['#ff7f0e','#2ca02c','#bcbd22']
    for i,delta in enumerate(deltas):
        supraL = get_supra_Laplacian(G,omegas[i],deltas[i])
        #         supraL = supraL.T
        lam,uu = eig(supraL.T)
        idd = argsort(real(lam))[0]
        lambda1 = lam[idd]
        x_final = abs(real(uu[:,idd]))
        x_final = x_final / sum(x_final)
        x_final = dot(x_final,x0)
        
        x = simulate_collective_dynamics(supraL,x0,dt,tt)
        
        ax[i].set_xlabel('time, $t$, ($\delta = %s $)' %delta)
        ax[i].set_ylabel('$x_p(t)$')
        ax[i].plot(dt*arange(tt),dt*arange(tt)*0,'k:')
            
        for k in range(n):
            ax[i].semilogx(dt*arange(tt),x.T[n+k] ,color=[.4,.7,.9],alpha=.7)      #blue
            ax[i].semilogx(dt*arange(tt),x.T[k] ,color=[.8,.1,.1],alpha=.5)       #red
        
        #ax[i].set_ylim([-.5,.5])
        y = x.copy()
        for t in range(tt):
            y[t] = x[t] - x_final #- np.mean(x[t])
            errors[t,i] = np.linalg.norm(y[t],2)
        ax[3].semilogy(dt*arange(tt),errors[:,i],color=colors[i])
    ax[3].set_ylabel('$||\mathbb{x}(t) - \overline{\mathbb{x}} ||$')
    ax[3].set_xlabel('time, $t$')
    ax[3].legend(['$\delta = '+str(deltas[i])+'$' for i in range(3)])
    return errors,x,x_final,x0,lam,uu,idd,supraL



def simulate_collective_dynamics(supraL,x0,dt,tt):
    x = matlib.repmat(x0,tt,1)
    for t in range(tt-1):
        x[t+1] = x[t] - dt * np.dot(supraL, x[t])
    #x[t+1] = x[t+1]/sum(x[t+1])
    return x



##########################################################
## fig 2 ##
##########################################################
def plot_lambda2_vs_omega_deltas(G,omegas,deltas,axis):
    """
        plotting the real part of second smallest eigenvalue of supra-Laplacian versus omega for different choices of delta.
        
        Parameters
        ----------
        G: a list of a network's properties
        omegas: a bunch of coupling strengths
        deltas: a bunch of asymmetry parameters
        axis
        """
    Levals = get_evals_across_omegas_deltas(G,omegas,deltas) # eigenvalues of supra-Laplacian
    for i_delta in range(len(deltas)):
        axis.semilogx(omegas,Levals[1,:,i_delta].real,alpha=.9,color=get_color1(i_delta,deltas) ,linewidth=1)



def plot_lambda2_vs_delta_omegas(Levals,omegas,deltas,axis):
    """
        plot the real part of second smallest eigenvalue of supra-Laplacian versus delta for different choices of omega.
        
        Parameters
        ----------
        Levals: the real part of second smallest eigenvalue of supra-Laplacian
        omegas: a bunch of coupling strengths
        deltas: a bunch of asymmetry parameters
        axis
        """
    
    for i_omega in range(len(omegas)):
        axis.plot(deltas,Levals[1,i_omega,:],alpha=.8,color=get_color2(i_omega,omegas) ,linewidth=1)



def plot_lambda2_vs_delta_thm(G,deltas,ax):
    """
        plot the real part of predicted second smallest eigenvalue of supra-Laplacian versus delta for large omega.
        
        Parameters
        ----------
        G: a list of a network's properties
        deltas: a bunch of asymmetry parameters
        """
    
    lambda2s_thm = get_lambda2_thm(G,deltas)
    ax.plot(deltas,lambda2s_thm[:,0].real,'k',linewidth=2)



######################################################
## fig 3A ##
######################################################
def plot_lambda2_vs_delta_chis(G,deltas,chis,plot_chis,axis):
    """
        plotting the predicted convergence rate versus delta for different choices of chis
        
        Parameters
        ----------
        G: a list of a network's properties
        deltas: a bunch of asymmetry parameters
        chis: a bunch of rate-scaling parameters
        plot_chis: a bunch of index of chis you want to plot
        axis
        """
    
    lambda2s_thm = get_lambda2_2layer_thm(G,deltas,chis) # delta*chi
    for i in plot_chis:
        axis.plot(deltas, lambda2s_thm[:,i].real,color=get_color1(i,chis))
    axis.set_xlabel('asymmetry, $\delta$')
    axis.set_aspect(1.0/axis.get_data_ratio(),adjustable='box')


######################################################
## fig 3B ##
######################################################
def plot_dlambda2_ddelta(G,omegas,deltas,chis,axis):
    """
        Parameters
        ----------
        G: a list of a network's properties
        omegas:a bunch of coupling strengths
        deltas: a bunch of asymmetry parameters
        chis: a bunch of rate-scaling parameters
        """
    
    lambda2s_exact = get_lambda2_2layer_exact(G,omegas,deltas,chis) #omega*delta*chi
    dlambda2_ddelta_exact = (lambda2s_exact[:,1:,:] - lambda2s_exact[:,0:-1,:]).T / (deltas[1] - deltas[0])
    dlambda2_ddelta_thm = get_dlambda2_ddelta_thm(G,chis) #chi*2
    
    for i in range(len(omegas)):
        axis.plot(chis,dlambda2_ddelta_exact[:,0,i].real,'-',alpha=.8,color=get_color2(i,omegas))

    for i in range(len(omegas)):
        axis.plot(chis,dlambda2_ddelta_exact[:,-1,i].real,'--',alpha=.8,color=get_color2(i,omegas))

    axis.set_xlabel('rate-scaling parameter, $\chi$')
    axis.plot(chis,chis*0,'k:')
    axis.set_aspect(1.0/axis.get_data_ratio(),adjustable='box')
#    index1 = arange(-2*(int(len(chis)/30)+1),2*(int(len(chis)/30)+1))
#    index2 = arange(-2*(int(len(chis)/20)+1),3*(int(len(chis)/20)+1))
    chi_1 = where(sign(dlambda2_ddelta_thm[:-1,0]) != sign(dlambda2_ddelta_thm[1:,0]))[0]
    chi_2 = where(sign(dlambda2_ddelta_thm[:-1,-1]) != sign(dlambda2_ddelta_thm[1:,-1]))[0]
    
    if chi_1.size != 0 and chi_2.size != 0:
        axis.plot(chis,dlambda2_ddelta_thm[:,0].real,'--k')
        axis.plot(chis,dlambda2_ddelta_thm[:,-1].real,'-k')
#         index1 = index1 + chi_1[0]
#         index2 = index2 + chi_2[0]
#         axis.plot(chis[index1],dlambda2_ddelta_thm[index1,0].real,'--k',linewidth=2.5)
#         axis.plot(chis[index2],dlambda2_ddelta_thm[index2,-1].real,'-k',linewidth=2.5)
        return chis[chi_2[0]],chis[chi_1[0]]

    elif chi_1.size != 0 and chi_2.size == 0:
        axis.plot(chis,dlambda2_ddelta_thm[:,0].real,'--k')
#         index1 = index1 + chi_1[0]
#         axis.plot(chis[index1],dlambda2_ddelta_thm[index1,0].real,'--k',linewidth=2.5)
    #         axis.plot(chis[chi_1[0]] + chis*0,linspace(-2,2,len(chis)),'r')
        return 0,chis[chi_1[0]]
    
    elif chi_1.size == 0 and chi_2.size != 0:
        axis.plot(chis,dlambda2_ddelta_thm[:,-1].real,'-k')
#         index2 = index2 + chi_2[0]
#         axis.plot(chis[index2],dlambda2_ddelta_thm[index2,-1].real,'-k',linewidth=2.5)
        #         axis.plot(chis[chi_2[0]] + chis*0,linspace(-2,2,len(chis)),'g')
        return chis[chi_2[0]],1
    
    else:
        return


######################################################
## fig 3C ##
######################################################
def plot_optima(G,omegas,deltas,chis,axis):
    optima_delta = get_optima(G,omegas,deltas,chis)
    axis.plot(chis,optima_delta)
    axis.plot([0,1],[0,0],':k')
    axis.set_xlabel('rate-scaling parameter, $\chi$')



def plot_contour(G,deltas,chis,axis):
    lambda2s_thm = get_lambda2_2layer_thm(G,deltas,chis) # delta*chi
    contours = axis.contour(chis, deltas, lambda2s_thm.real, 10, colors='black',linewidths=.5)
    axis.clabel(contours,
                fontsize=8)
    im = axis.imshow(lambda2s_thm.real,
                     extent=[0, 1, -1, 1],
                     origin='lower',
                     cmap='viridis',
                     vmin=0,vmax=1.5
                     )
    axis.set_xlabel('rate-scaling parameter, $\chi$')
    axis.set_aspect(1/2)
    return im



######################################################
## fig 3D ##
######################################################
def plot_Leval2_eachlayer(G,chis,axis):
    """
    Parameters
    ----------
    G: a list of a network's properties
    chis: a bunch of rate-scaling parameters
    axis
    
    Returns
    -------
    Ls: dictionary
        Laplacian matrix dictionary for different Laplacian matrices
        (intralayer and interlayer (supra)Laplacians)
    """
    Leval2_eachlayer,Ls = get_Leval2_2layer_eachlayer(G,chis)
    axis.plot(chis,Leval2_eachlayer[0,:].real)
    axis.plot(chis,Leval2_eachlayer[2,:].real,'--')
    axis.set_xlabel('rate-scaling parameter, $\chi$')
    return Ls


######################################################
## supplementary fig 2: bifurcation ##
######################################################
def plot_ReIm_lambda23(G,omegas,deltas):
    """
        Parameters
        ----------
        G: a list of a network's properties
        omegas: a bunch of coupling strengths
        deltas: a bunch of asymmetry parameters
        """
    
    G = get_network_properties(G)
    Levals = get_evals_across_omegas_deltas(G,omegas,deltas)
    lambda2s_thm = get_lambda2_thm(G,deltas)
    
#    S = StrongLimit(G,omegas,deltas)
    fig,ax = plt.subplots(1,2,figsize=(6,3),constrained_layout=True)
    for i_omega in range(len(omegas)):
        ax[0].plot(deltas,Levals[1,i_omega,:].real,'o-',fillstyle='none',label='Exact $Re(\lambda_2)$')
        ax[0].plot(deltas,lambda2s_thm[:,0].real,'k',label='Approx $Re(\lambda_2)$')
#        ax[0].plot(deltas,S['lambda1'][1,:].real,'k',label='Approx $Re(\lambda_2)$')
        ax[0].plot(deltas,Levals[2,i_omega,:].real,'<-',fillstyle='none',label='Exact $Re(\lambda_3)$')
        ax[0].plot(deltas,lambda2s_thm[:,1].real,'--k',label='Approx $Re(\lambda_2)$')
#        ax[0].plot(deltas,S['lambda1'][2,:].real,'--k',label='Approx $Re(\lambda_2)$')

        ax[1].plot(deltas,Levals[1,i_omega,:].imag,'o-',fillstyle='none',label='Exact $Re(\lambda_2)$')
        ax[1].plot(deltas,-abs(lambda2s_thm[:,0].imag),'k',label='Approx $Re(\lambda_2)$')
#        ax[1].plot(deltas,-abs(S['lambda1'][1,:].imag),'k',label='Approx $Re(\lambda_2)$')
        ax[1].plot(deltas,Levals[2,i_omega,:].imag,'<-',fillstyle='none',label='Exact $Re(\lambda_3)$')
        ax[1].plot(deltas,abs(lambda2s_thm[:,1].imag),'--k',label='Approx $Re(\lambda_3)$')
#        ax[1].plot(deltas,abs(S['lambda1'][2,:].imag),'--k',label='Approx $Re(\lambda_3)$')

    ax[0].set_xlabel('asymmetry, $\delta$',fontsize=10)
    ax[1].set_xlabel('asymmetry, $\delta$',fontsize=10)
    ax[0].set_ylabel('$Re(\lambda_2)$ and $Re(\lambda_3)$',fontsize=10)
    ax[1].set_ylabel('$Im(\lambda_2)$ and $Im(\lambda_3)$',fontsize=10)
    ax[0].text(-1.8, 1.2, '(A)',fontsize=14)
    ax[1].text(-1.9, .8, '(B)',fontsize=14)

    lines,labels = ax[0].get_legend_handles_labels()
    labels = ['Exact $\lambda_2$','Approx $\lambda_2$','Exact $\lambda_3$','Approx $\lambda_3$']
    lgd = fig.legend(lines,labels,loc='lower center',ncol=4,bbox_to_anchor=(.5, -.04))
    
    return  fig,ax,lgd
