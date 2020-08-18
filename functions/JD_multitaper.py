from scipy.io.matlab import mio
import numpy as np
import spectrum as sp
from scipy import stats


def pmtmPH(x,dt=1.,nw=3,nfft=None):
    """
    function [P,s,ci] = pmtmPH(x,dt,nw,qplot,nfft);
    Computes the power spectrum using the multi-taper method with adaptive weighting.
    Inputs:
    x      - Input data vector.
    dt     - Sampling interval, default is 1.
    nw     - Time bandwidth product, acceptable values are
    0:.5:length(x)/2-1, default is 3.  2*nw-1 dpss tapers
    are applied except if nw=0 a boxcar window is applied 
    and if nw=.5 (or 1) a single dpss taper is applied.
    qplot  - Generate a plot: 1 for yes, else no.  
    nfft   - Number of frequencies to evaluate P at, default is
    length(x) for the two-sided transform. 
    Outputs:
    P      - Power spectrum computed via the multi-taper method.
    s      - Frequency vector.
    ci     - 95% confidence intervals. Note that both the degrees of freedom
    calculated by pmtm.m and chi2conf.m, which pmtm.m calls, are
    incorrect.  Here a quick approximation method is used to
    determine the chi-squared 95% confidence limits for v degrees
    of freedom.  The degrees of freedom are close to but no larger
    than (2*nw-1)*2; if the degrees of freedom are greater than
    roughly 30, the chi-squared distribution is close to Gaussian.
    The vertical ticks at the top of the plot indicate the size of
    the full band-width.  The distance between ticks would remain
    fixed in a linear plot.  For an accurate spectral estimate,
    the true spectra should not vary abruptly on scales less than
    the full-bandwidth.
    Other toolbox functions called: dpps.m; and if nfft does not equal length(x)    , cz.m
    Peter Huybers
    MIT, 2003
    phuybers@mit.edu

    Adapted from Matlab to Python by Nicolas Barrier"""

    if nfft is None:
        nfft=len(x)

    nx=len(x)
    k=np.min([np.round(2*nw),nx])
    k=np.max([k-1,1])
    s=np.arange(0,1/dt,1/(nfft*dt));
    w=nw/(dt*nx) # half-bandwidth of the dpss
    
    E,V=sp.dpss(nx,NW=nw,k=k)
 
    if nx<=nfft:
        tempx=np.transpose(np.tile(x,(k,1)))
        Pk=np.abs(np.fft.fft(E*tempx,n=nfft,axis=0))**2
    else:
        raise IOError('Not implemented yet')
    
    #Iteration to determine adaptive weights:    
    if k>1:
        xmat=np.mat(x).T
        sig2 = xmat.T*xmat/nx; # power
        P    = (Pk[:,0]+Pk[:,1])/2.;   # initial spectrum estimate
        Ptemp= np.zeros(nfft);
        P1   = np.zeros(nfft);
        tol  = .0005*sig2/nfft;    
        a    = sig2*(1-V);
        while np.sum(np.abs(P-P1)/nfft)>tol:
            Pmat=np.mat(P).T
            Vmat=np.mat(V)
            amat=np.mat(a)
            temp1=np.mat(np.ones((1,k)))
            temp2=np.mat(np.ones((nfft,1)))
            b=(Pmat*temp1)/(Pmat*Vmat+temp2*amat); # weights
            temp3=np.mat(np.ones((nfft,1)))*Vmat
            temp3=np.array(temp3)
            b=np.array(b)
            wk=b**2*temp3       
            P1=np.sum(wk*Pk,axis=1)/np.sum(wk,axis=1)
            Ptemp=P1; P1=P; P=Ptemp;                 # swap P and P1

        #b2=b**2
        #temp1=np.mat(np.ones((nfft,1)))*V
        temp1=b**2
        temp2=np.mat(np.ones((nfft,1)))*Vmat
        num=2*np.sum(temp1*np.array(temp2),axis=1)**2
        
        temp1=b**4
        temp2=np.mat(np.ones((nfft,1)))*np.mat(V**2)
        den=np.sum(temp1*np.array(temp2),axis=1)
        v=num/den
        
    select=np.arange(0,(nfft+1)/2+1).astype(np.int64)
    P=P[select]
    s=s[select]
    v=v[select]

    temp1=1/(1-2/(9*v)-1.96*np.sqrt(2./(9*v)))**3
    temp2=1/(1-2/(9*v)+1.96*np.sqrt(2/(9*v)))**3

    ci=np.array([temp1,temp2])
    
    return P,s,ci

def JD_spectra(ts,dt,ax,f,nrj=0,nw=3,unit='unit',col='k'):

    """
    % calculates multitaper spectra using pmtmPH.m 
    % Author: J. Deshayes, CNRS IRD, March 2013

    Arguments:
    ts=time series whose spectrum to plot. 
    dt=time step of ts in seconds
    ax=the axes in which the drawing has to be done
    f=the frequency where to do draw the errorbar

    Optional arguments:
    nrj=0 for a variance spectrum, 1 for an energy spectrum (F*Px)
    nw=3 (multitaper used 2*nw-1 tapers)
    unit='unit': label of the yaxis
    col='k': color of the line

    Returns:
    Px=the vector of spectrum
    F=the vector of frequency
    PxC=the vector of error bar

    Adapted from Matlab to Python by Nicolas Barrier"""
    
    T=len(ts)
    ts=ts-np.mean(ts)
    [Px,F,Pxc]=pmtmPH(ts,nw=nw);

    F=F*365*86400/dt;    # to get the result in cpy
    Px=Px/(365*86400/dt);    # to get the result in cpy^{-1}
    Pxc=Pxc/(365*86400/dt);    # to get the result in cpy^{-1}
    barmax=Pxc[0,0]*1e2;
    barmin=Pxc[1,0]*1e2;

    if nrj==0:
        hh=ax.loglog(F,Px,color=col);
        ax.set_ylabel('power spectrum density ('+ unit+ '$^2$ cpy$^{-1}$)');
        ax.loglog([f, f],[barmin, barmax],color=col,marker='_');

    else:     
        hh=ax.loglog(F,F*Px,color=col);
        ax.set_ylabel('Energy power spectrum F*Px ('+ unit +'$^2$)');
        
    ax.set_xlabel('frequency (cpy)')
    ax.grid(True,which='both',axis='x')
    ax.grid(True,which='major',axis='y')

    return Px,F,Pxc

def JD_space_spectra(ts,dl,ax,f,nw=3,unit='unit',col='k'):

    """
    % calculates multitaper spectra using pmtmPH.m 
    % Author: J. Deshayes, CNRS IRD, March 2013

    Arguments:
    ts=time series whose spectrum to plot. 
    dl=spatial step of ts in km 
    ax=the axes in which the drawing has to be done
    f=the frequency where to do draw the errorbar

    Optional arguments:
    nw=3 (multitaper used 2*nw-1 tapers)
    unit='unit': label of the yaxis
    col='k': color of the line

    Returns:
    Px=the vector of spectrum
    F=the vector of frequency
    PxC=the vector of error bar

    Adapted from Matlab to Python by Nicolas Barrier"""

    T=len(ts)
    ts=ts-np.mean(ts)
    [Px,F,Pxc]=pmtmPH(ts,nw=nw);

    F=F/dl;    # to get the result in cpy
    Px=Px*dl;    # to get the result in cpy^{-1}
    Pxc=Pxc*dl;    # to get the result in cpy^{-1}
    barmax=Pxc[0,0]*1e2;
    barmin=Pxc[1,0]*1e2;

    hh=ax.loglog(F,Px,color=col);
    ax.set_ylabel('power spectrum density ('+ unit+ ' cpkm$^{-1}$)');
    ax.loglog([f, f],[barmin, barmax],color=col,marker='_');

    ax.set_xlabel('wavenumber (cpkm)')
    ax.grid(True,which='both',axis='x')
    ax.grid(True,which='major',axis='y')

    return Px,F,Pxc,hh

def plot_slope(Px,F,ax,fmin=None,fmax=None,col='k',lin='--',lw=2,offy=1):

    """ This function draws the slope of the spectrum
    Arguments:
    Px=the vector of spectrum
    F=the vector of frequency
    ax=the axes in which the slopes are plotted
    
    Optional arguments:
    fmin=None (if fmin=None, takes the min of F): frequency minimum where the slope is computed
    fmax=None (if fmax=None, takes the max of F): frequency maximum where the slope is computed
    col='k': color of the slope
    lin='--': linestyle of the slope
    lw=2: linewidth of the slope
    offy=1: y offset of the slope
    
    Outputs:
    The value of the slope

    Author: Nicolas Barrier
    """

    if fmin==None:
        fmin=F.min()
    if fmax==None:
        fmax=F.max()

    i=np.nonzero((F>=fmin)&(F<=fmax)&(F!=0))[0]
    fout=F[i]
    pxout=Px[i]

    y=np.log(pxout)
    x=np.log(fout)
    from scipy import stats
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

    confidence_interval = 2.58*std_err #99%
    
    p=np.polyfit(x,y,1)
    
    trend=np.exp(p[0]*x+p[1]+offy)
    
    ax.loglog(fout,trend,color=col,linestyle=lin,lw=lw)

    # add regression slope and confidence interval 


    return p[0]	, slope, confidence_interval
    
def plot_ref_slope(fmin,fmax,f,ax,kvec=[2],col='k',lw=2,ls='--'):

    """ This function draws reference slopes (k=-2, k=-4 for instance)
    Arguments:
    fmin=frequency where to start the reference the slope
    fmax=fmin=frequency where to end the reference the slope
    f= y intercept of the slope
    ax=the axes in which the slopes are plotted
    
    Optional arguments:
    kvec=[2]: list containing the reference values of k, whose slope to draw
    col='k': colors of the slopes
    lw=2: linewidths of the slopes
    ls=--: linestyles of the slopes

    Author: Nicolas Barrier
    """
	
    x=np.linspace(fmin,fmax,5)
    
    for p in range(0,len(kvec)):
        k=kvec[p]
        y=np.log(f)+k*(np.log(fmin)-np.log(x[:]))
        yout=np.exp(y)
        ax.loglog(x,yout,color=col,linewidth=lw,linestyle=ls)
        ax.text(x[-1],yout[-1],' k = -'+str(k),ha='center',va='center',color=col,
        bbox=dict(boxstyle="round", fc="0.9"))
