import numpy as np
import pylab as plt
from scipy import stats

def normalize(x,ddof=1):
    
    """ Returns the normalization of array 1D, ie removes the mean and divide the difference by the standard deviation
    
    Arguments: x (1D) array
    
    Optional arguments:
    ddof=1
    cf http://docs.scipy.org/doc/numpy/reference/generated/numpy.std.html
    
    Author: Nicolas Barrier
    """

    if x.ndim>1:
        raise IOError('x must be 1D (currently, ndim='+str(x.ndim)+')')

    try:
        mmm=x.mask
    except:
        x=np.ma.masked_where(x!=x,x)
        mmm=x.mask
    output=(x-np.mean(x))/(np.std(x,ddof=ddof))
    return output

def bretherton(a,b):
    res=(1-a*b)/(1+a*b);
    return res

def corlag(x,y,maxlag):

    """ Compute the lead-lag correlation between two time series, given a maximum lag.
    
    Arguments:
    x: 1D array (dominates at positive lags)
    y: 1D array
    maxlag: maximum number of lags (integer)
    
    Outputs:
    tuple(ccr,lag), with ccr the cross-correlation and lag the lag vectors
    
    Author: Nicolas Barrier
    """

    x=np.ma.array(x,mask=x!=x)
    
    y=np.ma.array(y,mask=y!=y)

    n=len(x)
    xleady=np.zeros(maxlag+1)
    yleadx=np.zeros(maxlag+1)
   
    #eval(['print'+y+' leads for negative lags'])
    
    lag=np.arange(0,maxlag+1,1)

    for l in lag:
        
        xint=x[l:]
        yint=y[0:n-l]
        iok=np.nonzero((xint.mask==False)&(yint.mask==False))[0]
        if len(iok)>0:
            yleadx[l]=corr_JD.correlation(xint[iok],yint[iok])
  
        yint=y[l:]
        xint=x[0:n-l]
        iok=np.nonzero((xint.mask==False)&(yint.mask==False))[0]
        if len(iok)>0:
            xleady[l]=corr_JD.correlation(xint[iok],yint[iok])[0,1]
  
    ccr=np.empty(2*maxlag+1)
    ccr[0:maxlag]=yleadx[1:][::-1]
    ccr[maxlag:]=xleady
    ccr=np.ma.masked_where(ccr==0,ccr)

    lag=np.arange(-maxlag,maxlag+1,1)
    return ccr,lag

def covlag(x,y,maxlag):

    """ Compute the lead-lag covariance between two time series, given a maximum lag.
    
    Arguments:
    x: 1D array (dominates at positive lags)
    y: 1D array
    maxlag: maximum number of lags (integer)
    
    Outputs:
    tuple(ccv,lag), with ccv the cross-covariance and lag the lag vectors
    
    Author: Julie Deshayes
    """

    x=np.ma.array(x,mask=x!=x)
    
    y=np.ma.array(y,mask=y!=y)

    n=len(x)
    xleady=np.zeros(maxlag+1)
    yleadx=np.zeros(maxlag+1)
    
    lag=np.arange(0,maxlag+1,1)

    for l in lag:
        
        xint=x[l:]
        yint=y[0:n-l]
        iok=np.nonzero((xint.mask==False)&(yint.mask==False))[0]
        if len(iok)>0:
            yleadx[l]=corr_JD.covariance(xint[iok],yint[iok])
  
        yint=y[l:]
        xint=x[0:n-l]
        iok=np.nonzero((xint.mask==False)&(yint.mask==False))[0]
        if len(iok)>0:
            xleady[l]=corr_JD.covarianve(xint[iok],yint[iok])
  
    ccv=np.empty(2*maxlag+1)
    ccv[0:maxlag]=yleadx[1:][::-1]
    ccv[maxlag:]=xleady
    ccv=np.ma.masked_where(ccv==0,ccv)

    lag=np.arange(-maxlag,maxlag+1,1)
    return ccv,lag

def correlation(x,y):

    """ Compute the 0-lag correlation between two time series
    
    Arguments:
    x: 1D array 
    y: 1D array
    
    Outputs:
    corr which is the 0-lat correlation 
    
    Author: Julie Deshayes (CNRS-IRD)
    """

    x=np.ma.array(x,mask=x!=x)
    
    y=np.ma.array(y,mask=y!=y)

    if (x.mask.any()==True):
       raise IOError('there should be no masked value in time series') 
    
    if (y.mask.any()==True):
       raise IOError('there should be no masked value in time series') 
    
    if len(x)==1:
        x=np.squeeze(np.transpose(x))

    if len(y)==1:
        y=np.squeeze(np.transpose(y))

    if len(np.shape(x))>1:
        x=np.squeeze(x)
    
    if len(np.shape(y))>1:
        y=np.squeeze(y)

    n=len(x)
    if len(x)!=len(y):
        raise IOError('x and y must have same length') 

    corr=1./(len(x)-1)*np.sum(normalize(x)*normalize(y))    
    return corr

def covariance(x,y):

    """ Compute the 0-lag covariance between two time series
    
    Arguments:
    x: 1D array 
    y: 1D array
    
    Outputs:
    cov which is the 0-lat covariance
    
    Author: Julie Deshayes (CNRS-IRD)
    """

    x=np.ma.array(x,mask=x!=x)
    
    y=np.ma.array(y,mask=y!=y)

    if (x.mask.any()==True):
       raise IOError('there should be no masked value in time series') 
    
    if (y.mask.any()==True):
       raise IOError('there should be no masked value in time series') 
    
    if len(x)==1:
        x=np.squeeze(np.transpose(x))

    if len(y)==1:
        y=np.squeeze(np.transpose(y))

    if len(np.shape(x))>1:
        x=np.squeeze(x)
    
    if len(np.shape(y))>1:
        y=np.squeeze(y)

    n=len(x)
    if len(x)!=len(y):
        raise IOError('x and y must have same length') 
    
    cova=1./(len(x)-1)*np.sum((x-np.mean(x))*(y-np.mean(y)))
    return cova

def JD_significant(Tm,maxlag,proba=0.95):

    """
    calculates the level of significance for a correlation
    
    Arguments:    
    Tm is the length of the time series
    maxlag: the maximum lag
    
    Optional argument:
    proba=0.95 the significance interval (0.95 -> 95%)
    
    Author: J. Deshayes, CNRS IRD, March 2013
    
    Adapted from Matlab to Python by Nicolas Barrier    
    """

    Nstep=1;
    vectemps=np.arange(1,Tm+1)
    lagf=np.arange(-maxlag,maxlag+Nstep,Nstep)

    rsign=np.zeros(len(lagf))

    for lag in lagf:
        veclag=vectemps+lag;
        I1=np.nonzero(veclag>=1)[0];
        I1=I1[0];
        I2=np.nonzero(veclag<=Tm)[0];
        I2=I2[-1];
        nptcom=I2-I1+1;
        tlim=np.abs(stats.t.ppf(proba, nptcom-2))
        rsign[lag/Nstep+maxlag/Nstep]=tlim/np.sqrt(nptcom-2+tlim*tlim);
    return rsign,lagf

def JD_significant_bretherton(Tm,maxlag,breth,mdf,proba=0.95):

    """
    calculates the level of significance for a correlation following Bretherton (1999)
    
    Arguments:      
    Tm is the length of the time series
    maxlag: the maximum lag
    breth is the coefficient computed using the bretherton(a,b) function
    mdf=the number of degrees of freedom to be removed
        
    Optional arguments:
    proba=0.95 the significance interval (0.95 -> 95%)

    Author: J. Dehsayes, CNRS IRD, March 2013
    
    Adapted from Matlab to Python by Nicolas Barrier"""
 
    Nstep=1;
    vectemps=np.arange(1,Tm+1)
    lagf=np.arange(-maxlag,maxlag+1,Nstep)
    rsign=np.zeros(len(lagf))

    for lag in lagf:
            veclag=vectemps+lag;
            I1=np.nonzero(veclag>=1)[0];
            I1=I1[0];
            I2=np.nonzero(veclag<=Tm)[0];
            I2=I2[-1];
            nptcom=I2-I1+1; 
            dl=(nptcom-mdf)*breth; # effective number of degree of freedom
            tlim=np.abs(stats.t.ppf(proba, nptcom-2))
            rsign[lag/Nstep+maxlag/Nstep]=tlim/np.sqrt(dl-2+tlim*tlim);
    return rsign,lagf


