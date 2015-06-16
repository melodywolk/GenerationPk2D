import numpy as np

def lnrhostar(Nca, n12):
    """
    returns the value of the density rhobar as function of s2a and nbar
    s2a is the variance of the log, nbar the mean density in a cell
    n12 is (N-1/2)/Nbar
    Nca is s2a*Nbar
    Uses newtonraphson
    """

    # starting point
    y0 = np.log(np.abs(n12))
#    print 'yo %f' % y0
    #if N ne 0 then y0 =alog(N/Nbar) else y0 = 1.d0
    err = np.exp(y0)+ y0/Nca - n12
#    print 'err %f' % err
    I = 0
    STEPMAX = 1000
    TOL = 10**(-12)
    while ((np.abs(err) > TOL) and (I < STEPMAX)):
         y1 = y0 - (err)/(np.exp(y0) + 1.0/Nca)
         err = np.exp(y1)+ y1/Nca - n12
         y0 = y1
         I += 1

    if I > STEPMAX:
       print 'NewtonRaphson failed'
    else:   
       return y0
    

