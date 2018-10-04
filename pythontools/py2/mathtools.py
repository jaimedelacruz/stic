"""
Mathtools: Bezier intepolation, 1D FFT-based convolution, parabolic fit to 3 points, etc...
Author: Jaime de la Cruz Rodriguez (ISP-SU 2015)
"""
import numpy as np

def convolve(var, tr):
    # Dimensions
    n = len(var)
    n1 = len(tr)
    npad = n + n1
    
    if((n1//2)*2 != n1):
        npad -= 1
        off = 1
    else:
        off = 0
    
    # Pad arrays using wrap around effect
    pvar = np.empty(npad, dtype='float64')
    pvar[0:n] = var
    pvar[n:n+n1//2] = var[-1]
    pvar[n+n1//2::] = var[0]
    
    ptr = np.zeros(npad, dtype = 'float64')
    ptr[0:n1] = tr / np.sum(tr)
    ptr = np.roll(ptr, -n1//2 + off)
    
    # FFT, convolve and FFT back
    return( (np.fft.irfft(np.fft.rfft(pvar) * np.fft.rfft(ptr)))[0:n])


def cent_der(x, y):
    n = len(y)
    yp = np.zeros(n, dtype = 'float64')

    dx0 = x[1::] - x[0:-1]
    der0 = (y[1:-1] - y[0:-2]) / dx0[0:-1]

    dx1 = x[2::] - x[1:-1]
    der1 = (y[2::] - y[1:-1]) / dx1
    
    idx = (np.where(der0*der1 > 0.0))[0]
    if(len(idx) > 0):
        lamb = (1.0 + dx1[idx] / (dx1[idx] + dx0[idx])) / 3.0
        yp[idx+1] = (der0[idx] * der1[idx] ) / (lamb * der1[idx] + (1.0 - lamb) * der0[idx])

    yp[0] = der0[0]
    yp[-1] = der1[-1]

    return(yp, dx0)


def bezier2(x, y, xp, linear = False):
    # Centered derivatives
    n = len(x)
    n1 = len(xp)
    yp, dx = cent_der(x,y)

    # compute interpolated values
    res = np.zeros(n1)
    
    for k in range(n-1):
        k1 = k+1
        idx = np.where((xp >= x[k]) & (xp < x[k1]))[0]
        if(len(idx) == 0):
            continue

        # Control point
        cntr = 0.5 * (y[k] + dx[k] * 0.5 * yp[k] + y[k1] - dx[k] * 0.5 * yp[k1])

        # Result
        u = (xp[idx] - x[k]) / dx[k]
        u1 = 1.0 - u
        res[idx] = y[k] * u1**2 + y[k1] * u**2 + 2.0 * cntr * u * u1

    # Values outside the x, y domain -> Extrapolate or repear?
    ma = np.max(x)
    mi = np.min(x)
    idx0 = np.where(xp <  mi)[0]
    idx1 = np.where(xp >= ma)[0]

    if(len(idx0)):
        if(linear):
            b = y[0] - yp[0]*x[0]
            res[idx0] = yp[0]*xp[idx0] + b
        else:
            res[idx0] = y[0]

    if(len(idx1)):
        if(linear):
            b = y[-1] - yp[-1]*x[-1]
            res[idx1] = yp[-1]*xp[idx1] + b
        else:
            res[idx1] = y[-1]
          
    return(res)


def bezier3(x, y, xp, linear = False):
    # Centered derivatives
    n = len(x)
    n1 = len(xp)
    yp, dx = cent_der(x,y)

    # compute interpolated values
    res = np.zeros(n1)
    
    for k in range(n-1):
        k1 = k+1
        idx = np.where((xp >= x[k]) & (xp < x[k1]))[0]
        if(len(idx) == 0):
            continue
        
        # Control points
        kdx = dx[k] / 3.0
        cntr  =  y[k]  + kdx * yp[k] 
        cntr1 =  y[k1] - kdx * yp[k1]

        # Result
        u = (xp[idx] - x[k]) / dx[k]
        u1 = 1.0 - u
        res[idx] = y[k] * u1**3 + y[k1] * u**3 + 3.0 * cntr * u * u1**2 + 3.0 * cntr1 * u**2 * u1

    # Values outside the x, y domain -> Extrapolate or repear?
    ma = np.max(x)
    mi = np.min(x)
    idx0 = np.where(xp <  mi)[0]
    idx1 = np.where(xp >= ma)[0]

    if(len(idx0)):
        if(linear):
            b = y[0] - yp[0]*x[0]
            res[idx0] = yp[0]*xp[idx0] + b
        else:
            res[idx0] = y[0]

    if(len(idx1)):
        if(linear):
            b = y[-1] - yp[-1]*x[-1]
            res[idx1] = yp[-1]*xp[idx1] + b
        else:
            res[idx1] = y[-1]
          
    return(res)

def parab_fit(x, y):
    cf = np.empty(3, dtype='float64')

    d = x[0]
    e = x[1]
    f = x[2]
  
    yd = y[0]
    ye = y[1]
    yf = y[2]
  
    cf[1] = ((yf - yd) - (f**2 - d**2) * ((ye - yd) / (e**2 - d**2)))/ \
      ((f - d) - (f**2 - d**2) * ((e - d) / (e**2 - d**2)))
    #
    cf[2] = ((ye - yd) - cf[1] * (e - d)) / (e**2 - d**2)
    #
    cf[0] = yd - cf[1] * d - cf[2] * d**2
    #
    return( cf )

def gaussian(x, p):
    sig = p[2] / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    z = (x - p[1]) / sig
    return(p[0] * np.exp(-0.5 * z**2))
