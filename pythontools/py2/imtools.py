"""
IMTOOLS module contains routines that are usually used to manipulate images.
Author: J. de la Cruz Rodriguez (ISP-SU 2015)
"""
import numpy as np
import ipdb as db
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

def histo_opt(img, cutoff = 1.e-3, top_only = False, bottom_only = False, get_idx = False):
    """
    HISTO_OPT, adapted from P. Suetterlin's IDL histo_opt routine.
    Performs a histogram optimization (not equalization) of a given variable.
    Particularly useful to clip outliers in images.

    INPUT: 
       img = [numpy array]
    
    OPTIONAL KEYWORDS:
       cutoff =      [float 0 < cutoff < 1.0]. The cutoff value, default 1.e-3
       top_only =    [True / False (default)]. Only apply correction to the top values.
       bottom_only = [True / False (default)]. Only apply correction to the bottom values.
       get_idx     = [True / False (default)]. Also return the indexes of the clipped values.

    DEPENDENCIES: Numpy

    AUTHOR: Originally written by P. Sutterlin in IDL, 
            adapted to python by J. de la Cruz Rodriguez (ISP-SU 2015)

    """
    # Compute histogram depending on the data type
    if(img.dtype.name.find('int') == -1): # Image is float/double
        ma = img.max()
        mi = img.min()
        nbins = 10001
        fak = float(nbins-1) / (ma - mi)
        hi, xl = np.histogram(((img-mi)*fak).astype('int32'), bins = nbins)
    else: # Image is integer
        fak = 1
        mi = img.min()
        ma = img.max()
        nbins = ma - mi
        hi, xl = np.histogram(img , bins = nbins)

    nh = hi.size
    hi = np.asarray(hi)

    # Cumulative sum and normalize total power
    for i in range(1,nh,1): hi[i] += hi[i-1]
    hi = np.float32(hi) / hi[-1]
    
    # Locate the value that is below the threshold
    cmin = np.where(hi <= cutoff)[0]
    if(len(cmin) == 0):
        cmin = mi
    else:
        cmin = float(cmin.max()) / fak + mi
        
    cmax = np.where(hi >= (1.0 - cutoff))[0]
    if(len(cmax) == 0):
        cmax = ma
    else:
        cmax = float(cmax.min()) / fak + mi
        
        
    # Return the clipped image
    if(top_only == True and bottom_only == False):
        res =  img.clip(img.min(), cmax)
    elif(bottom_only == True and top_only == False):
        res =  img.clip(cmin, img.max())
    else:
        res =  img.clip(cmin, cmax)

    # Returned indexes of clipped values?
    if(get_idx == True):
        idx = np.where(((img <= cmin) | (img >= cmax)))
        return(res, idx)
    else: return(res)

def f_scale(im, newy, newx):
    """
    f_scale: performs a image scaling based on FFT transform (adding zero-ed high-frequencies).
    Converted from P. Suetterlin's f_scale.

    Input:
         im   = 2D image
         newy = new Y dimension (int)
         newx = new X dimension (int)

    """
    
    sx = im.shape[1]
    sy = im.shape[0]
    sx2 = im.shape[1]/2
    sy2 = im.shape[0]/2
    newimg = np.zeros((newy, newx), dtype=im.dtype) + im.mean()

    sx21 = sx2 + 0
    sy21 = sy2 + 0

    if(sx2*2 != sx): sx2 +=1
    if(sy2*2 != sy): sx2 +=1

    fi = np.fft.fft2(im)
    nfi = np.zeros((newy, newx), dtype='complex128')
    
    if(newx > im.shape[1] and newy > im.shape[0]):
        nfi[0:sy21,0:sx21] = fi[0:sy21, 0:sx21]      # 1st
        nfi[newy-sy2::,0:sx21] = fi[sy-sy2::, 0:sx21]     # 3rd
        nfi[0:sy21,newx-sx2::] = fi[0:sy21,sx-sx2::]    # 2nd
        nfi[newy-sy2::,newx-sx2::] = fi[sy-sy2::,sx-sx2::] # 4th
        
        newimg[:,:] = np.fft.ifft2(nfi * float(newx*newy)/float(sx*sy)).real
        
    return newimg
    


def mirror_image(im):
    sx = im.shape[1]
    sy = im.shape[0]
    nim = np.zeros((2*sy, 2*sx), dtype=im.dtype, order='c')
    nim[0:sy, 0:sx] = im
    nim[0:sy,sx::] = im[:,::-1]
    nim[sy::,:] = nim[0:sy,:][::-1,:]
    return nim

def m_scale(im_in, ynew, xnew):

    im = mirror_image(im_in)
    im1 = f_scale(im, 2*ynew, 2*xnew)[0:ynew, 0:xnew]
    
    
    return(im1)
    
    


def fftconvol2d(im, psf, padding = 1, no_shift = False):
    """
    FFTCONVOL2D: performs the convolution of an image with a PSF, 
                 using FFTs. The rutine arranges the PSF ordering.
    INPUT:
       im  = image (2d numpy array)
       psf = psf (2d numpy array). Should have smaller dims than im (?).

    OPTIONAL KEYWORDS: 
       padding = [int]
          0: Pad extra values with the mean of the image
          1: Pad extra values with the image iteself mirrored (default)
          2: Pad with zeroes
       no_shift = [True/False (default)], don't shift the PSF by half it's size.

    DEPENDENCIES: Numpy

    AUTHOR: J. de la Cruz Rodriguez (ISP-SU 2015)
    """
    # Check dims
    n = np.asarray(im.shape)
    n1 = np.asarray(psf.shape)
    #
    if(len(n) != 2 or len(n1) != 2):
        print("fftconvol2d: ERROR, images must be 2 dimensions:")
        print(" IMG -> ", n)
        print(" PSF -> ", n1)
        return(im)

    # Get padded dims
    npad = n + n1
    off = np.zeros(2, dtype = 'int16')
    #
    for ii in xrange(2):
        if((n1[ii]/2)*2 != n1[ii]):
            npad[ii] -= 1
            off[ii] = 1

    # Create padded arrays
    pim = np.zeros(npad, dtype='float64')
    ppsf= np.zeros(npad, dtype='float64')

    npsf = npad - n
    
    # Copy data to padded arrays
    if(padding == 0): # Pad with the mean of the image
        me = np.mean(im)
        pim[0:n[0], 0:n[1]] = im - me
    elif(padding == 1): # Pad by mirroring a reversed version of the image
        me = 0.0

        pim[0:n[0], 0:n[1]] = im
        
        pim[n[0]:n[0]+npsf[0]/2,0:n[1]] = im[n[0]-1:n[0]-npsf[0]/2-1:-1]
        pim[n[0]+npsf[0]/2:npad[0], 0:n[1]] = im[npsf[0]/2:0:-1,0:n[1]]
        
        pim[:,n[1]:n[1] + npsf[1]/2] = pim[:,n[1]-1:n[1]-npsf[1]/2-1:-1]
        pim[:,n[1]+npsf[1]/2::] = pim[:,npsf[1]/2:0:-1]
        
    elif(padding == 2):
        me = 0.0
        pim[0:n[0], 0:n[1]] = im
    else:
        print("fftconvol2d: ERROR, padding can take values:")
        print("   -> 0: pad the arrays with the average value of the array")
        print("   -> 1: pad the arrays with the same image mirrored (default)")
        print("   -> 2: pad the arrays with zeroes")
        return(0)
    
    # Pad psf and shift
    ppsf[0:n1[0], 0:n1[1]] = psf
    if(not no_shift): 
        ppsf = np.roll(np.roll(ppsf, -n1[0]/2 + off[0] ,axis = 0), -n1[1]/2 + off[1] ,axis = 1)

    
    # Convolve & return
    return(np.fft.irfft2(np.fft.rfft2(pim) * np.fft.rfft2(ppsf))[0:n[0], 0:n[1]] + me)


def rebin(a, new_shape):
    """
    Resizes a 2d array by averaging or repeating elements, 
    new dimensions must be integral factors of original dimensions
    Parameters
    ----------
    a : array_like
        Input array.
    new_shape : tuple of int
        Shape of the output array
    Returns
    -------
    rebinned_array : ndarray
        If the new shape is smaller of the input array, the data are averaged, 
        if the new shape is bigger array elements are repeated
    See Also
    --------
    resize : Return a new array with the specified shape.
    Examples
    --------
    >>> a = np.array([[0, 1], [2, 3]])
    >>> b = rebin(a, (4, 6)) #upsize
    >>> b
    array([[0, 0, 0, 1, 1, 1],
           [0, 0, 0, 1, 1, 1],
           [2, 2, 2, 3, 3, 3],
           [2, 2, 2, 3, 3, 3]])
    >>> c = rebin(b, (2, 3)) #downsize
    >>> c
    array([[ 0. ,  0.5,  1. ],
           [ 2. ,  2.5,  3. ]])
    """
    M, N = a.shape
    m, n = new_shape
    if m<M:
        return a.reshape((m,M/m,n,N/n)).mean(3).mean(1)
    else:
        return np.repeat(np.repeat(a, m/M, axis=0), n/N, axis=1)

def gauss2d(npsf, fwhm):
    sig  = (fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0))))**2
    psf = np.zeros((npsf,npsf), dtype='float64', order='c')
    npsf2 = float(npsf/2)
    for yy in xrange(npsf):
        for xx in xrange(npsf):
            psf[yy,xx] = np.exp(-0.5 * ((xx-npsf2)**2 + (yy-npsf2)**2) / sig)
    return(psf)


def align(a, b):        

    if(a.shape[0] != b.shape[0] or a.shape[1] != b.shape[1]):
        print("align: ERROR, both images must have the same size")
        return(0.0,0.0)
    
    fa = np.fft.fft2(a)
    fb = np.fft.fft2(b)

    cc = np.roll(np.roll(np.real(np.fft.ifft2(fa.conjugate() * fb)), -int(fa.shape[0]/2), axis=0), -int(fa.shape[1]/2), axis=1)
    
    mm = np.argmax(cc)
    xy = ( mm / fa.shape[1], mm % fa.shape[1])

    cc = cc[xy[0]-1:xy[0]+2, xy[1]-1:xy[1]+2]
    y = 2.0*cc[1,1]
    x = (cc[1,0]-cc[1,2])/(cc[1,2]+cc[1,0]-y)*0.5
    y = (cc[0,1]-cc[2,1])/(cc[2,1]+cc[0,1]-y)*0.5

    x += xy[1] - fa.shape[1]/2
    y += xy[0] - fa.shape[0]/2

    return(y,x)
