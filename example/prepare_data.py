"""
STIC inversion example

In this example, we prepare a dataset consisting of Ca II K, Ca II 8543 and Fe I 6301/6302
that was recorded with CRISP and CHROMIS.

The assumption is that the data has been calibrated to absolute CGS intensity units.
All dependencies should be available in stic/pythontools/py2/

Modifications:
                  2017-03-27, JdlCR: Created!

"""
import numpy as np
import matplotlib.pyplot as plt
import sparsetools as sp
import mfits as mf
import prepare_data as S
import satlas
import crisp as fpi
import chromis as cr
from netCDF4 import Dataset as nc

#-----------------------------------------------------------------------------------------
# SOME DEFINITIONS
#-----------------------------------------------------------------------------------------

def findgrid(w, dw, extra=8):
    w1=np.round(w*1000).astype('int32')
    dw1 = int(dw*1000)

    w2 = w/dw
    w2 = np.round(w2)
    
    idx = np.zeros(w.size, dtype='int32')

    np0 = w2[-1] - w2[0] + extra
    wn = np.arange(np0, dtype='float64')*dw + w[0] -extra/2*dw

    for ii in range(w.size):
        idx[ii] = np.argmin(np.abs(wn-w[ii]))
    
    return wn, idx

#-----------------------------------------------------------------------------------------

def getCont(lam):
    s = satlas.satlas()
    x,y,c = s.getatlas(lam-0.1,lam+0.1, cgs=True)
    return np.median(c)

#-----------------------------------------------------------------------------------------

def writeInstProf(oname, var, pref=None):
    ncfile1 = nc(oname,'w', format='NETCDF4')
    ncfile1.createDimension('wav',var.size)
    par1 = ncfile1.createVariable('iprof','f8',('wav'))
    par1[:] = var


    if(len(pref) == 3):
        ncfile1.createDimension('np',len(pref))
        par2 = ncfile1.createVariable('pref','f8',('np'))
        par2[:] = np.float32(pref)

    ncfile1.close()
    
#-----------------------------------------------------------------------------------------

#
# MAIN PROGRAM
#
if __name__ == "__main__":

    #read data ([wav, StokesI, StokesQ, StokesU, StokesV])
    
    fe  = mf.readfits('feobs.fits')
    ca8= mf.readfits('ca8obs.fits')
    ck  = mf.readfits('cakobs.fits')


    # Get the continuum level for all lines using the FTS atlas

    cw = np.asarray([4000., 8542., 6302.])
    cont = []
    for ii in cw: cont.append(getCont(ii))
    
    # These observations are not acquired in a regular wavelength grid
    # Find a finer grid that contains all wavelength points for each line
    # This grid should be at least the FWHM of the instrumental profile

    
    # 6301 and 6302, each in a different region. The CRISP FWHM is
    # ~30 mA at 6302/
    # We try to find a grid of ~10 mA, that should be able to fit all points
    # The original grid is multiple of ~40 mA. wfe1 is the new grid and ife1
    # is the location of the observed points in the new grid
    
    wfe1, ife1 = findgrid(fe[0,0:9], (fe[0,4]-fe[0,3])*0.25, extra=8)
    wfe2, ife2 = findgrid(fe[0,9::],  (fe[0,4]-fe[0,3])*0.25, extra=12)
    wfe = np.append(wfe1,wfe2)
    ife   = np.append(ife1, ife2+wfe1.size)
    
    # 8542, the observations where performed in a grid multiple of ~85 mA.
    # The FWHM of CRISP at 8542 is ~55 mA so 1/2 of the original grid should do.

    wc8, ic8 = findgrid(ca8[0,:], (ca8[0,10]-ca8[0,9])*0.5, extra=8)
    
    # Ca II K, the observations are recorded in a grid of 3 km/s + 2 external points
    # These profiles are in theory critically sampled, but we need to add extra points
    # for the outer points and for the continuum (last point in the array).

    wck, ick = findgrid(ck[0,0:39], (ck[0,10]-ck[0,9]), extra=8)
    
    
    #
    # Now we create a container for each spectral region in STIC format
    # We will add the fine grid, but all points that were not observed will
    # be given weight zero, so they don't contribute to the inversion.
    # 

    fe_1 = sp.profile(nx=1, ny=1, ns=4, nw=wfe.size)
    ca_8 = sp.profile(nx=1, ny=1, ns=4, nw=wc8.size)
    ca_k = sp.profile(nx=1, ny=1, ns=4, nw=wck.size+1) # add 1 continuum point!


    # Now fill in the profiles, weights and wavelength

    fe_1.wav[:] = wfe[:]
    ca_8.wav[:] = wc8[:]
    ca_k.wav[0:-1] = wck[:]
    ca_k.wav[-1]    = ck[0,-1] # Continuum point

    # Fill arrays with the observed points. The rest can be left zeroed.
    # To scale the problem to numbers close to 1, the code allows to
    # choose a normalization factor. I normally choose the quiet-Sun
    # continuum. It does not really matter as long as we tell the code
    # about this normalization in the input file. The normalization can
    # be different for each region.
    
    for ii in range(4): # Fill all stokes parameters for 8542 and 6302
        fe_1.dat[0,0,0,ife,ii] = fe[ii+1,:] /cont[2] 
        ca_8.dat[0,0,0,ic8,ii] = ca8[ii+1,:] / cont[1]

        
    # Only fill Stokes I for Ca II K
    ca_k.dat[0,0,0,ick,0] = ck[1,0:39] / cont[0]
    ca_k.dat[0,0,0,-1,0] = ck[1,39]  / cont[0]


    # Now fill-in weights for the inversion. In practice this should be the noise level
    # But we also use it to mask all the points that were not observed.
    fe_1.weights[:,:] = 1.e16 # Very high value means weight zero
    fe_1.weights[ife,:] = 0.005
    fe_1.weights[ife,1:3] /= 4.5 # Some more weight for Q&U
    fe_1.weights[ife,3] /= 3.5    # Some more weight for V

    ca_8.weights[:,:] = 1.e16 # Very high value means weight zero
    ca_8.weights[ic8,:] = 0.004
    ca_8.weights[ic8,1:3] /= 7.0 # Some more weight for Q&U
    ca_8.weights[ic8,3] /= 4.0    # Some more weight for V
    ca_8.weights[ic8[9:12],0] /= 2.0
    
    ca_k.weights[:,:] = 1.e16 # Very high value means weight zero
    ca_k.weights[ick,0] = 0.002
    ca_8.weights[ick[19:22],0] /= 2.0
    ca_k.weights[-1,0] = 0.004 # Continuum point
    

    # Now combine all regions in one object.
    # The "sum" operator is overloaded to merge regions
    # Write the profiles and weights to HD
    
    sp_all = ca_k + ca_8 + fe_1
    sp_all.write('observed.nc')
    

    #
    # Now print the regions for the config file
    #
    lab = "region = {0:10.5f}, {1:8.5f}, {2:3d}, {3:e}, {4}"
    print(" ")
    print("Regions information for the input file:" )
    print(lab.format(ca_k.wav[0], ca_k.wav[1]-ca_k.wav[0], ca_k.wav.size-1, cont[0], 'fpi, 3934.nc'))
    print(lab.format(ca_k.wav[-1], ca_k.wav[1]-ca_k.wav[0], 1, cont[0], 'none, none'))
    print(lab.format(ca_8.wav[0], ca_8.wav[1]-ca_8.wav[0], ca_8.wav.size, cont[1], 'fpi, 8542.nc'))
    print(lab.format(wfe1[0], wfe1[1]-wfe1[0], wfe1.size, cont[2], 'fpi, 6302.nc'))
    print(lab.format(wfe2[0], wfe2[1]-wfe2[0], wfe2.size, cont[2], 'fpi, 6302.nc'))
    print("(w0, dw, nw, normalization, degradation_type, instrumental_profile file)")
    print(" ")

    #
    # Generate Instrumental profiles and save them
    #

    # Ca K region
    dw =  ca_k.wav[1]-ca_k.wav[0]
    ntw= 25 # Always an odd number < ca_k.wav.size
    tw1 = (np.arange(ntw)-ntw/2)*dw + 3934.0
    tr1 = cr.dual_fpi(tw1, erh = -0.1)
    tr1 /= tr1.sum()
    # Stores the FPI profile and the parameters of the prefilter
    writeInstProf('3934.nc', tr1, [ca_k.wav[ick[39/2]], 4.5, 3.0]) 

    
    # Ca II 8542
    dw =  ca_8.wav[1]-ca_8.wav[0]
    ntw= 25
    f=fpi.crisp(8542.0)
    tw = (np.arange(ntw)-ntw/2)*dw 
    tr = f.dual_fpi(tw, erh = -0.02)
    tr /= tr.sum()
    writeInstProf('8542.nc', tr,  [8542.091, 9.0, 2.0])

    
    # 6301/6302, we will use the same profile, so it should not have more points than any
    # of the 2 regions
    dw =  fe_1.wav[1]-fe_1.wav[0]
    ntw= 45 
    f=fpi.crisp(6302.0)
    tw = (np.arange(ntw)-ntw/2)*dw 
    tr = f.dual_fpi(tw, erh=-0.01)
    tr /= tr.sum()
    writeInstProf('6302.nc', tr,  [6302.1, 4.4, 2.0])
    
    
    #
    # Init the input model, all quantities in CGS units!
    #

    # First create a tau scale
    
    taumin = -8.6
    taumax= 0.8
    dtau = 0.18
    ntau = int((taumax-taumin)/dtau) + 1
    tau = np.arange(ntau, dtype='float64')/(ntau-1.0) * (taumax-taumin) + taumin

    # Now create a smooth temperature profile
    temp = np.interp(tau, np.asarray([-8.6, -6.0, -4.0, -2.0 , 0.8]), np.asarray([20000., 8000., 4000., 4800., 7000.]))
    
    # Fill in the model
    m = sp.model(nx=1, ny=1, nt=1, ndep=ntau)
    m.ltau[0,0,0,:] = tau
    m.temp[0,0,0,:] = temp

     # The inversion only needs to know the gas pressure at the upper boundary. FALC has Pgas[top] ~ 0.3, but
     # this value is for quiet-Sun. Active regions can have up to Pgas[top] = 10.
     
    m.pgas[0,0,0,:] = 3.0    

    # Fill in initial B field and velovity (optional)
    m.vturb[0,0,0,:] = 1.e5
    m.vlos[0,0,0,:] = 0.5e5 # cm/s
    m.B[0,0,0,:] = 900.
    m.inc[0,0,0,:] = 80. * 3.14159 / 180.
    m.azi[0,0,0,:] = 100. * 3.14159 / 180.

    # Write to HD
    m.write('modelin.nc', write_all=True)
    
