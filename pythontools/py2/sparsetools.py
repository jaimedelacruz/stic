import numpy as np
from netCDF4 import Dataset as nf
import os
import sys
#import mcmath as mt
import imtools as im
import scipy.ndimage.filters as fil
import scipy.interpolate as inte
import mathtools as mt
#
#import ipdb as db
import matplotlib.pyplot as plt
import sparsetools as sp
from scipy.interpolate import interp2d

#

def nodes(ma, mi, n, ex):
    kk = np.arange(n, dtype='float64')**ex
    kk /= kk.max()
    return(kk*(ma-mi) + mi)
#
def is_odd(num):
    if(num & 0x1): return True
    else: return False

class model:
    def __init__(self,  file=None, nx=0, ny=0, ndep=0, nt=0, verb = True):
        self.dtype = 'float64'
        self.vnames = np.asarray(['ltau500', 'z', 'temp', 'pgas', 'rho', 'vlos', \
                                  'nne', 'pel', 'blong', 'bhor', 'azi', 'vturb', 'B', 'inc'])
        if(file != None):
            self.read(file)
        else:
            self.setSize(nx=nx, ny=ny, ndep=ndep, nt=nt)


    def setSize(self, nx=0, ny=0, ndep=0, nt=1, verb=True):
        if(nx==0 and ny==0 and ndep==0):
            nt = 0
            verb = False
        
        self.nx = nx
        self.ny = ny
        self.ndep = ndep
        self.nt = nt
        
        self.ltau =  np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.temp =  np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.vlos =  np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.vturb = np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.Bln =     np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.Bho =   np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.azi =   np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.pgas =  np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.rho =   np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.nne =   np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.z   =   np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)
        self.cmass  =   np.zeros((nt, ny, nx, ndep), order='c', dtype=self.dtype)

                
        if(verb):
            print("model::setSize: nx = {0}, ny = {1}, ndep = {2}, nt = {3}".\
              format(nx, ny, ndep, nt))


            
    def extract(self, x0=0, x1=-1, y0=0, y1=-1, t0=0, t1=-1, z0=0, z1=-1):
        if(x0 < 0): x0 = 0
        if(x1 > self.nx or x1 < 0): x1 = self.nx
        if(y0 < 0): y0 = 0
        if(y1 > self.ny or y1 < 0): y1 = self.ny
        if(z0 < 0): z0 = 0
        if(z1 >self.ndep or z1 < 0): z1 = self.ndep
        if(t0 < 0): t0 = 0
        if(t1 >self.nt or t1 < 0): t1 = self.nt

        nx = x1 - x0
        ny = y1 - y0
        nt = t1 - t0
        ndep = z1 - z0
        
        print('model::extract: x=({0},{1}), y=({2},{3}), z=({4},{5}), t=({6},{7})'.\
              format(x0,x1,y0,y1,z0,z1,t0,t1))

        m = sp.model(nx = nx, ny=ny, nt=nt, ndep=ndep)
        #
        m.z[:,:,:,:] = self.z[t0:t1,y0:y1,x0:x1,z0:z1]
        m.ltau[:,:,:,:] = self.ltau[t0:t1,y0:y1,x0:x1,z0:z1]
        m.temp[:,:,:,:] = self.temp[t0:t1,y0:y1,x0:x1,z0:z1]
        m.rho[:,:,:,:] = self.rho[t0:t1,y0:y1,x0:x1,z0:z1]
        m.pgas[:,:,:,:] = self.pgas[t0:t1,y0:y1,x0:x1,z0:z1]
        m.nne[:,:,:,:] = self.nne[t0:t1,y0:y1,x0:x1,z0:z1]
        m.vlos[:,:,:,:] = self.vlos[t0:t1,y0:y1,x0:x1,z0:z1]
        m.vturb[:,:,:,:] = self.vturb[t0:t1,y0:y1,x0:x1,z0:z1]
        m.Bln[:,:,:,:] = self.Bln[t0:t1,y0:y1,x0:x1,z0:z1]
        m.Bho[:,:,:,:] = self.Bho[t0:t1,y0:y1,x0:x1,z0:z1]
        m.azi[:,:,:,:] = self.azi[t0:t1,y0:y1,x0:x1,z0:z1]
        m.cmass[:,:,:,:] = self.cmass[t0:t1,y0:y1,x0:x1,z0:z1]
        #
        return(m)

    def _intVar(self, v, nx=1, ny=1, fac = 1):

        if(nx==self.nx and self.ny == ny): return v
            
        facx = np.arange(nx)[0::fac][-1]
        facy = np.arange(ny)[0::fac][-1]
        print("model::_intVar: facy={0}, facx={1}, ny={2}, nx={3}".format(facy, facx, ny-1, nx-1))
        
        xx = np.arange(self.nx, dtype='float32') / (self.nx-1.0)*facx
        yy = np.arange(self.ny, dtype='float32') / (self.ny-1.0)*facy
        
        xx1 = np.arange(nx, dtype='float32')
        yy1 = np.arange(ny, dtype='float32')
        

        m = np.zeros((self.nt, ny, nx, self.ndep), dtype='float32')
        
        for tt in range(self.nt):
            for kk in range(self.ndep):
                inte = interp2d(xx,yy,v[tt,:,:,kk].squeeze(), kind='linear')
                m[tt,:,:,kk] = inte(xx1, yy1)
        return m
    
    def scale(self, ny=1, nx=1, fac=1):

        
        m = model(nx=nx, ny=ny, nt=self.nt, ndep=self.ndep)

        m.ltau[:,:,:,:] = self._intVar(self.ltau, nx=nx, ny=ny, fac=fac)
        m.z[:,:,:,:] = self._intVar(self.z*1.e-5, nx=nx, ny=ny, fac=fac)*1.e5
        m.temp[:,:,:,:] = self._intVar(self.temp, nx=nx, ny=ny, fac=fac)
        m.vlos[:,:,:,:] = self._intVar(self.vlos*1.e-5, nx=nx, ny=ny, fac=fac)*1.e5
        m.vturb[:,:,:,:] = self._intVar(self.vturb*1.e-5, nx=nx, ny=ny, fac=fac)*1.e5
        m.Bln[:,:,:,:] = self._intVar(self.Bln, nx=nx, ny=ny, fac=fac)
        m.Bho[:,:,:,:] = self._intVar(self.Bho, nx=nx, ny=ny, fac=fac)
        m.azi[:,:,:,:] = self._intVar(self.azi, nx=nx, ny=ny, fac=fac)
        m.pgas[:,:,:,:] = self._intVar(self.pgas, nx=nx, ny=ny, fac=fac)
        m.nne[:,:,:,:] = self._intVar(self.nne, nx=nx, ny=ny, fac=fac)
        m.rho[:,:,:,:] = self._intVar(self.rho, nx=nx, ny=ny, fac=fac)
        m.cmass[:,:,:,:] = self._intVar(self.cmass, nx=nx, ny=ny, fac=fac)

        return m

    
    def read(self, file, verb = True, t0 = 0, t1 = -1, x0 = 0, x1 = -1 , y0 = 0, y1 = -1):
        inam =  'model::read: '
        if(not os.path.isfile(file)):
            print(inam+'Error, file {0} not found'.format(file))
            return

        # Check atmos type
        self.type = self._getType(file)
        vtypes = ['Unknown', 'depth-stratified', 'Compressed at nodes']
        if(self.type == 0):
            print(inam + 'ERROR, unknown file type')
            return


        # Check dimensions
        self.nx, self.ny, self.ndep, self.nt = \
          self._getDimensions(file, atype=self.type)

        
        # Printout
        if(verb):
            print(inam+'atmos type -> '+vtypes[self.type])
            print(inam+'nx={0}, ny={1}, ndep={2}, nt={3}'.\
              format(self.nx, self.ny, self.ndep, self.nt))


        # Check sections
        if(t1 == -1): t1 = self.nt
        if(t0 < 0): t0 = 0
        if(t1 > self.nt): t1 = self.nt
        if(x0 < 0): x0 = 0
        if(x1 == -1): x1 = self.nx
        if(x1 > self.nx): x1 = self.nx
        if(y1 == -1): y1 = self.ny
        if(y0 < 0): y0 = 0
        if(y1 > self.ny): y1 = self.ny

            
        # Allocate
        self.setSize(nx=x1-x0, ny=y1-y0, ndep = self.ndep, nt = t1 - t0)
        dims, var = self._getContent(file)

        
        # read vars
        if(self.type == 1):
            if(len(np.where(var == 'temp')[0]) == 1):
                self.temp[:] = self._readVar1(file, "temp", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'vlos')[0]) == 1):
                self.vlos[:] = self._readVar1(file, "vlos", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'vturb')[0]) == 1):
                self.vturb[:] = self._readVar1(file, "vturb", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'blong')[0]) == 1):
                self.Bln[:] = self._readVar1(file, "blong", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'bhor')[0]) == 1):
                self.Bho[:] = self._readVar1(file, "bhor", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'azi')[0]) == 1):
                self.azi[:] = self._readVar1(file, "azi", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'nne')[0]) == 1):
                self.nne[:] = self._readVar1(file, "nne", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'pgas')[0]) == 1):
                self.pgas[:] = self._readVar1(file, "pgas", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'rho')[0]) == 1):
                self.rho[:] = self._readVar1(file, "rho", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'ltau500')[0]) == 1):
                self.ltau[:] = self._readVar1(file,"ltau500",x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'z')[0]) == 1):
                self.z = self._readVar1(file, "z", x0, x1, y0, y1, t0, t1)
            if(len(np.where(var == 'cmass')[0]) == 1):
                self.cmass[:] = self._readVar1(file, "cmass", x0, x1, y0, y1, t0, t1)

            # Backwaards compatibility
            
            if((len(np.where(var == 'b')[0]) == 1) and (len(np.where(var == 'inc')[0]) == 1)):
                tmp = self._readVar1(file, "b", x0, x1, y0, y1, t0, t1)
                tmp1= self._readVar1(file, "inc", x0, x1, y0, y1, t0, t1)
                
                self.Bln[:] = tmp * np.cos(tmp1)
                self.Bho[:] = tmp * np.sin(tmp1)
                
            
    def _readVar1(self, file, vname, x0, x1, y0, y1, t0, t1):
        nc_fid = nf(file, 'r')
        nc_vars = np.asarray([var for var in nc_fid.variables])
        inam = "model::_readVar1: "

        dum = np.where(nc_vars == vname)
        res = None
        if(len(dum[0]) == 1):
            res = nc_fid.variables[vname][t0:t1, y0:y1, x0:x1, :]
        nc_fid.close()
        return(res)
        
    def _getDimensions(self, file, atype=0):
        nc_fid = nf(file, 'r')

        if(atype == 0):
            return(0,0,0,0)
        elif(atype == 1 or atype == 2):
            nx = len(nc_fid.dimensions['x'])
            ny = len(nc_fid.dimensions['y'])
            ndep = len(nc_fid.dimensions['ndep'])
            nt = len(nc_fid.dimensions['time'])
            nc_fid.close()
            return(nx, ny, ndep, nt)
     
    def _getContent(self, file, verb=False):
        inam = "model::_getContent: "
        nc_fid = nf(file, 'r')

        nc_dims = [dim for dim in nc_fid.dimensions]
        nc_vars = np.asarray([var for var in nc_fid.variables])

        if(verb):
            dstr = '' 
            dvar = ' '
            ivar = ''
            ii = 0
    
            for dim in nc_dims:
                dstr += str(dim)+'={0} '.format(len(nc_fid.dimensions[dim]))
            print(inam + 'Dimensions:')
            print('   '+dstr)
            print(' ')
            print(inam + 'Variables:')
            
            for var in nc_vars:
                ivar = ' '
                tmp = nc_fid.variables[var].dimensions
                for d in tmp:
                    ivar += str(d)+' '
                print('   [{0}] '.format(ii)+'{:>10}'.format(var)+' ['+ivar+']')
                ii += 1
            print(" ")

        nc_fid.close()
        return(nc_dims, nc_vars)
    
    def _getType(self, file):
        nc_fid = nf(file, 'r')
        inam = "model::getType: "
        nc_dims = [dim for dim in nc_fid.dimensions]
        nc_vars = np.asarray([var for var in nc_fid.variables])
        nc_fid.close()


        dum = np.where(nc_vars == 'temp')
        if(len(dum[0]) == 1):
            return 1

        dum = np.where(nc_vars == 'model')
        if(len(dum[0]) == 1):
            return 2


        return 0
    def _placenodes(self, tau, nodes):
        nn = len(nodes)
        if(nn == 0): return None

        idx = np.zeros(nn, dtype='int32')
        for ii in range(nn):
            idx[ii] = np.argmin(np.abs(tau-np.asarray(nodes[ii])))
        return(idx)
    
    
    def _equidist(self, tau, n):
        if(n == 1): return(0)

        res = np.zeros(n, dtype='int32', order='c')
        mi = tau.min()
        ma = tau.max()

        res[0] = 0
        res[-1] = tau.size-1

        
        for ii in range(1, n-1):
            mmin = 1.e10
            idx = 0
            xloc = (ma-mi)*ii / (n-1.0) + mi
                    
            for kk in range(1,tau.size-1):
                val = np.abs(tau[kk] - xloc)
                if(mmin > val):
                    mmin = val
                    idx = kk
            res[ii] = idx
        
        return(res)
    
    def _smoothVar(self, nnodes, var, tau, psf, median, t0=0, t1=-1, convolve=False):

        if(t1 == -1): t1 = self.nt

        if(type(nnodes) == list):
            idx = self._placenodes(tau, nnodes)
            numnodes = len(nnodes)
        else:
            idx = self._equidist(tau, nnodes)
            numnodes = nnodes
            
        itau = np.copy(tau[idx].astype(self.dtype), order='c')
        res = var.copy(order='c')
        for tt in range(t0,t1):
            ivar = np.asarray(var[tt,:,:,idx], order='c').reshape((numnodes, self.ny, self.nx))            
               
            
            # Smooth with PSF
            for ii in range(numnodes):
                if(median > 0):
                     ivar[ii,:,:] = fil.median_filter( ivar[ii,:,:].squeeze(), size=median).reshape((ivar.shape[1::]))
                if(len(ivar[ii].squeeze().shape) == 1):
                    idx = np.where(psf == psf.max())
                    ipsf = psf[idx[0], :].squeeze()
                    if(convolve):
                        ivar[ii,:,:] = mt.convolve(ivar[ii].squeeze(), ipsf/ipsf.sum()).reshape(ivar.shape[1::])
                else:
                    if(convolve):
                        ivar[ii,:,:] = im.fftconvol2d(ivar[ii,:,:].squeeze(), psf, padding = 1).reshape(ivar.shape[1::])
            ivar = np.copy(ivar.transpose((1,2,0)), order='c')

            # interpolate depth
            if(numnodes == 1):
                ivar = ivar.squeeze()
                for kk in range(self.ndep):
                    res[tt,:,:,kk] = ivar.reshape(res.shape[1:3])
            elif(numnodes == 2):
                for yy in range(self.ny):
                    for xx in range(self.nx):
                        res[tt,yy,xx,:] = np.interp(tau, itau, ivar[yy,xx,:])
            else:
                for yy in range(self.ny):
                    for xx in range(self.nx):
                        res[tt,yy,xx,:] = np.interp(tau, itau, ivar[yy,xx,:]) #mt.hermpol(itau, ivar[yy,xx,:], tau)
        return(res)
    
    def smooth(self, ntemp=0, nvlos=0, nvturb=0, nBln=0, nBho=0, nazi=0, npgas=0, fwhm = 0.0, t0=0, t1=-1, median = -1):
        
        # Get Gaussian PSF to smooth the vars
        if(fwhm == 0):
            convolve=False
            psf = np.zeros((1,1), dtype='float64')+1.0
        else:
            convolve = True
            npsf = min(int(fwhm*1.5), 3)
            if(not is_odd(npsf)): npsf -= 1
            
            npsf2 = float(npsf/2)
            sig2  = (fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0))))**2
            psf = np.zeros((npsf,npsf), dtype='float64', order='c')
            for yy in range(npsf):
                for xx in range(npsf):
                    psf[yy,xx] = np.exp(-0.5 * ((xx-npsf2)**2 + (yy-npsf2)**2) / sig2)
            psf /= np.sum(psf)
        tau = self.ltau[0,0,0,:]

        # Smooth vars
        if(ntemp > 0):
            self.temp[:] = self._smoothVar(ntemp, self.temp, tau, psf, median, t0=t0, t1=t1, convolve=convolve)
        if(nvlos > 0):
            self.vlos[:] = self._smoothVar(nvlos, self.vlos, tau, psf, median, t0=t0, t1=t1, convolve=convolve)   
        if(nvturb > 0):
            self.vturb[:] = self._smoothVar(nvturb, self.vturb, tau, psf, median, t0=t0, t1=t1, convolve=convolve)
        if(nBln > 0):
            self.Bln[:] = self._smoothVar(nBln, self.Bln, tau, psf, median, t0=t0, t1=t1, convolve=convolve)
        if(nBho > 0):
            self.Bho[:] = self._smoothVar(nBho, self.Bho, tau, psf, median, t0=t0, t1=t1, convolve=convolve)
        if(nazi > 0):
            self.azi[:] = self._smoothVar(nazi, self.azi, tau, psf, median, t0=t0, t1=t1, convolve=convolve)
        if(npgas > 0):
            self.pgas[:] = 10.0**self._smoothVar(npgas, np.log10(self.pgas), tau, psf, median, t0=t0, t1=t1, convolve=convolve)

    def write(self, filename, write_all=True, t0 = 0, \
              t1 = -1, x0=0, x1=-1, y0=0, y1=-1, z0 = 0, z1 = -1):

              
        ncfile = nf(filename,'w', format='NETCDF4', clobber=True)


        # Check ranges
        if(t1 == -1): t1 = self.nt
        if(t1 > self.nt): t1 = self.nt
        if(x1 == -1): x1 = self.nx
        if(y1 == -1): y1 = self.ny
        if(x1 > self.nx): x1 = self.nx
        if(y1 > self.ny): y1 = self.ny
        if(z1 == -1): z1 = self.ndep
        if(z1 > self.ndep): z1 = self.ndep

        nx = x1 - x0
        ny = y1 - y0
        nz = z1 - z0
        nt = t1 - t0
        print("model::write: {4} -> nx={0}, ny={1}, ndep={2}, nt={3}".\
          format(nx, ny, nz, nt, filename))
        sys.stdout.flush()
        
        # Create dimensions
        ncfile.createDimension('y',ny)
        ncfile.createDimension('x',nx)
        ncfile.createDimension('ndep', nz)
        ncfile.createDimension('time')

        # Create vars and write
        par0 = ncfile.createVariable('ltau500','f4',('time','y','x','ndep'))
        par1 = ncfile.createVariable('z','f4',('time','y','x','ndep'))
        par2 = ncfile.createVariable('temp','f4',('time','y','x','ndep'))
        par3 = ncfile.createVariable('pgas','f4',('time','y','x','ndep'))
        par4 = ncfile.createVariable('vlos','f4',('time','y','x','ndep'))
        par5 = ncfile.createVariable('blong','f4',('time','y','x','ndep'))
        par6 = ncfile.createVariable('bhor','f4',('time','y','x','ndep'))
        par7 = ncfile.createVariable('azi','f4',('time','y','x','ndep'))
        par8 = ncfile.createVariable('vturb','f4',('time','y','x','ndep'))
        par9 = ncfile.createVariable('cmass','f4',('time','y','x','ndep'))

        if(write_all):
            par10  = ncfile.createVariable('nne','f4',('time','y','x','ndep'))
            par11 = ncfile.createVariable('rho','f4',('time','y','x','ndep'))


        for tt in range(t0,t1):
            
            par0[tt] = self.ltau[tt,y0:y1, x0:x1, z0:z1]
            par1[tt] = self.z[tt,y0:y1, x0:x1, z0:z1]
            par2[tt] = self.temp[tt,y0:y1, x0:x1, z0:z1]
            par3[tt] = self.pgas[tt,y0:y1, x0:x1, z0:z1]
            par4[tt] = self.vlos[tt,y0:y1, x0:x1, z0:z1]
            par5[tt] = self.Bln[tt,y0:y1, x0:x1, z0:z1]
            par6[tt] = self.Bho[tt,y0:y1, x0:x1, z0:z1]
            par7[tt] = self.azi[tt,y0:y1, x0:x1, z0:z1]
            par8[tt] = self.vturb[tt,y0:y1, x0:x1, z0:z1]
            par9[tt] = self.cmass[tt,y0:y1, x0:x1, z0:z1]
            if(write_all):
                par10[tt] = self.nne[tt,y0:y1, x0:x1, z0:z1]
                par11[tt]= self.rho[tt,y0:y1, x0:x1, z0:z1]

        ncfile.close()


    def resize(self, nx=-1, ny=-1, ndep=-1):
        if(nx < 0): nx = self.nx
        if(ny < 0): ny = self.ny
        if(ndep <0):ndep = self.ndep
        nt = self.nt
            
        nm = sp.model(nx=nx, ny=ny, nt=nt, ndep=ndep)

        xx = np.arange(self.nx, dtype='float64')
        yy = np.arange(self.ny, dtype='float64')
        xx1 = np.arange(nx, dtype='float64') * (self.nx-1) / np.maximum(1, nx-1.0)
        yy1 = np.arange(ny, dtype='float64') * (self.ny-1) / np.maximum(1, ny-1.0)

        db.set_trace()
        for tt in range(nt):
            for kk in range(ndep):
                nm.temp[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.temp[tt,:,:,kk], kind='linear'))(yy1,xx1)
                nm.vlos[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.vlos[tt,:,:,kk], kind='linear'))(yy1,xx1)
                nm.vturb[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.vturb[tt,:,:,kk], kind='linear'))(yy1,xx1)
                nm.Bln[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.Bln[tt,:,:,kk], kind='linear'))(yy1,xx1)
                nm.Bho[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.Bho[tt,:,:,kk], kind='linear'))(yy1,xx1)
                nm.azi[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.azi[tt,:,:,kk], kind='linear'))(yy1,xx1)
                nm.pgas[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.pgas[tt,:,:,kk], kind='linear'))(yy1,xx1)
                nm.rho[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.rho[tt,:,:,kk], kind='linear'))(yy1,xx1)
                nm.nne[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.nne[tt,:,:,kk], kind='linear'))(yy1,xx1) 
                nm.ltau[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.ltau[tt,:,:,kk], kind='linear'))(yy1,xx1)
                nm.z[tt,:,:,kk] = \
                  (inte.interp2d(yy,xx,self.z[tt,:,:,kk], kind='linear'))(yy1,xx1) 
                  
        return(nm)


    def rebin(self, fac=1):

        nx = self.nx/fac
        ny = self.ny/fac

        #if((nx/fac) * fac)

        
        m2 = sp.model(nx = nx, ny = ny, ndep = self.ndep, nt = self.nt)
        for tt in range(self.nt):
            for ii in range(self.ndep):
                m2.temp[tt,:,:,ii] = im.rebin(self.temp[tt,:,:,ii], (ny, nx))
                m2.vlos[tt,:,:,ii] = im.rebin(self.vlos[tt,:,:,ii], (ny, nx))
                m2.vturb[tt,:,:,ii] = im.rebin(self.vturb[tt,:,:,ii], (ny, nx))
                m2.Bln[tt,:,:,ii] = im.rebin(self.Bln[tt,:,:,ii], (ny, nx))
                m2.Bho[tt,:,:,ii] = im.rebin(self.Bho[tt,:,:,ii], (ny, nx))
                m2.azi[tt,:,:,ii] = im.rebin(self.azi[tt,:,:,ii], (ny, nx))
                m2.pgas[tt,:,:,ii] = im.rebin(self.pgas[tt,:,:,ii], (ny, nx))
                m2.rho[tt,:,:,ii] = im.rebin(self.rho[tt,:,:,ii], (ny, nx))
                m2.nne[tt,:,:,ii] = im.rebin(self.nne[tt,:,:,ii], (ny, nx))
                m2.ltau[tt,:,:,ii] = im.rebin(self.ltau[tt,:,:,ii], (ny, nx))
                m2.z[tt,:,:,ii] = im.rebin(self.z[tt,:,:,ii], (ny, nx))



                
        return(m2)



class profile:
    def scale(self, ny=1, nx = 1):
        m = sp.profile(nt=self.nt, ny=ny, nx=nx, nw=self.nw, ns=self.ns)
        
        xx = np.arange(self.nx, dtype='float32')
        yy = np.arange(self.ny, dtype='float32')

        xx1 = np.arange(nx, dtype='float32')/np.maximum(1, nx-1.0) * self.nx
        yy1 = np.arange(ny, dtype='float32')/np.maximum(1, ny-1.0) * self.ny

        me = np.mean(self.dat[:,:,:,:,0])
        if(me < 1.e-19): me = 1.0
        me1 = 1.0 / me
        
        
        for tt in range(self.nt):
            inte = interp2d(xx,yy,self.pweights[tt,:,:], kind='linear')
            m.pweights[yy,:,:] = inte(xx1,yy1)
             
            for ww in range(self.nw):
                for ss in range(self.ns):
                    inte = interp2d(xx,yy,self.dat[tt,:,:,ww,ss].squeeze()*me1, kind='linear')
                    m.dat[tt,:,:,ww,ss] = inte(xx1,yy1)*me

        m.wav[:] = self.wav
        m.weights[:,:] = self.weights[:,:]
        
                     
        return m

    def skip(self, fx=2, fy=2):
        nt,ny,nx,nw,ns = self.dat[:,0::fy, 0::fx, :,:].shape[:]
        m = sp.profile(nx=nx, ny=ny, nt=nt, nw=nw, ns=ns)
        m.dat[:,:,:,:,:] = self.dat[:,0::fy,0::fx,:,:]
        m.weights[:,:] = self.weights[:,:]
        m.wav[:] = self.wav
        m.pweights[:] = self.pweights[:,0::fy,0::fx]
        return m
                    
    def __init__(self, filename = "", nx = 1, ny = 1, nw = 1, ns = 4, nt = 1, \
                 dtype = 'float64', verbose=True):
                 
        self.dtype = dtype
        self.verbose = verbose

        
        if(filename != ""):
            self.read(filename)
        else:
            self.setsize(nx, ny, nw, ns, nt)

            
    def setsize(self, nx = 1, ny = 1, nw = 1, ns = 4, nt = 1):
        
        self.dat = np.zeros((nt, ny, nx, nw, ns), order='c', dtype=self.dtype)
        self.wav = np.zeros((nw), order='c', dtype='float64')
        self.weights = np.ones((nw, ns), order='c', dtype='float64') 
        self.pweights = np.ones((nt, ny, nx), order='c', dtype='float64')

        self.nx = nx
        self.ny = ny
        self.nw = nw
        self.ns = ns
        self.nt = nt
        
        if(self.verbose):
            print( "profile::setsize: nx={0}, ny={1}, nw={2}, ns={3}, nt={4}"\
                   .format(nx, ny, nw, ns, nt))
                   

    def _getContent(self, filename):
        nc_fid = nf(filename, 'r')

        nc_dims = [dim for dim in nc_fid.dimensions]
        nc_vars = np.asarray([var for var in nc_fid.variables])

        nc_fid.close()
        return(nc_dims, nc_vars)
        
                
    def read(self, filename):

        # Get content of file
        nc_dims, nc_vars = self._getContent( filename)
        nc_fid = nf(filename, 'r')

        # Set dims
        nx = 0
        ny = 0
        nw = 0
        ns = 0
        nt = 0
        for ii in nc_dims:
            if(ii == 'x'): nx = len(nc_fid.dimensions['x'])
            elif(ii == 'y'): ny = len(nc_fid.dimensions['y'])
            elif(ii == 'wav'): nw = len(nc_fid.dimensions['wav'])
            elif(ii == 'stokes'): ns = len(nc_fid.dimensions['stokes'])
            elif(ii == 'time'): nt = len(nc_fid.dimensions['time'])
        
        self.setsize(nx, ny, nw, ns, nt)

        read =""
        
        # Readvars
        if(len(np.where(nc_vars == 'profiles')[0]) == 1):
            self.dat[:,:,:,:,:] = nc_fid.variables['profiles'][:]
            read += "[profiles]"
        if(len(np.where(nc_vars == 'wav')[0]) == 1):
            self.wav[:] = nc_fid.variables['wav'][:]
            read += "[wav]"
        if(len(np.where(nc_vars == 'weights')[0]) == 1):
            self.weights[:,:] = nc_fid.variables['weights'][:]
            read+= "[weights]"
        if(len(np.where(nc_vars == 'pixel_weights')[0]) == 1):
            self.pweights[:,:,:] = nc_fid.variables['pixel_weights'][:]
            read+="[pixel_weights]"
        print("profile::read: "+read)
        nc_fid.close()
        
    def write(self, filename):
        nc_fid = nf(filename, 'w')

        # Create dims
        nc_fid.createDimension('x',self.nx)
        nc_fid.createDimension('y',self.ny)
        nc_fid.createDimension('stokes',self.ns)
        nc_fid.createDimension('wav',self.nw)
        nc_fid.createDimension('time')
        
        if(self.dat.dtype == 'float64'): prec = 'f8'
        else: prec = 'f4'
        
        prof = nc_fid.createVariable('profiles', prec ,('time','y', 'x', 'wav', 'stokes'))
        wav = nc_fid.createVariable('wav', 'f8' , ('wav'))
        wei = nc_fid.createVariable('weights', 'f4' , ('wav','stokes'))
        pwe = nc_fid.createVariable('pixel_weights', 'f4' , ('time','y','x'))
        print("profile::write: saving data [{0}]".format(filename))

        # Write data
        for ii in range(self.nt):
            prof[ii,:,:,:,:] = self.dat[ii,:,:,:,:]
            pwe[ii] = self.pweights[ii,:,:]
        
        wav[:] = self.wav
        wei[:] = self.weights

        nc_fid.close()

    def __add__(self, o):
        
        if(self.nt == o.nt and\
           self.nx == o.nx and\
           self.ny == o.ny and\
           self.ns == o.ns):
           nw = o.nw + self.nw
        else: return None

            
        n = sp.profile(nx=self.nx, ny=self.ny, nw = nw, ns = self.ns, nt = self.nt)

        for ss in range(self.ns):
            n.weights[:,ss] = np.append(self.weights[:,ss].squeeze(), o.weights[:,ss].squeeze()) 
            for tt in range(self.nt):
                for yy in range(self.ny):
                    for xx in range(self.nx):
                        n.dat[tt,yy,xx,:,ss] = np.append(self.dat[tt,yy,xx,:,ss].squeeze().copy(), o.dat[tt,yy,xx,:,ss].squeeze().copy())

        n.wav[:] = np.append(self.wav, o.wav)
        n.pweights[:,:,:] = 0.5 * (self.pweights + o.pweights)

        return(n)

    def averageSpectrum(self):
        return(self.dat.sum(axis=(0,1,2))/(self.nx*self.ny*self.nt))
        
    
    def extractWav(self, w0 = 0, w1 = -1):
        if(w1 == -1): w1 = self.nw
            
        nw = w1 - w0
        
        n = sp.profile(nx = self.nx,\
                       ny = self.ny,\
                       nw = nw,\
                       ns = self.ns,\
                       nt = self.nt)
        
        n.dat[:,:,:,:,:] = self.dat[:,:,:,w0:w1,:].copy()
        n.wav[:] = self.wav[w0:w1].copy()
        n.weights[:,:] = self.weights[w0:w1,:].copy()
        n.pweights[:,:,:] = self.pweights.copy()

        return(n)

    def extractPix(self, x0 = 0, x1 = -1, y0 = 0, y1 = -1, t0=0, t1=-1):
        if(x1 == -1): x1 = self.nx
        if(y1 == -1): y1 = self.ny
        if(t1 == -1): t1 = self.nt

        nx = x1 - x0
        ny = y1 - y0
        nt = t1 - t0

        n = sp.profile(nx = nx, ny = ny, nt = nt, nw = self.nw)

        n.dat[:,:,:,:,:] = self.dat[t0:t1,y0:y1, x0:x1, :, :]
        n.wav[:] = self.wav
        n.weights[:,:] = self.weights
        n.pweights[:,:,:] = self.pweights[t0:t1, y0:y1, x0:x1]

        return(n)
    
    def splitRegions(self):
        nReg = 0
        w = self.wav
        dw = w[1]-w[0]
        w0 = 0
        m1 = []
        odw = dw
        nwav = w.size-1
        ii = 1
        while(ii < nwav):
            odw = dw
            dw = w[ii]-w[ii-1]
            if(np.abs(odw-dw) > 1.e-4):
                ddw = w[ii+1]-w[ii]
                m1.append(self.extractWav(w0=w0,w1=ii))
                w0=ii
                dw = ddw

                if(ddw > 1):
                    m1.append(self.extractWav(w0=ii,w1=ii+1))
                    ii+=1
                    w0=ii
                    dw = w[ii+1]-w[ii]
                ii+=1
            
            else:
                ii+=1
            
        m1.append(self.extractWav(w0=w0,w1=w.size))
        return(m1)
    
def writeInstProf(oname, var, pref=[]):
    ncfile1 = nf(oname,'w', format='NETCDF4')
    ncfile1.createDimension('wav',var.size)
    par1 = ncfile1.createVariable('iprof','f8',('wav'))
    par1[:] = var


    if(len(pref) == 3):
        ncfile1.createDimension('np',len(pref))
        par2 = ncfile1.createVariable('pref','f8',('np'))
        par2[:] = np.float32(pref)

    ncfile1.close()
    
def readVALD(filename, writeto='vald_lines.cfg', stellar=True, width = 3.0, verbose=False):
    f = open(filename, 'r')
    li = []

    nreg = 4
    if(stellar):
        nli = int(f.readline().split(',')[2])
        print(nli)
        
        for ii in range(0,(nli)*nreg+2):
            line = f.readline()
            if(ii < 2): continue
            li.append(line)

    else:
        for line in f:
            if(line[0] == "'"):
                li.append(line)

    nli = len(li)
    lines = np.asarray(li).reshape((nli/nreg, nreg))
    #Elm Ion      WL_air(A)   log gf* E_low(eV) J lo  E_up(eV) J up  lower  upper   mean   Rad.  Stark  Waals
    ## Label    	Elem	ion	anum	wav		loggf	j_low	j_up	g_low	g_up	e_low	gam_rad		gam_strk	vdW/Bark	wid

    id_el=np.asarray(['H ','He','Li','Be','B ','C ','N ','O ','F ',
                      'Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ',
                      'Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu',
                      'Zn','Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ',
                      'Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In',
                      'Sn','Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr',#$  ; 5x
                      'Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm',#$  ; 6x
                      'Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au',#$  ; 7x
                      'Hg','Tl','Pb','Bi','Po','At','Rn','Fa','Ra','Ac',#$  ; 8x
                      'Th','Pa','U '])
    
    fo = open(writeto,'w')
    line = "{:<20}{:>9}{:>7}{:>13}{:>11}{:>7}{:>7}{:>9}{:>9}{:>10}{:>10}{:>10}{:>10}{:>9}\n".\
           format('# Label', 'Elem+Ion', 'A_num','wl_air', 'loggf','j_l','j_u', 'g_l', 'g_u', 'e_low', 'Grad', 'Gstark', 'GvdW', 'Width')
    fo.write(line)
    
    
    for ii in range(nli/nreg):
        l = lines[ii,:]
        l0 = l[0].split(',')
        anum = np.where(id_el == "{:<2}".format(l0[0][1:-1].split()[0]))[0] + 1
        label = "_".join(l0[0][1:-1].split())+'_'+l0[1].strip()
        if(len(anum) == 0): continue
        
        line = "{:<20}{:>9}{:>7}{:>13}{:>11}{:>7}{:>7}{:>9}{:>9}{:>10}{:>10}{:>10}{:>10}{:>9}\n".\
               format(label,l0[0][1:-1], anum[0], l0[1].strip(), l0[2].strip(), l0[4].strip(),
                      l0[6].strip(), l0[7].strip(), l0[8].strip(), l0[3].strip(),\
                      l0[10].strip(),  l0[11].strip(), l0[12].strip(), width)
        
        
        fo.write(line)
        if(verbose): print(line)
    fo.close()
    
    return lines
