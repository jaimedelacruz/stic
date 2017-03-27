"""
CRISP class: Generates a theoretical CRISP tranmission profile based on different approximations.
Author: Jaime de la Cruz Rodriguez (USP-SU 2015)
Dependencies: crisp_ref.txt

Usage: 

import crisp as c
fpi = c.crisp(6302)
tw = (np.arange(101)-50)*0.01
tr = fpi.dual_fpi(tw)      # Assumes parallel rays.
tr = fpi.dual_fpi_conv(tw) # includes the convergence of the beam at f165.

"""
import numpy as np
import matplotlib.pyplot as plt
import os

class crisp:
    def __init__(self, w0):
        inam = "crisp::__init__: "
        
        # get reflectivity files
        this_dir, this_filename = os.path.split(__file__)
        self.dfile = os.path.join(this_dir, "crisp_data","crisp_ref.txt")
        if(not os.path.isfile(self.dfile)):
            print inam + "ERROR, reflectivities file not found in "+ self.dfile
            return 0
        self.w0 = float(w0)
        self.read_reflectivity(self.w0)
        

        # Init cavity separations
        shr = 787.e4
        nhr = long(0.5 + shr / (self.w0 * 0.5)) 
        self.hc = nhr * self.w0 * 0.5
        #
        slr = 295.5e4
        nlr = long(0.5 + slr / (self.w0 * 0.5)) 
        self.lc = nlr * self.w0 * 0.5
  
        # other parameters
        self.pi2 = 2.0 * 3.1415926535897932384626433832
        self.calp = np.cos(np.sqrt(np.asarray((0.04691007703067,0.23076534494716, 0.50,0.76923465505284,0.95308992296933), dtype='float64')*(0.5/165.)**2)) * self.pi2
        self.wng = np.asarray((0.11846344252809,0.23931433524968,0.28444444444444,0.23931433524968,0.11846344252809), dtype='float64')
        

        # FSR
        self.hfsr = self.w0**2 / (2.0 * self.hc + self.w0)
        self.lfsr = self.w0**2 / (2.0 * self.lc + self.w0)

        
    def read_reflectivity(self, w0):
        inam = "crisp::read_reflectivity: "
        w, l, h = np.loadtxt(self.dfile, dtype='float64', unpack=True)

        ma = np.max(w)
        mi = np.min(w)
        
        if((w0 > ma) or (w0<mi)):
            print(inam + "Warning, w0={0} is outside the valid range [{1}, {2}], taking value at the edge".format(w0, int(mi), int(ma)))
            w0 = np.max((np.min((w0, ma)), mi))
            
            
        self.lr = np.interp(w0, w, l)
        self.hr = np.interp(w0, w, h)
        print(inam + "(RL,RH)[{0}] = ({1}, {2})".format(w0, self.lr, self.hr))
        


    def dual_fpi(self, wav, ang = 0.0, erh = 0.0, erl = 0.0, ech = 0.0, ecl = 0.0):
        # Total reflectivity
        thr = self.hr  + erh
        tlr = self.lr  + erl


        # Finesse
        fhr = 4.0 * thr / (1.0 - thr)**2
        flr = 4.0 * tlr / (1.0 - tlr)**2


        # Phase
        ca = self.pi2 * np.cos(ang)
        phr = self.hc * ca
        plr = self.lc * ca

  
        # Transmission profiles
        return(1.0 / (1.0 + flr * np.power(np.sin(plr / (wav + ecl + self.w0)),2)) * \
               1.0 / (1.0 + fhr * np.power(np.sin(phr / (wav + ech + self.w0)),2)))


    def getwidth(self):
        n = 2001
        tw = (np.arange(n)-n/2) * 0.001 
        tr = self.dual_fpi(tw)

        dum = np.argmax(tr)
        return( abs(np.interp(0.5, tr[0:dum+1], tw[0:dum+1]) - \
                    np.interp(0.5, tr[-1:dum:-1], tw[-1:dum:-1])) )


    
    def dual_fpi_conv(self, wav, erh = 0.0, erl = 0.0, ech = 0.0, ecl = 0.0):
        # Total reflectivity
        thr = self.hr  + erh
        tlr = self.lr  + erl


        # Finesse
        fhr = 4.0 * thr / (1.0 - thr)**2
        flr = 4.0 * tlr / (1.0 - tlr)**2


        # Phase
        trans = np.zeros(len(wav), dtype='float64')
  
        phr = self.calp * self.hc
        plr = self.calp * self.lc

        wav1 = wav + self.w0

        # Integrate over the beam
        for n in xrange(self.calp.size):
            trans += (1.0 / (1.0 + fhr * np.power(np.sin(phr[n] / (wav1 + ech)),2))) * \
              (1.0 / (1.0 + flr * np.power(np.sin(plr[n] / (wav1 + ecl)),2))) * \
              self.wng[n]
                     
        return(trans)


    
def time2double( v):
    tmp = np.float64(v.split(':'))
    return 3600. * tmp[0] + 60. * tmp[1] + tmp[2]

def double2time( v):
    hh = np.int32(np.fix(v/3600.))
    mm = np.int32(np.fix((v-hh*3600) / 60.))
    ss = v - hh*3600. - mm*60.
    return "%02d"%hh+":%02d"%mm+":%07.4f"%ss
