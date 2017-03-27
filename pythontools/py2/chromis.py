import numpy as np


def dual_fpi(wav, ang = 0.0, erh = 0.0, erl = 0.0):
    
    
    w0 = np.median(wav)
    
    # Reflectivities
    if(w0 < 4010):
        thr = 0.90 + erh
        tlr = 0.80 + erl
    else:
        thr = 0.91 + erh
        tlr = 0.80 + erl
  
  
    # Fix cavity separation
    shr = 357.8e4
    nhr = long(0.5 + shr / (w0 * 0.5))
    hc = nhr * w0 * 0.5
    #
    slr = 136.9e4
    nlr = long(0.5 + slr / (w0 * 0.5)) 
    lc = nlr * w0 * 0.5

    # Finesse
    fhr = 4.0 * thr / (1.0 - thr)**2
    flr = 4.0 * tlr / (1.0 - tlr)**2

    # Phase
    ca = 6.28318530779  * np.cos(ang)
    phr = hc * ca
    plr = lc * ca
  
    # Transmission profiles
    hre = 1.0 / (1.0 + fhr * (np.sin(phr / (wav)))**2)
    lre = 1.0 / (1.0 + flr * (np.sin(plr / (wav)))**2)
  
    return lre * hre

