from netCDF4 import Dataset as nf

def loadvar(filename, varname):
    ff = nf(filename, 'r')
    dat = ff.variables[varname][:]
    ff.close()
    
    return dat

