import numpy as np
import os
import sys
import ipdb as ip

def lphead(name, verbose = True):
    inam = 'lphead:'
            
    # Open file
    datfil = open(name, 'rb')

    # get header and extract dimensions
    head = (np.fromfile(datfil, dtype=np.dtype('a512'), count=1))[0]
    dum = head.split()
    datfil.close()
    
    ndims = 0
    dtype = 0
    nx = 0
    ny = 0
    nt = 0
    nstokes = 1
    for it in dum:
        du1 = it.split("=")
        if("dims" in du1[0]):
            ndims = int((du1[1].split(','))[0])
        elif("datatype" in du1[0]):
            dtype = int(du1[1].split(',')[0])
        elif("nx" in du1[0]):
            nx = int(du1[1].split(',')[0])
        elif("ny" in du1[0]):
            ny = int(du1[1].split(',')[0])
        elif("nt" in du1[0] and len(du1[0]) == 2):
            nt = int(du1[1].split(',')[0])
        elif("stokes" in du1[0]):
            du2 = du1[1].split(']')
            nstokes = int(np.size(du2[0].split(',')))
            

    if(dtype == 1):
        dtype = np.dtype('b')
    elif(dtype == 2):
        dtype = np.dtype('h')
    elif(dtype == 3):
        dtype = np.dtype('i')
    elif(dtype == 4):
        dtype = np.dtype('f')
    elif(dtype == 5):
        dtype = np.dtype('d')
    else:
        print inam, 'Warning, dtype={0} not supported!'.format(dtype)

    if(verbose):
        print inam, "[dtype={0}, ndims={1}, nx={2}, ny={3}, nt={4}, nstokes={5}] -> {6}".format(dtype, ndims, nx, ny, nt, nstokes, os.path.basename(name))
    return nx, ny, nt, nstokes, dtype, ndims  


def mk_lpheader(nx, ny, nt, dtype = 'float32', stokes = False):
    if(dtype == 'float32'):
        idt = 4
        stringtype = '(float)'
    elif(dtype == 'int8'):
        idt = 1
        stringtype = '(byte)'
    elif(dtype == 'int16'):
        idt = 2
        stringtype = '(integer)'
    elif(dtype == 'int32'):
        idt = 3
        stringtype = '(long)'
    elif(dtype == 'float64'):
        idx = 5
        stringtype = '(double)'
    else:
        print "mk_pol_lpheader: ERROR, supported types are: int8, int16, int32, float32, float64"
        return 0

    # Full-Stokes ?
    if(stokes): header = 'stokes=[I,Q,U,V], ns=4 :  '
    else: header = ''

    # Common part
    header += 'datatype={0}'.format(idt)+' '+stringtype
    header += ', dims={0}'.format(3)
    header += ', nx={0}'.format(nx)
    header += ', ny={0}'.format(ny)
    header += ', nt={0}'.format(nt)
    if(sys.byteorder == 'little'): endian = 'l'
    else: endian = 'b'
    header += ', endian='+endian

    # Encode to byte array
    head = np.zeros(512, dtype='int8')
    head[0:len(header)] = map(ord, header.encode('utf-8'))
    
    return(head)

