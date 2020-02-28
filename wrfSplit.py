#!/usr/bin/env python3
import argparse, configparser, os, sys
from netCDF4 import Dataset
import numpy as np 
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'Split WRF into n files by latitude chunks')
    parser.add_argument('--in', help = 'path to ncdiag input', required = True, dest = 'input')
    parser.add_argument('--output-dir', help = 'path to ncdiag output', required = True, dest = 'outdir')
    parser.add_argument('--chunks', help = 'number of chunks to split it into.',required = True, dest = 'chunks')
    a = parser.parse_args()

    # basic error checking to kick out, if your IO is busted.
    if(not os.path.isfile(a.input)):
        sys.exit("Can't find input file: {}".format(a.input))
    if(not os.path.isdir(a.outdir)):
        sys.exit("Can't find output dir: {}".format(a.outdir))
    
    ncIn = Dataset(a.input,'r')

    SN = len(ncIn.dimensions['south_north'])
    SN_stag = len(ncIn.dimensions['south_north_stag'])
    WE = len(ncIn.dimensions['west_east'])
    chnk = int(SN)//int(a.chunks)
    chnk_stag = int(SN_stag)//int(a.chunks)
    remainder = int(SN)%int(a.chunks)
    remainder_stag = int(SN_stag)%int(a.chunks) 
    i = 0
    ii = 0 
    for c in range(int(a.chunks)):
        #zz.append(z[i:i+chnk])
    
        outName = os.path.join(a.outdir, os.path.basename(a.input)+'_section_{:02d}'.format(c+1))
        ncOut = Dataset(outName, "w", format="NETCDF4")
        # copy all attributes
        ncOut.setncatts(ncIn.__dict__)
        # copy all dimensions if time, make it just 1 'cause that's what we're doing.
        for name, dimension in ncIn.dimensions.items():
            if name == 'south_north':
                ncOut.createDimension(
                    name, (chnk))
            elif name == 'south_north_stag':
                ncOut.createDimension(
                    name, (chnk_stag))
 
            else:
                ncOut.createDimension(
                    name, (len(dimension) if not dimension.isunlimited() else 1))
            
        for name, variable in ncIn.variables.items():
            x = ncOut.createVariable(name, variable.datatype, variable.dimensions)
            if('south_north' in variable.dimensions):
                sn_pos = np.where(np.asarray(variable.dimensions)=='south_north')[0][0]
                print(name, variable.dimensions, np.where(np.asarray(variable.dimensions)=='south_north')[0][0] )
                if(sn_pos == 1):
                    ncOut[name][:,:,:] = ncIn[name][:,i:i+chnk,:]
                elif(sn_pos == 2):
                    ncOut[name][:,:,:,:] = ncIn[name][:,:,i:i+chnk,:]
            elif ('south_north_stag' in variable.dimensions):
                sn_pos = np.where(np.asarray(variable.dimensions)=='south_north_stag')[0][0]
                print(name, variable.dimensions, np.where(np.asarray(variable.dimensions)=='south_north_stag')[0][0] )
                if(sn_pos == 1):
                    ncOut[name][:,:,:] = ncIn[name][:,ii:ii+chnk_stag,:]
                elif(sn_pos == 2):
                    ncOut[name][:,:,:,:] = ncIn[name][:,:,ii:ii+chnk_stag,:]
            else:
                ncOut[name][:] = ncIn[name][:] 
            # copy variable attributes all at once via dictionary
            ncOut[name].setncatts(ncIn[name].__dict__)
        ncOut.close()
        i = i+chnk
        ii = ii+chnk_stag
        if c == int(a.chunks) - 2:
            chnk = chnk + remainder
            chnk_stag = chnk_stag + remainder_stag
    
    ncIn.close() 
