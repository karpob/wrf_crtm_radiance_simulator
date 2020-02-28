#!/usr/bin/env python3
import argparse, configparser, os, sys
from netCDF4 import Dataset
import numpy as np



def main(input, outputDir, nchunks):
    sectionPosition = input.find('section_01')
    if(sectionPosition == -1):
        print("This wasn't a run based off my splitter!")
        sys.exit()
    filesToProcess = []
    filesToProcess.append(input)
    for i in range(1,nchunks):
        filesToProcess.append( input.replace('section_01','section_{:02d}'.format(i+1)) )
    print("combining these files:")
    southNorthSizes = [] 
    for f in filesToProcess:
        v = Dataset(f,'r')
        print(f)
        southNorthSizes.append( v.getncattr('SOUTH-NORTH_PATCH_END_UNSTAG') )
    if(any(southNorthSizes != southNorthSizes[0])):
        print('wrong north south dimension!')
        sys.exit()
    # create a netcdf file removing "section_NN"
    ncOut = Dataset(input.replace('_section_01',''), 'w')
    # copy all attributes
    ncIn = Dataset(a.input,'r')
    ncOut.setncatts(ncIn.__dict__)
    #carry forward all dimensions, except south_north, revert that back to original wrf file.
    for name, dimension in ncIn.dimensions.items():
        if name == 'south_north':
            #put back original size
            ncOut.createDimension(
                name, southNorthSizes[0])
        else:
            ncOut.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else 1))
    i = 0
    ncIn.close()
    for f in filesToProcess:
        ncIn = Dataset(f,'r')
        for name, variable in ncIn.variables.items():
            if ( i == 0 ): x = ncOut.createVariable(name, variable.datatype, variable.dimensions)
            chnk = len(ncIn.dimensions['south_north'])
            if('south_north' in variable.dimensions):
                sn_pos = np.where(np.asarray(variable.dimensions)=='south_north')[0][0]
                if('nchannels' in variable.dimensions):
                    ncOut[name][:,i:i+chnk,:,:] = ncIn[name][:,:,:,:]
                else:
                    ncOut[name][:,i:i+chnk,:] = ncIn[name][:,:,:]
            else:
                ncOut[name][:] = ncIn[name][:] 
        i += chnk
        ncIn.close()
    ncOut.close()
    print("Done!")
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'Combine Brightness temperature files into one contigous file aligned with grid.')
    parser.add_argument('--in', help = 'path to first input file', required = True, dest = 'input')
    parser.add_argument('--output-dir', help = 'path to ncdiag output', required = True, dest = 'outdir')
    parser.add_argument('--chunks', help = 'number of chunks to split it into.',required = True, dest = 'chunks')
    a = parser.parse_args()
    main(a.input, a.outdir, int(a.chunks))


