#!/usr/bin/env python3
import matplotlib
import argparse, configparser, os, sys, h5py, glob
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from netCDF4 import Dataset
# my libraries
from lib.graphics.maps import plotMap
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'Run a wrf simulation through CRTM.')
    parser.add_argument('--in', help = 'path to ncdiag input', required = True, dest = 'input')
    parser.add_argument('--out', help = 'path to ncdiag input', required = True, dest = 'output')
    a = parser.parse_args()
    inDir = a.input
    outDir = a.output
    if(not os.path.isdir(inDir)):
        sys.exit('Input Directory not present')
    if(not os.path.isdir(outDir)):
        sys.exit('Output Directory not present')
    filez = glob.glob( os.path.join(inDir,'wrf*crtm_brightness_temperature*.nc4') )
    filez.sort()
    for i,f in enumerate(filez):
        print(f)
        ncIn = Dataset(f,'r')
        Tb = ncIn['Tb']
        wavenumbers = ncIn['wavenumbers']
        lat = np.asarray(ncIn['XLAT']).flatten()
        lon = np.asarray(ncIn['XLONG']).flatten()
        for idx in range(len(wavenumbers)):
            plotMap( lat, lon, np.asarray(Tb[0,:,:,idx]).flatten(),\
                    'BT map {} cm'.format(wavenumbers[idx]), os.path.join(outDir,'step{:04d}_Tb{:04d}'.format(i,idx)),\
                    units='K', colorScheme = 'viridis',
                    extent = [ min(lon), max(lon), min(lat), max(lat) ])
