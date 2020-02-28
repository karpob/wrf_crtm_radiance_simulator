#!/usr/bin/env python3
import IPython
import argparse, configparser, os, sys, h5py
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
import numpy as np
from netCDF4 import Dataset
from wrf import getvar
import gc
# my libraries
from lib.graphics.maps import plotMap
from lib.graphics.profile import plotContour, plotLines
from lib.graphics.linePlots import basicLine
from lib.pycrtm.pyCRTM import pyCRTM, profilesCreate
from lib.pycrtm.crtm_io import readTauCoeffODPS 
from lib.pycrtm.units import  waterPpmvDry2GmKgDry, gasPpmvMoistToDry, gasKgKgMoistToDry
from lib.pycrtm.interpolation import profileInterpolate 

def logLinear(x, xo, yo):
    """
    Do a log-linear interpolation.
    """
    logX =  np.log(x)
    logXo = np.log(xo)
    logYo = np.log(yo)
    return np.exp(np.interp(logX, logXo, logYo))

def pressureFromEta(pTop,pSurf,eta):
    """
    Get Level/Layer pressure from model top pressure, surface pressure, eta value
    pTop -- top/lowest pressure in model (interface top)-  units- whatever you like
    pSurf -- surface pressure (bottom) -  units - whatever you like
    eta -- eta values from model - units - none
    returns Pressure associated with Eta Level/Layer
    """

    Pout = pTop + eta.T*(pSurf - pTop)
    Pout = Pout.T
    return Pout

def readProfileItemsH5( filename, additionalItems = []):
    """
    Read an RTTOV-style atmosphere profile.
    In: filename to hdf5
    Out: Pressure, Temperature, CO2, O3 [nprofiles,nlevels]
    Out: Gas_Units (mass, ppmv dry, ppmv moist)
    """
    items = ['T','Q','O3']
    if(len(additionalItems)>0):
        for i in additionalItems:
            items.append(i)
    h5 = h5py.File( filename )
    groups = list(h5['PROFILES'].keys())
    nprofiles = len(groups)
    nlevs, = np.asarray( h5['PROFILES'][groups[0]]['P'] ).shape 
    P = np.zeros([nprofiles,nlevs])
    itemsOut = {}
    for i in items: itemsOut[i] = np.zeros([nprofiles,nlevs])
    for i,g in enumerate(groups):
        P[i,:] = np.asarray(h5['PROFILES'][g]['P'])
        for ii in items:
            itemsOut[ii][i,:] = np.asarray(h5['PROFILES'][g][ii])
        GasUnits = int(np.asarray(h5['PROFILES'][g]['GAS_UNITS']))
    
    return P, itemsOut, GasUnits 
def getProfilesWrfStep(nc, step):
    p = {}
    Play = np.asarray(getvar(nc, 'pressure', timeidx=step))
    Play = np.flip( Play.reshape(Play.shape[0],Play.shape[1]*Play.shape[2]), axis=0)
    p['P'] = Play.T # need to transpose to make nprofiles, nlay

    # Psfc,Ptop outside wrf interface units convert from Pa to mbar
    Psfc = np.asarray(getvar(nc,'PSFC', timeidx=step))/100.0
    Ptop = np.asarray(nc.variables['P_TOP'][step])/100.0 #P_TOP made wrf interface angry...just access it through regular NC interface
    # eta on interfaces.
    ZNW = np.asarray(getvar(nc, 'ZNW', timeidx = step))
    Psfc = Psfc.reshape([Psfc.shape[0]*Psfc.shape[1],1])
    Plev = pressureFromEta(Ptop,Psfc,ZNW)
    Plev = np.flip(Plev, axis=0)
    p['Pi'] = Plev.T 
    print('fixing Q')
    Q = np.asarray(getvar(nc, 'QVAPOR', timeidx=step)) #kg/kg moist air  
    Q[np.where(Q<0.0)] = 0.0

    Q = 1000.0*gasKgKgMoistToDry(Q,Q) # to get 1000 g/kg which CRTM loves, probably going all the way back to lblrtm / Bill Gallery's lblatm trainwreck.
    Q[np.where(Q<0.0)] = 0.0
    print('Done fixing Q!')
    Q = np.flip(Q.reshape(Q.shape[0],Q.shape[1]*Q.shape[2]), axis=0)
    p['Q'] = Q.T
   
    T = np.asarray(getvar(nc, 'tk', timeidx = step))
    T = np.flip(T.reshape(T.shape[0],T.shape[1]*T.shape[2]), axis=0)
    p['T'] = T.T

    s = {}

    windSpeed10m, windDirection10m =  getvar(nc, 'wspd_wdir10', timeidx=step) 
    windSpeed10m = np.asarray(windSpeed10m)
    s['windSpeed10m'] = windSpeed10m.reshape(windSpeed10m.shape[0]*windSpeed10m.shape[1])

    U10 = np.asarray( getvar(nc,'U10') ) 
    V10 = np.asarray( getvar(nc,'V10') ) 
    spd = np.sqrt(U10**2+V10**2)
    windDirection10m = np.degrees(np.arctan2(U10/spd,V10/spd)) % 360
    s['windDirection10m'] = windDirection10m.reshape(windDirection10m.shape[0]*windDirection10m.shape[1])

    sst = np.asarray( getvar(nc,'SST', timeidx=step) )
    sst = sst.reshape(sst.shape[0]*sst.shape[1])
    # surface temps go land, water, snow, ice
    s['surfaceTemperatures'] = np.zeros([sst.shape[0],4])
    s['surfaceTemperatures'][:,0] = sst 
    s['surfaceTemperatures'][:,1] = sst 
    s['surfaceTemperatures'][:,2] = sst 
    s['surfaceTemperatures'][:,3] = sst 
    
    a = {}
    lon = np.asarray( getvar(nc,'lon',timeidx=step) )
    a['lon'] = lon.reshape(lon.shape[0]*lon.shape[1])
     
    lat = np.asarray( getvar(nc,'lat',timeidx=step) )
    a['lat'] = lat.reshape(lat.shape[0]*lat.shape[1])

    cloudFraction = np.asarray( getvar(nc,'cloudfrac', timeidx=step) ) 
    cloudFraction = cloudFraction.reshape(cloudFraction.shape[0],cloudFraction.shape[1]*cloudFraction.shape[2])
    a['cloudFraction'] =  cloudFraction # documentation says it's low, mid, high clouds (whatever that means)
    return p, s, a

def setProfilesCRTM(h5_ppmv, nc, step):
    # First read in Standard profiles to get climatological value for Ozone (not in WRF Files). 
    # Can we run CRTM without specifying ozone? I vaguely remember reading somewhere you have to always
    # specify ozone so we will here. 
    Pi, profileItems, gas_units = readProfileItemsH5(h5_ppmv)

    for i in list(profileItems.keys()):
        if ( i != 'T' ):
            profileItems[i] = gasPpmvMoistToDry(profileItems[i], profileItems['Q']) 
    climatologyO3 = profileItems['O3'][0,:]
    climatologyP = Pi[0,:]
    Pi = []
    profileItems = []
    #####
    # End Ozone climatology read in. Will need to interpolate later on here.
    #####
    prof, surf, anc = getProfilesWrfStep(nc, step)
    nprof = prof['P'].shape[0]
    nlay  = prof['P'].shape[1]
    profilesCRTM = profilesCreate( nprof, nlay, nAerosols=0, nClouds=0 )
    o3Out = np.zeros([nprof,nlay])
    for i in range(nprof):
        o3Out[i,:] = logLinear(prof['P'][i,:], climatologyP, climatologyO3)

    profilesCRTM.O3[:,:] = o3Out
    for i in list(prof.keys()):
        exec( "profilesCRTM.{}[:,:] = prof['{}']".format(i,i) )
    for i in list(surf.keys()):
        if surf[i].ndim>1: 
            exec( "profilesCRTM.{}[:,:] = surf['{}']".format(i,i) )
        else:   
            exec( "profilesCRTM.{}[:] = surf['{}']".format(i,i) )
    profilesCRTM.Angles[:,:] = 0.0
    profilesCRTM.Angles[:,2] = 100.0  # Solar Zenith Angle 100 degrees zenith below horizon.

    profilesCRTM.DateTimes[:,0] = 2017
    profilesCRTM.DateTimes[:,1] = 9
    profilesCRTM.DateTimes[:,2] = 1

    # Turn off Aerosols and Clouds
    #profilesCRTM.aerosolType[:] = -1
    #profilesCRTM.cloudType[:] = -1

    profilesCRTM.surfaceFractions[:,:] = 0.0
    profilesCRTM.surfaceFractions[:,1] = 1.0 # all water!
    profilesCRTM.S2m[:,1] = 35.0 # just use salinity out of S2m for the moment.

    # land, soil, veg, water, snow, ice
    profilesCRTM.surfaceTypes[:,3] = 1
    #print('min O3', 'min H2O') 
    #print(profilesCRTM.O3.min().min(), profilesCRTM.Q.min().min()) 
    return profilesCRTM, anc
if __name__ == "__main__":
    parser = argparse.ArgumentParser( description = 'Run a wrf simulation through CRTM.')
    parser.add_argument('--in', help = 'path to ncdiag input', required = True, dest = 'input')
    parser.add_argument('--output-dir', help = 'path to ncdiag output', required = True, dest = 'outdir')
    parser.add_argument('--instrument',help = 'instrument name to process', required = False, dest='instrument', default='iasi_metop-a')
    parser.add_argument('--nthreads', help="number of threads", dest='nthreads', type=int, required=False, default=16 )
    parser.add_argument('--idx', help="time idx", dest='idx', type=int, required=True )
    a = parser.parse_args()

    # basic error checking to kick out, if your IO is busted.
    if(not os.path.isfile(a.input)):
        sys.exit("Can't find input file: {}".format(a.input))
    if(not os.path.isdir(a.outdir)):
        sys.exit("Can't find output dir: {}".format(a.outdir))
    

    #################################################################################################
    # Get installed path to coefficients from pycrtm submodule install (crtm.cfg in pycrtm directory)
    # load stuff we need for CRTM coefficients
    #################################################################################################
    pathToThisScript = os.path.dirname(os.path.abspath(__file__))
    pathInfo = configparser.ConfigParser()
    pathInfo.read( os.path.join(pathToThisScript,'lib','pycrtm','crtm.cfg') )
    coefficientPathCrtm = pathInfo['CRTM']['coeffs_dir']
    #################################################################################################
    # Get CRTM coefficient interface levels, and pressure layers 
    # Pull pressure levels out of coefficient
    # get pressures used for profile training in CRTM. 
    #################################################################################################
    tauCoeff = os.path.join(coefficientPathCrtm,'iasi_metop-a.TauCoeff.bin'.format(a.instrument))
    if ( not os.path.isfile(tauCoeff)):
        sys.exit( "Can't find sensor TauCoeff file:".format(tauCoeff) )
    crtmTauCoef, _ = readTauCoeffODPS( tauCoeff )
    coefLevCrtm = np.asarray(crtmTauCoef['level_pressure'])
    layerPressuresCrtm = np.asarray(crtmTauCoef['layer_pressure'])

    ##########################
    # Set Profiles
    ##########################
    rttovPath = 'rttov/rttov/' 
    h5_ppmv =  os.path.join(rttovPath,'rttov_test','profile-datasets-hdf','standard101lev_allgas.H5')
    print('input file {}'.format(a.input)) 
    
    ncIn = Dataset(a.input,'r')
    sarr = []
    for b in ncIn.variables['Times'][a.idx]:
        sarr.append(b.decode("utf-8"))
    timestamp = ''.join(sarr)
    #outName = os.path.join(a.outdir, os.path.basename(a.input)+'_'+timestamp+"_crtm_brightness_temperature_stiched")
    outName = os.path.join(a.outdir, os.path.basename(a.input)+'_crtm_brightness_temperature')
    ncOut = Dataset(outName, "w", format="NETCDF4")
    # copy all attributes
    ncOut.setncatts(ncIn.__dict__)
    # copy all dimensions if time, make it just 1 'cause that's what we're doing.
    for name, dimension in ncIn.dimensions.items():
        ncOut.createDimension(
            name, (len(dimension) if not dimension.isunlimited() else 1))
    for name, variable in ncIn.variables.items():
        if name in ['Times','XLAT','XLONG']:
            x = ncOut.createVariable(name, variable.datatype, variable.dimensions)
            print(name, ncIn[name].shape,ncOut[name].shape)
            if(ncOut[name].ndim == 2):
                ncOut[name][:,:] = ncIn[name][a.idx,:]
            else:ncOut[name][:,:,:] = ncIn[name][a.idx,:,:]
            # copy variable attributes all at once via dictionary
            ncOut[name].setncatts(ncIn[name].__dict__)
    SN = len(ncIn.dimensions['south_north'])
    WE = len(ncIn.dimensions['west_east'])



    crtmOb = pyCRTM()
    crtmOb.StoreTrans = False
    crtmOb.coefficientPath = coefficientPathCrtm 
    crtmOb.sensor_id = 'iasiB2_metop-a'
    crtmOb.loadInst()
    crtmOb.computeTauLevels = False
    crtmOb.nThreads = int(a.nthreads)
    profilesCRTM, anc = setProfilesCRTM(h5_ppmv, ncIn, a.idx)
    ncIn.close()
    ncIn = []
    gc.collect()
    crtmOb.profiles = profilesCRTM
    crtmOb.runDirect()
    Tb1 =  np.asarray(crtmOb.Bt)
    Wn1 = np.asarray(crtmOb.Wavenumbers)
    crtmOb = []
    gc.collect()
    idx, = np.where(Wn1 >= 1700)
    Tb1 = Tb1[:,idx]
    Wn1 = Wn1[idx]
    crtmOb2 = pyCRTM()
    crtmOb2.StoreTrans = False
    crtmOb2.coefficientPath = coefficientPathCrtm 
    crtmOb2.sensor_id = 'iasiB3_metop-a' 
    crtmOb2.loadInst()
    crtmOb2.computeTauLevels = False
    crtmOb2.nThreads = int(a.nthreads)
    crtmOb2.profiles = profilesCRTM
    crtmOb2.runDirect()
    profilesCRTM = []
    gc.collect()
    Tb2 =  np.asarray(crtmOb2.Bt)
    Wn2 = np.asarray(crtmOb2.Wavenumbers)
    idx, = np.where(Wn2 <= 2500)
    Tb2 = Tb2[:,idx]
    Wn2 = Wn2[idx]
 
    crtmOb2 = []
    Wn = np.concatenate([Wn1, Wn2])

    print('done direct')
    print('Writing Brightness Temps.')

 
    ncOut.createDimension('nchannels', len(Wn) )
    outTb = ncOut.createVariable('Tb',np.float32,('Time','south_north', 'west_east','nchannels') ) 
    #lowFrac = ncOut.createVariable('cloudFracLow',np.float32,('Time','south_north','west_east') )
    #midFrac = ncOut.createVariable('cloudFracMid',np.float32,('Time','south_north','west_east') )
    #highFrac = ncOut.createVariable('cloudFracHigh',np.float32,('Time','south_north','west_east') )
    wavenumbersOut = ncOut.createVariable('wavenumbers',np.float32,('nchannels'))
    Tb1 = Tb1[:,:].reshape([SN,WE,len(Wn1)])
    Tb2 = Tb2[:,:].reshape([SN,WE,len(Wn2)])
    outTb[0,:,:,0:len(Wn1)] = Tb1[:,:,:]
    outTb[0,:,:,len(Wn1):len(Wn1)+len(Wn2)] = Tb2[:,:,:]
    """
    for i in range(len(Wn1)):
        outTb[:,:,:,i] = Tb1[:,i].reshape([SN,WE])
        ncOut.sync()
    offset = len(Wn1)
    Tb1 = []
    gc.collect()
    for i in range(len(Wn2)):
        outTb[:,:,:,offset+i] = Tb2[:,i].reshape([SN,WE])
        ncOut.sync()
    Tb2 = []
    """
    gc.collect()

    print("Writing Ancillary Data.")
    wavenumbersOut[:] = Wn
    #lowFrac[:,:,:], midFrac[:,:,:], highFrac[:,:,:] = getvar(ncIn,'cloudfrac',timeidx=a.idx)
    print("Done!")
