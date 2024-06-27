#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 16:20:50 2024

@author: ashleybellas
"""

import numpy as np
import xarray as xr
from scipy.interpolate import griddata

def IsolateRegionalTimeseriesAR6_5PointPolygon(AR6case,regions):
    #%% load data 
    print('Begin regional avg of %s' % AR6case)
    
    # AR6 projections
    data0 = xr.open_dataset('data/AR6/total_'+AR6case+'_medium_confidence_values.nc')
    sla   = data0['sea_level_change'].values  # [mm]
    lat   = data0['lat'].values
    lon   = data0['lon'].values
    qnt   = data0['quantiles'].values
    t     = data0['years'].values
    # translate from -180:180 to 0:360
    lon[lon<0] +=360
    
    # regions
    box  = np.load('data/'+regions+'.npz')
    lat0 = box['lat0']
    lat1 = box['lat1']
    lat2 = box['lat2']
    lat3 = box['lat3']
    lat4 = box['lat4']
    lon0 = box['lon0']
    lon1 = box['lon1']
    lon2 = box['lon2']
    lon3 = box['lon3']
    lon4 = box['lon4']
    Nbox = lat0.size
    
    
    #%% take mean across quantiles for preferred value (quantiles have equal probability)
    sla1  = np.mean(sla,axis=0)
    
    #%% compute 90% confidence bounds as distance between median and 0.95 or .05 quantiles 
    
    sig_95 = np.zeros((t.size, lat.size))
    sig_05 = np.zeros((t.size, lat.size))
    sig_83 = np.zeros((t.size, lat.size))
    sig_17 = np.zeros((t.size, lat.size))
    for i in range(t.size):
        for j in range(lat.size):
            sig_05[i,j] = sla[qnt==0.5,i,j] - sla[qnt==0.05,i,j]
            sig_95[i,j] = sla[qnt==0.95,i,j] - sla[qnt==0.5,i,j]
            sig_17[i,j] = sla[qnt==0.5,i,j] - sla[qnt==0.17,i,j]
            sig_83[i,j] = sla[qnt==0.83,i,j] - sla[qnt==0.5,i,j]
     
    #%% interpolate onto a 2d grid
    
    grid_lat, grid_lon = np.mgrid[-89.5:89.5:180j, 0.5:359.5:360j]
    
    pnts      = np.zeros((lat.size,2))
    pnts[:,0] = lat
    pnts[:,1] = lon
    
    sla_grd    = np.zeros((t.size,grid_lat.shape[0],grid_lat.shape[1]))
    sig_05_grd = np.zeros((t.size,grid_lat.shape[0],grid_lat.shape[1]))
    sig_95_grd = np.zeros((t.size,grid_lat.shape[0],grid_lat.shape[1]))
    sig_17_grd = np.zeros((t.size,grid_lat.shape[0],grid_lat.shape[1]))
    sig_83_grd = np.zeros((t.size,grid_lat.shape[0],grid_lat.shape[1]))
    for i in range(t.size):
        sla_grd[i, :,:]    = griddata(pnts, np.array(sla1[i,:]),   (grid_lat, grid_lon), method='linear')
        sig_05_grd[i, :,:] = griddata(pnts, np.array(sig_05[i,:]), (grid_lat, grid_lon), method='linear')
        sig_95_grd[i, :,:] = griddata(pnts, np.array(sig_95[i,:]), (grid_lat, grid_lon), method='linear')
        sig_17_grd[i, :,:] = griddata(pnts, np.array(sig_17[i,:]), (grid_lat, grid_lon), method='linear')
        sig_83_grd[i, :,:] = griddata(pnts, np.array(sig_83[i,:]), (grid_lat, grid_lon), method='linear')
    
    lat = grid_lat[:,0]
    lon = grid_lon[0,:]
    
    #%% reshape if necessary (function Reshape will decide)
 
    temp, temp, sla_grd            = Reshape(grid_lat, grid_lon, sla_grd, lon1)
    temp, temp, sig_05_grd         = Reshape(grid_lat, grid_lon, sig_05_grd, lon1)
    temp, temp, sig_95_grd         = Reshape(grid_lat, grid_lon, sig_95_grd, lon1)
    temp, temp, sig_17_grd         = Reshape(grid_lat, grid_lon, sig_17_grd, lon1)
    grid_lat, grid_lon, sig_83_grd = Reshape(grid_lat, grid_lon, sig_83_grd, lon1)
    
    #%% isolate indices within a 5-point polygon
    
    ind = np.zeros((Nbox,lat.size,lon.size),dtype=int)
    
    for n in range(Nbox):
    # for n in range(4,5):
        ind_lat = np.where((lat>=lat0[n]) * (lat<=lat3[n]))[0]        # true for ptype = 1 2 or 3
        # query the polygon type and define left & right boundaries ~~~~~~~~~~~~~~#
        if lat2[n]==lat3[n]:
            if lat4[n]==lat0[n] and lon4[n]==lon0[n]:
                print('ptype=1')
                lb = (lon3[n]-lon4[n])/(lat3[n]-lat4[n])*(lat[ind_lat]-lat4[n]) + lon4[n]
                rb = (lon2[n]-lon1[n])/(lat2[n]-lat1[n])*(lat[ind_lat]-lat1[n]) + lon1[n]
            else:
                print('ptype=3')
                latlb = lat[ind_lat]
                latlba = latlb[latlb<lat4[n]]
                latlbb = latlb[latlb>lat4[n]]
                lba = (lon4[n]-lon0[n])/(lat4[n]-lat0[n])*(latlba-lat0[n]) + lon0[n]
                lbb = (lon3[n]-lon4[n])/(lat3[n]-lat4[n])*(latlbb-lat4[n]) + lon4[n]
                lb = np.hstack((lba,lbb))
                rb = (lon2[n]-lon1[n])/(lat2[n]-lat1[n])*(lat[ind_lat]-lat1[n]) + lon1[n]
        elif lat3[n]==lat4[n]:
            print('ptype=2')
            lb = (lon4[n]-lon0[n])/(lat4[n]-lat0[n])*(lat[ind_lat]-lat0[n]) + lon0[n]
            latrb = lat[ind_lat]
            latrba = latrb[latrb<lat2[n]]
            latrbb = latrb[latrb>lat2[n]]
            rba = (lon2[n]-lon1[n])/(lat2[n]-lat1[n])*(latrba-lat1[n]) + lon1[n] 
            rbb = (lon3[n]-lon2[n])/(lat3[n]-lat2[n])*(latrbb-lat2[n]) + lon2[n] 
            rb = np.hstack((rba,rbb))
        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
        # construct an n x m index matrix for points that are within the boundaries 
        k = 0
        for i in ind_lat:
            ind[n,i,:] = (lon>lb[k])*(lon<rb[k])
            k+=1
            
    #%% compute mean times series in each region
    
    dlat = (lat[1] - lat[0])*np.deg2rad(1)
    dlon = (lon[1] - lon[0])*np.deg2rad(1)
    sbar = np.zeros((Nbox, t.size))
    ebar_05 = np.zeros((Nbox, t.size))
    ebar_95 = np.zeros((Nbox, t.size))
    ebar_17 = np.zeros((Nbox, t.size))
    ebar_83 = np.zeros((Nbox, t.size))
    #
    for n in range(Nbox):
        #
        for i in range(t.size):
            dA        = 1.*dlat * 1.*np.cos(grid_lat*np.deg2rad(1))*dlon  * (1.-np.isnan(sla_grd[i,:,:])) * ind[n,:,:]
            A         = np.sum(dA)  
            sbar[n,i] = np.nansum(dA*sla_grd[i,:,:])/A
            ebar_05[n,i] = np.nansum(dA*sig_05_grd[i,:,:])/A
            ebar_95[n,i] = np.nansum(dA*sig_95_grd[i,:,:])/A
            ebar_17[n,i] = np.nansum(dA*sig_17_grd[i,:,:])/A
            ebar_83[n,i] = np.nansum(dA*sig_83_grd[i,:,:])/A
    
    
    #%% save extracted timeseries
    
    np.savez('data/ts_'+regions+'_AR6'+AR6case+'.npz', t=t, sbar=sbar, ebar_05=ebar_05, ebar_95=ebar_95, ebar_17=ebar_17, ebar_83=ebar_83, lat=lat, lon=lon)
    
    print('Done with %s' % AR6case)
    
    
    
def Reshape(grid_lat, grid_lon, data, lon1):
    
    ### reshape if lon bounds of regions exceed 360
    if np.max(lon1>360):
        data2 = np.copy(data)
        data2[:,grid_lon<=360-(np.max(lon1)-360)] = data[:,grid_lon>=np.max(lon1)-360]
        data2[:,grid_lon>360-(np.max(lon1)-360)]  = data[:,grid_lon<np.max(lon1)-360]
        data = np.copy(data2)
        grid_lat, grid_lon = np.mgrid[np.min(grid_lat):np.max(grid_lat):180j,  np.max(lon1)-360: 360+(np.max(lon1)-360):360j]
        
    return grid_lat, grid_lon, data