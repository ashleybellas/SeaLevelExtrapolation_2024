#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:02:48 2023

@author: ashleybellas
"""

########################################################################################
# Isolates MEaSUREs data within previously defined regions
# Computes the mean timeseries of sea level in each region
# Saves timeseries to file 


# INPUTS

    # data:     tells what corrections have been applied to the MEaSUREs data 
    # regions:  predefined regions of interest
    # plot:     whether to plot the isolated region outlines     
########################################################################################

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from   mpl_toolkits.basemap import Basemap

def IsolateTimeseries_5PointPolygonRegions(data, regions, plot):
    #%% load data 
    
    # sea level anomaly data
    file = np.load('data/ssh_grids_v2205_timeseries_'+data+'_monthly.npz')
    lat  = file.f.lat
    lon  = file.f.lon
    t    = file.f.t
    sla  = np.float64(file.f.sla)
    print('N time: '+str(np.shape(t)))
    print('N lat: '+str(np.shape(lat)))
    print('N lon: '+str(np.shape(lon)))
    print('N SL: '+str(np.shape(sla)))
    #
    grid_lat, grid_lon = np.mgrid[np.min(lat):np.max(lat):480j,np.min(lon): np.max(lon):1080j]
    
    #%% load regions
    
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
    lbls = box['labels']
    Nbox = lat0.size
    
    #%% reshape if necessary 
    
    if np.max(lon1>360):
        sla2 = np.copy(sla)
        sla2[:, grid_lon< 360-(np.max(lon1)-360)] = sla[:, grid_lon>=np.max(lon1)-360]
        sla2[:, grid_lon>=360-(np.max(lon1)-360)] = sla[:, grid_lon< np.max(lon1)-360]
        sla  = np.copy(sla2)

        grid_lat, grid_lon = np.mgrid[np.min(lat):np.max(lat):480j,np.max(lon1)-360: 360+(np.max(lon1)-360):1080j]
        lon = grid_lon[0,:]

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
    sla_ = np.zeros((Nbox, t.size))
    
    for n in range(Nbox):
        for i in range(t.size):
            dA         = 1.*dlat * 1.*np.cos(grid_lat*np.deg2rad(1))*dlon  * (1.-np.isnan(sla[i,:,:])) * ind[n,:,:]
            A          = np.sum(dA)
            sla_[n,i]  = np.nansum(dA*sla[i,:,:])/A
        print('done with box %i' % n)

    #%% plot to check regions are extracted correctly 

    if plot==1:
        for n in range(Nbox):
            
            plt.rcParams['figure.dpi'] = 600
            plt.rcParams['savefig.dpi'] = 600
            plt.figure(figsize=(3,4))   # initialize figure
            plt.rcParams.update({'font.size': 10})
    
            cmap = plt.cm.seismic
            cmap.set_bad(color='grey')
            levels = np.arange(-5, 5.1, 0.5)
            norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=False)
            
            map = Basemap(projection='merc',llcrnrlat=-70,urcrnrlat=70, llcrnrlon=20,urcrnrlon=380,lat_ts=20,resolution='c')
            map.drawcoastlines(linewidth=0.1,color='gray')
            map.drawmeridians(np.arange(90,380,90),linewidth=0.3,color='white',labels=[True,True,False,True],fontsize=8.)
            map.drawparallels(np.arange(-60,61,30),linewidth=0.3,color='white',labels=[True],fontsize=8.)
    
            x,y = map(grid_lon,grid_lat)
            map.pcolormesh(x,y,sla[n,:,:]*ind[n,:,:],cmap=cmap,norm=norm,shading='auto')
    
            x,y = map(lon0[n],lat0[n])
            map.scatter(x,y,5.,zorder=4)
            x,y = map(lon1[n],lat1[n])
            map.scatter(x,y,5.,zorder=4)
            x,y = map(lon2[n],lat2[n])
            map.scatter(x,y,5.,zorder=4)
            x,y = map(lon3[n],lat3[n])
            map.scatter(x,y,5.,zorder=4)
            x,y = map(lon4[n],lat4[n])
            map.scatter(x,y,5.,zorder=4)
            
    #%% save extracted timeseries to file 
    
    np.savez('data/ssh_grids_v2205_'+data+'_timeseries_'+regions+'.npz', t=t, sla_=sla_, labels=lbls)
    
    return t, sla_, lbls
    
