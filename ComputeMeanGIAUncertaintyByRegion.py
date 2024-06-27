#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:19:22 2024

@author: ashleybellas
"""


import numpy as np
from   scipy.interpolate import griddata
from LoadPresentDayOceanMask import LoadPresentDayOceanMask

def ComputeMeanGIAUncertaintyByRegion(regions):
    
    
    #%% load data 

    file0 = str(("data/GIA_maps_Caron_Ivins_2019.txt"))
    f0    = open(file0,"r").readlines()
    n     = int(sum(1 for line in f0))
    hdr   = 6
    sig_dg = np.zeros((n-hdr))
    colat  = np.zeros((n-hdr))
    lon    = np.zeros((n-hdr))
    pnts   = np.zeros((n-hdr,2))
    for k in range(n-hdr):
        colat[k]     = f0[k+hdr].split('\t')[0]
        lon[k]       = f0[k+hdr].split('\t')[1]
        sig_dg[k]    = f0[k+hdr].split('\t')[5]
        pnts[k,0]    = colat[k]
        pnts[k,1]    = lon[k]
    # put points onto grid 
    grid_colat, grid_lon = np.mgrid[0.: 179.:180j, 0.: 359.:360j]
    sig_dg_grd = griddata(pnts, np.array(sig_dg), (grid_colat, grid_lon), method='linear')
    # set any values on land = NaN
    grid_lat = 90.-grid_colat
    lat = grid_lat[:,0]
    lon = grid_lon[0,:]
    mask = LoadPresentDayOceanMask(grid_lat, grid_lon)   # 1 in ocean, 0 on land 
    mask[mask==0] = np.nan
    sig_dg_grd *= mask
    
    #%% load regions 
    
    box  = np.load('data/'+regions+'.npz')
    lat0 = box['lat0']
    lat1 = box['lat1']
    lon0 = box['lon0']
    lon1 = box['lon1']
    lbls = box['labels']
    Nbox = lat0.size
    # check for 2nd order box (e.g., Kiribati)
    if len(box.files)>=7:
        lat2 = box['lat2']
        lon2 = box['lon2']
        # check for 5-point polygon (e.g., boxetup7)
    if len(box.files)==11:
        lat3 = box['lat3']
        lat4 = box['lat4']
        lon3 = box['lon3']
        lon4 = box['lon4']
    
    ### reshape if lon bounds of regions exceed 360
    if np.max(lon1>360):
        sig_dg_grd2 = np.copy(sig_dg_grd)
        sig_dg_grd2[grid_lon<360-(np.max(lon1)-360)] = sig_dg_grd[grid_lon>=np.max(lon1)-360]
        sig_dg_grd2[grid_lon>=360-(np.max(lon1)-360)]  = sig_dg_grd[grid_lon<np.max(lon1)-360]
        sig_dg_grd = np.copy(sig_dg_grd2)
        grid_colat, grid_lon = np.mgrid[0.: 179.:180j, np.max(lon1)-360: 360+(np.max(lon1)-360):360j]
        lon = grid_lon[0,:]
    
    
    #%% compute mean GIA in each region 
    
    dlat = 1.*np.deg2rad(1)
    dlon = 1.*np.deg2rad(1)
    sig_dg_rgnl = np.zeros((Nbox))
    
    # 1st or 2nd order box ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    if len(box.files)<=7:  
        for n in range(Nbox):
            # 1st order box
            if len(box.files)==5:
                ind     = (grid_lat>=lat0[n])*(grid_lat<=lat1[n]) * (grid_lon>=lon0[n])*(grid_lon<=lon1[n])
            # 2nd order box 
            elif len(box.files)==7:
                ind     = ((grid_lat>=lat0[n])*(grid_lat<=lat1[n]) * (grid_lon>=lon0[n])*(grid_lon<=lon1[n])) * ~((grid_lat>=lat2[n])*(grid_lat<=lat1[n]) * (grid_lon>=lon2[n])*(grid_lon<=lon1[n]))
            # computen mean GIA uncertatinty    
            dA         = 1.*dlat * 1.*np.cos(grid_lat[ind]*np.deg2rad(1))*dlon  * (1.-np.isnan(sig_dg_grd[ind]))
            A          = np.sum(dA)
            sig_dg_rgnl[n]  = np.nansum(dA*sig_dg_grd[ind])/A

    # 5-point polygon ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    elif len(box.files)==11:  # 5-point polygon 
        ind = np.zeros((Nbox,lat.size,lon.size),dtype=int)
        for n in range(Nbox):
            ind_lat = np.where((lat>=lat0[n]) * (lat<=lat3[n]))[0]        # true for ptype = 1 2 or 3
            # query the polygon type and define left & right boundaries 
            if lat2[n]==lat3[n]:
                if lat4[n]==lat0[n] and lon4[n]==lon0[n]:
                    # print('ptype=1')
                    lb = (lon3[n]-lon4[n])/(lat3[n]-lat4[n])*(lat[ind_lat]-lat4[n]) + lon4[n]
                    rb = (lon2[n]-lon1[n])/(lat2[n]-lat1[n])*(lat[ind_lat]-lat1[n]) + lon1[n]
                else:
                    # print('ptype=3')
                    latlb = lat[ind_lat]
                    latlba = latlb[latlb<lat4[n]]
                    latlbb = latlb[latlb>=lat4[n]]
                    lba = (lon4[n]-lon0[n])/(lat4[n]-lat0[n])*(latlba-lat0[n]) + lon0[n]
                    lbb = (lon3[n]-lon4[n])/(lat3[n]-lat4[n])*(latlbb-lat4[n]) + lon4[n]
                    lb = np.hstack((lba,lbb))
                    rb = (lon2[n]-lon1[n])/(lat2[n]-lat1[n])*(lat[ind_lat]-lat1[n]) + lon1[n]
            elif lat3[n]==lat4[n]:
                # print('ptype=2')
                lb = (lon4[n]-lon0[n])/(lat4[n]-lat0[n])*(lat[ind_lat]-lat0[n]) + lon0[n]
                latrb = lat[ind_lat]
                latrba = latrb[latrb<lat2[n]]
                latrbb = latrb[latrb>=lat2[n]]
                rba = (lon2[n]-lon1[n])/(lat2[n]-lat1[n])*(latrba-lat1[n]) + lon1[n] 
                rbb = (lon3[n]-lon2[n])/(lat3[n]-lat2[n])*(latrbb-lat2[n]) + lon2[n] 
                rb = np.hstack((rba,rbb))
            # construct an n x m index matrix for points that are within the boundaries 
            k = 0
            for i in ind_lat:
                ind[n,i,:] = (lon>lb[k])*(lon<rb[k])
                k+=1
            # compute mean GIA uncertainty in each region
            dA         = 1.*dlat * 1.*np.cos(grid_lat*np.deg2rad(1))*dlon  * (1.-np.isnan(sig_dg_grd)) * ind[n,:,:]
            A          = np.sum(dA)
            # print(A)
            sig_dg_rgnl[n]  = np.nansum(dA*sig_dg_grd)/A

    #%% write to .txt
    
    file1   = str('data/mean_GIA_std_'+regions+'.txt')
    f1=open(file1,"w")
    f1.writelines(['The average std of the GIA correction, in terms of mm/year EWH :\n' ])
    for i in range(Nbox):
        f1.writelines(['%s\t%.3E\n' % (lbls[i], sig_dg_rgnl[i]) ])
    f1.close()
    
    
