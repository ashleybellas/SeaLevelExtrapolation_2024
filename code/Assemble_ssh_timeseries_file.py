#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 11:50:35 2024

@author: ashleybellas
"""

import numpy as np
import xarray as xr
import os
import datetime


# define how many lat/lat cells to skip when concatenating data from the native 1/4 deg grid
dx      = 2
# download grid at first time step
fpath   = '/Users/ashleybellas/Documents/Education & Research/CUB Postdoc Semester 1&2/SLA_gridded/jpl_gridded_ssha_2205/'
string1 = 'ssh_grids_v2205_'
ds      = xr.open_dataset(fpath+string1+'1992101012.nc')
lat     = ds.Latitude.data[::dx]
lon     = ds.Longitude.data[::dx]
xr.Dataset.close(ds)

# loop through data to concatenate grids, avg monthly
years  = np.arange(1993,2023,1)
months = np.arange(1,13,1)
days   = np.arange(1,32,1)
ii=0
jj=0
for year in years:
    print(np.round(100 * ii / len(years),1), '% complete')
    for month in months:
        kk=0.
        sea_level_monthly = 0
        for day in days:
            file_mn_yr = string1 + str(int(year)) + '{:02d}'.format(month) + '{:02d}'.format(day) +'12.nc'
            FN = os.listdir(fpath)
            if file_mn_yr in FN:
                ds   = xr.open_dataset(fpath+file_mn_yr)
                dtii = datetime.datetime.strptime(str(int(year))+'-'+'{:02d}'.format(int(month))+'-15', '%Y-%m-%d').timetuple().tm_yday
                sea_level_monthly += np.expand_dims(ds.SLA.data[0,::dx,::dx],axis= 0)
                time_monthly = year+(dtii/365.25)
                kk+=1.
                xr.Dataset.close(ds)
        sea_level_monthly/=kk
        if jj==0:
            sea_level = sea_level_monthly
            time = time_monthly
        else:
            sea_level = np.concatenate((sea_level,sea_level_monthly), axis =0 )
            time = np.hstack((time,time_monthly))
        jj+=1
    ii+=1

# convert units to [mm]
sea_level*=1.e3
print('N time: '+str(np.shape(time)))
print('N lat: '+str(np.shape(lat)))
print('N lon: '+str(np.shape(lon)))
print('N SL: '+str(np.shape(sea_level)))

np.savez('data/ssh_grids_v2205_timeseries_case1_monthly.npz', lat=lat, lon=lon, t=time, sla=sea_level)

