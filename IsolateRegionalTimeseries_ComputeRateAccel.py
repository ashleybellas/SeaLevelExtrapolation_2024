#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:31:06 2024

@author: ashleybellas
"""


import numpy as np 
import matplotlib.pyplot as plt
from IsolateTimeseries_5PointPolygonRegions import IsolateTimeseries_5PointPolygonRegions
from ComputeRateAccel_RegionalTimeseries import ComputeRateAccel_RegionalTimeseries

#%% Isolate regional timeseries

firstrun = 1                    # change this to 0 if you already ran this segment
data     = 'CSEOFtest3mode123'
regions  ='boxsetup1'
plotcheckregions = 1            # change this to 0 if you don't want to plot

if firstrun == 1:
    t, sla_, lbls = IsolateTimeseries_5PointPolygonRegions(data, regions, plotcheckregions)

else:
    file = np.load('data/ssh_grids_v2205_'+data+'_timeseries_'+regions+'.npz')
    t    = file['t']
    sla_ = file['sla_']
    lbls = file['labels']

Nbox = lbls.size

#%% compute rate & accel, error, etc. 

t0=2007.5     # solve rate at the this time
sla_dszn, fit_quad, rate, accel = ComputeRateAccel_RegionalTimeseries(data,regions,t0)


#%% plot

plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
plt.figure(figsize=(15,5))
plt.rcParams.update({'font.size': 14})
    
for n in range(Nbox):
     
    ax=plt.subplot(2,4,n+1)  
    plt.scatter(t, sla_dszn[n,:], 3.5, color='red', marker='o', label='data')
    # plt.scatter(t, sla_[n,:], 3.5, color='red', marker='o', label='data')
    plt.plot(t, fit_quad[n,:], color='black',label='a+bt+$\frac{1}{2}$c$t^2$')
    plt.grid('on')
    
    # annot.
    plt.title(str(lbls[n]))
    plt.annotate(str('rate = %.1f mm/yr ' % (rate[n])), color='black', xy=(0.025, .85), xycoords='axes fraction', horizontalalignment='left')
    plt.annotate(str('accel = %.3f mm/yr$^2$' % (accel[n])), color='black', xy=(0.025, .70), xycoords='axes fraction', horizontalalignment='left')
    #axes
    ax.set_xticks([2000,2010,2020])
    if n>2:
        plt.xlabel('Year')
        ax.set_xticklabels(['2000','2010','2020'])
    else:ax.set_xticklabels(['','',''])
    plt.ylim([-50,100])
    if n==0 or n==4: plt.ylabel('mm')
    else: ax.set_yticklabels([])
    
# save file
plt.savefig('plots/MEaSUREs_'+data+'_ts_'+regions+'.jpeg', bbox_inches='tight', format='jpeg')