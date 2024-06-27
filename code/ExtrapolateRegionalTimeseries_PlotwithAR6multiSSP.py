#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 16:12:52 2024

@author: ashleybellas
"""

import numpy as np
import matplotlib.pyplot as plt

from ConstructRateAccelEnsemble import ConstructRateAccelEnsemble
from ComputeExtrapolations_Preferred_and_ConfidenceLimits import ComputeExtrapolations_Preferred_and_ConfidenceLimits
from IsolateRegionalTimeseriesAR6_5PointPolygon import IsolateRegionalTimeseriesAR6_5PointPolygon
from SetExtrapolationReferenceYear import SetExtrapolationReferenceYear


#%% region and corrections applied to the data

regions = 'boxsetup1'
data    = 'case1'

#%% compute errors 

firstrun = 1      # set to 0 to load rather than recompute   
N_dist   = 10000  # no. of ensemble members

# compute an ensemble of rate and accel based on errors from GIA, measurement, and serial correlation
rdist, adist, lbls = ConstructRateAccelEnsemble(data,regions,N_dist,firstrun)


#%% compute preferred and 90% confidence extrpolations based on rate & accel ensemble


file            = np.load('data/ssh_grids_v2205_timeseries_'+data+'_monthly.npz')   # orig. gridded data set from which rate&accel were computed
t               = file.f.t
t0              = np.mean(t)                                                        # the time at which the rate was solved

confidencelimit = 0.9                                                               # desired confidence limit in decimal units

sla_pref, r_uconflim, r_lconflim, t = ComputeExtrapolations_Preferred_and_ConfidenceLimits(data,regions,N_dist,rdist,adist,confidencelimit,t0)


#%% load AR6 data 

firstrun = 1

# compute regional average timeseries and save to files
if firstrun==1:
    IsolateRegionalTimeseriesAR6_5PointPolygon('ssp126',regions)
    IsolateRegionalTimeseriesAR6_5PointPolygon('ssp245',regions)
    IsolateRegionalTimeseriesAR6_5PointPolygon('ssp370',regions)
    IsolateRegionalTimeseriesAR6_5PointPolygon('ssp585',regions)

# load the regional AR6 timeseries    
data4     = np.load('data/ts_'+regions+'_AR6ssp126.npz')
t_ssp     = data4['t']
s_ssp126  = data4['sbar']
e05ssp126 = data4['ebar_05']
# e95ssp126 = data4['ebar_95']

data5     = np.load('data/ts_'+regions+'_AR6ssp245.npz')
s_ssp245  = data5['sbar']
# e05ssp245 = data5['ebar_05']
# e95ssp245 = data5['ebar_95']

data6     = np.load('data/ts_'+regions+'_AR6ssp370.npz')
s_ssp370  = data6['sbar']
# e05ssp370 = data6['ebar_05']
# e95ssp370 = data6['ebar_95']

data7     = np.load('data/ts_'+regions+'_AR6ssp585.npz')
s_ssp585  = data7['sbar']
# e05ssp585 = data7['ebar_05']
e95ssp585 = data7['ebar_95']

#%% translate extrapolations to be relative to a desired year

YEAR0 = 2020

sla_prefYEAR0   = SetExtrapolationReferenceYear(t,sla_pref,YEAR0)

r_uconflimYEAR0 = SetExtrapolationReferenceYear(t,r_uconflim,YEAR0)
r_lconflimYEAR0 = SetExtrapolationReferenceYear(t,r_lconflim,YEAR0)

s_ssp126YEAR0   = SetExtrapolationReferenceYear(t_ssp,s_ssp126,YEAR0)
s_ssp245YEAR0   = SetExtrapolationReferenceYear(t_ssp,s_ssp245,YEAR0)
s_ssp370YEAR0   = SetExtrapolationReferenceYear(t_ssp,s_ssp370,YEAR0)
s_ssp585YEAR0   = SetExtrapolationReferenceYear(t_ssp,s_ssp585,YEAR0)

e05ssp126YEAR0  = SetExtrapolationReferenceYear(t_ssp,e05ssp126,YEAR0)
e95ssp585YEAR0  = SetExtrapolationReferenceYear(t_ssp,e95ssp585,YEAR0)


#%% Plot relative to 2020 in 2 columns 

plt.rcParams['figure.dpi'] = 600
plt.rcParams['savefig.dpi'] = 600
plt.figure(figsize=(12,14))
plt.rcParams.update({'font.size': 18})

Nbox = sla_prefYEAR0.shape[0]
for n in range(Nbox):
    
    ax = plt.subplot(4,2,n+1)
    
    # altimeter extrapolation
    plt.plot(t, sla_prefYEAR0[n,:], color='black', label='extrapolated satellite altimetry',linewidth=4.0,zorder=4)
    ax.fill_between(t, (sla_prefYEAR0[n,:]+r_uconflimYEAR0[n,:]), (sla_prefYEAR0[n,:]-r_lconflimYEAR0[n,:]), color='black', alpha=0.25)
    
    # AR6 projections 
    plt.plot(t_ssp, s_ssp126YEAR0[n,:]-e05ssp126YEAR0[n,:], color='blue',linewidth=1.0,linestyle='--',label='min(AR6(SSP-1,2.6),90%)')
    
    plt.plot(t_ssp, s_ssp126YEAR0[n,:], color='blue',  linewidth=2.0,linestyle='-',label='AR6(SSP-1,2.6)')
    plt.plot(t_ssp, s_ssp245YEAR0[n,:], color='lime',  linewidth=2.0,linestyle='-', label='AR6(SSP-2,4.5)')
    plt.plot(t_ssp, s_ssp370YEAR0[n,:], color='orange',linewidth=2.0,linestyle='-', label='AR6(SSP-3,7.0)')
    plt.plot(t_ssp, s_ssp585YEAR0[n,:], color='red',   linewidth=2.0,linestyle='-', label='AR6(SSP-5,8.5)')
    
    plt.plot(t_ssp, s_ssp585YEAR0[n,:]+e95ssp585YEAR0[n,:], color='red',linewidth=1.0,linestyle='--',label='max(AR6(SSP-5,8.5),90%)')
    
    # annot.
    plt.title(lbls[n], color='black')
    plt.grid('on', alpha=0.5)
    plt.xlim([2020,2050])
    ax.set_xticks([2020,2025,2030,2035,2040,2045,2050],['2020','','2030','','2040','','2050'])
    plt.ylim([-10, 300])
    ax.set_yticks([0,50,100,150,200,250,300])
    #
    # if n==6: plt.legend(bbox_to_anchor=(1.55, 0.6))
    if n==6: plt.legend(bbox_to_anchor=(1.15, 0.85), loc=2, borderaxespad=0.)
    if n==0 or n==2 or n==4 or n==6:plt.ylabel('mm')
    else:ax.set_yticklabels([])
    if n<5:ax.set_xticklabels([])
    else:plt.xlabel('Year')
    
plt.savefig('plots/extrpltSLA_MEaSUREs_'+data+'_multiSSP_'+regions+'_t0'+str(YEAR0)+'.jpeg', bbox_inches='tight', format='jpeg')



