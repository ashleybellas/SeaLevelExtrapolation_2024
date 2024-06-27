#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:30:24 2024

@author: ashleybellas
"""

#########################################################################################
# Computes an ensemble of rates and accelerations based on errors from
#   1. GIA, 
#   2. serially correlated formal errors, and 
#   3. measurement errors
# 
# INPUTS 
    # data:     tells what corrections have been applied to the MEaSUREs data 
    # regions:  predefined regions of interest
    #  N_dist:  no. of members in the ensemble 
    # firstrun: tells whether first run through or not, option to load versus compute 
#########################################################################################

import numpy as np

from ComputeSerialCorrelationScaleFactor import ComputeSerialCorrelationScaleFactor
from ComputeMeanGIAUncertaintyByRegion import ComputeMeanGIAUncertaintyByRegion
from ComputeMeasurementErrors import ComputeMeasurementErrors


def ConstructRateAccelEnsemble(data,regions,N_dist,firstrun):
    
    #%% serially correlated formal errors
    
    print('Begin serially correlated formal errors.')
    if firstrun==1:# compute
        Fi   = ComputeSerialCorrelationScaleFactor(data,regions)
    else:# load from saved file
        data1 = np.load('data/ScaleFMaul_MEaSUREs_'+data+'_'+regions+'.npz')
        Fi   = data1['Fi']
    
    # load and scale the formal errors for lag-1 serial correlation 
    data0    = np.load('data/ssh_grids_v2205_'+data+'_rate_accel_'+regions+'.npz')
    rate     = data0['rate']      # [mm]
    accel    = data0['accel']
    errR_srl = data0['errR']*Fi   # 1 sigma
    errA_srl = data0['errA']*Fi
    lbls     = data0['labels']
    Nbox     = rate.size
    
    # construct an ensemble of errors from serially correlated formal errors
    rdist_srl = np.zeros((Nbox, N_dist))
    adist_srl = np.zeros((Nbox, N_dist))
    for i in range(Nbox):
        rdist_srl[i,:] = np.random.normal(0.,errR_srl[i],N_dist)
        adist_srl[i,:] = np.random.normal(0.,errA_srl[i],N_dist)
    
    print('End serially correlated formal errors.')
    
    #%% measurement errors 
    
    print('Begin measurement errors.')
    if firstrun==1:         
        inputfile  = 'ssh_grids_v2205_GMSLtimeseries_monthly_'+data+'.npz'
        inputdir   = ''
        outputfile = 'global_uniform_measurement_error_perturbations.npz'
        ComputeGMSLTimeseries(data)
        trend_ptb_ensemble_mat, trend_stats = ComputeMeasurementErrors(inputfile,inputdir,outputfile,N_dist)
    else:
        data2 = np.load('data/global_uniform_measurement_error_perturbations.npz')
        trend_ptb_ensemble_mat = data2['trend_ptb_ensemble_mat']
    
    rdist_msr = trend_ptb_ensemble_mat[:,1]
    adist_msr = trend_ptb_ensemble_mat[:,2]
    
    print('End measurement errors.')
    
    #%% GIA
    
    print('Begin GIA errors.')
    # compute and/or load the regionally averaged GIA errors
    if firstrun==1:
        ComputeMeanGIAUncertaintyByRegion(regions)
    data3 = 'data/mean_GIA_std_'+regions+'.txt'
    f0   = open(data3,"r").readlines()
    N    = int(sum(1 for line in f0))
    errR_gia  = np.zeros((Nbox))
    for n in range(1,N):
        errR_gia[n-1] = f0[n].split('\t')[1]
    
    # construct an ensemble of errors from GIA 
    rdist_gia = np.zeros((Nbox, N_dist))
    for i in range(Nbox):    
        rdist_gia[i,:] = np.random.normal(0.,errR_gia[i],N_dist)
    
    print('End GIA errors.')    
    
    #%% construct 10k pairs of rate & accel by adding perturbations to the mean values 
    
    print('Begin construct full error ensemble.')
    rdist = np.zeros([Nbox, N_dist])
    adist = np.zeros([Nbox, N_dist])
    for i in range(Nbox):
        rdist[i,:] = rate[i] + rdist_srl[i,:] + rdist_gia[i,:] + rdist_msr
        adist[i,:] = accel[i] + adist_srl[i,:] + adist_msr
    
    print('End construct full error ensemble.')
    
    #%% outputs 
    
    return rdist, adist, lbls


def ComputeGMSLTimeseries(data):
    
    data1 = np.load('data/ssh_grids_v2205_timeseries_'+data+'_monthly.npz')
    lon   = data1['lon']*np.deg2rad(1)
    lat   = data1['lat']*np.deg2rad(1)
    time  = data1['t']
    sla   = data1['sla']
    
    # average
    sla_ = np.zeros(time.size)
    dlat = (lat[1] - lat[0])
    dlon = (lon[1] - lon[0])
    sla_ = np.zeros((time.size))
    A    = 0
    #
    for i in range(time.size):
        for j in range(lat.size):
            dA       = 1.*dlat * 1.*np.cos(lat[j])*dlon  * (1.-np.isnan(sla[i,j,:]))
            sla_[i] += np.nansum(dA*sla[i,j,:])
            if i==0:
                A += np.sum(dA)
    sla_/=A
    np.savez('data/ssh_grids_v2205_GMSLtimeseries_monthly_'+data+'.npz', sla_=sla_, t=time, lat=lat, lon=lon)