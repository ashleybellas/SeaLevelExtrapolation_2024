#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 14:55:18 2024

@author: ashleybellas
"""

########################################################################################
# Computes the rate and acceleration based on MEaSUREs data isolated within predefined regions 

# INPUTS

    # data:     tells what corrections have been applied to the MEaSUREs data 
    # regions:  predefined regions of interest
    # t0:       time at which rate is solved
########################################################################################

import numpy as np 

def ComputeRateAccel_RegionalTimeseries(data,regions,t0):
    
    file = np.load('data/ssh_grids_v2205_'+data+'_timeseries_'+regions+'.npz')
    t    = file['t']
    sla_ = file['sla_']
    lbls = file['labels']
    Nbox = lbls.size

    #%% use least squares to fit rate&accel, formal errors
    
    nt = t.size
    pi = np.pi
    
    t -= t0
    t1 = np.zeros((t.size,1))
    t1[:,0] = t
    
    H = np.hstack([np.ones((nt,1)), t1, 0.5*t1**2, np.sin(2*pi*t1), np.cos(2*pi*t1), np.sin(4*pi*t1), np.cos(4*pi*t1)])
    
    # fit parameters
    rate     = np.zeros((Nbox))
    accel    = np.zeros((Nbox))
    errR     = np.zeros((Nbox))
    errA     = np.zeros((Nbox))
    sla_dszn = np.zeros((Nbox, t.size))
    dsla_    = np.zeros((Nbox, t.size))
    fit_quad = np.zeros((Nbox, t.size))
    for n in range(Nbox):
        # compute rate & accel
        c        = np.matmul(np.linalg.pinv(H),sla_[n,:])
        rate[n]  = c[1]
        accel[n] = c[2]
        
        # compute formal error
        fit_full      = c[0] + (c[1]*t) + (0.5*c[2]*t**2) + (c[3]*np.sin(2*t*np.pi)) + (c[4]*np.cos(2*t*np.pi)) + (c[5]*np.sin(4*t*np.pi)) + (c[6]*np.cos(4*t*np.pi))
        eps     = np.nanstd(sla_[n,:]-fit_full)
        Hi      = np.linalg.inv(H.T.dot(H))
        err     = np.sqrt((eps**2)*np.diag(Hi))
        errR[n],errA[n]    = err[1],err[2]
        
        # save the residuals so you can scale formal error to total error (separate)
        dsla_[n,:]    = sla_[n,:]-fit_full
        fit_quad[n,:] = c[0] + (c[1]*t) + (0.5*c[2]*t**2)
        
        # deseason
        sznl = np.matmul(H[:,3:7],c[3:7])
        sla_dszn[n,:] = sla_[n,:] - sznl 

    #%% save extracted rate & accel
    
    np.savez('data/ssh_grids_v2205_'+data+'_rate_accel_'+regions+'.npz', rate=rate, accel=accel, errR=errR, errA=errA, ds=dsla_,labels=lbls)
    
    return sla_dszn, fit_quad, rate, accel 