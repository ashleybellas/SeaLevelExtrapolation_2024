#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 11:45:23 2024

@author: ashleybellas
"""

####################################################################################
# Computes the extrapolated sea level curves associated with the preferred trajectory,
# and the X% confidence limits above and below this preferred trajectory based
# on the ensemble of rate and accel produced by considering errors from various sources
# (separate). 
# 
#   INPUTS
#   confidencelimit: confidence limit in decimal units (i.e., for 90% confience limits, input 0.9)
####################################################################################
import numpy as np

def ComputeExtrapolations_Preferred_and_ConfidenceLimits(data,regions,N_dist,rdist,adist,confidencelimit,t0):

    #%% extrapolate the N_dist rate & accel pairs 
    
    data0    = np.load('data/ssh_grids_v2205_'+data+'_rate_accel_'+regions+'.npz')
    rate     = data0['rate']      # [mm]
    accel    = data0['accel']
    Nbox     = rate.size
    
    t        = np.arange(t0,2101)
    sla      = np.zeros((Nbox, N_dist, t.size))
    sla_pref = np.zeros((Nbox, t.size))
    for i in range(Nbox):
        # compute the preferred extrapolation         
        sla_pref[i] = rate[i]*(t-np.min(t)) + 1/2*accel[i]*(t-np.min(t))**2
        for j in range(N_dist):
            # compute the ensemble extrapolation         
            sla[i,j,:] = rdist[i,j]*(t-np.min(t)) + 1/2*adist[i,j]*(t-np.min(t))**2
    
    #%% compute the confidence limits as f(t), 
    ### based on the percentage of cases that fall between the mean and an upper (r_uconflim) or lower (r_lconflim) limit
    
    N_conflim = confidencelimit*N_dist/2
    #
    r_uconflim = 0.01*np.ones((Nbox,t.size))      # init. 
    r_lconflim = 0.01*np.ones((Nbox,t.size))      # init. 
    for it in range(1,t.size):
        # 
        n_conflimu=np.zeros((Nbox))
        n_lconflim=np.zeros((Nbox))
        for i in range(Nbox):
            # upper confidence limit bound
            while n_conflimu[i]<N_conflim:
                n_conflimu[i] = np.sum( (sla[i,:,it]>sla_pref[i,it]) * (sla[i,:,it]<=sla_pref[i,it]+r_uconflim[i,it]) )    
                r_uconflim[i,it]*=1.01
            # lower confidence limit bound
            while n_lconflim[i]<N_conflim:
                n_lconflim[i] = np.sum( (sla[i,:,it]<sla_pref[i,it]) * (sla[i,:,it]>=sla_pref[i,it]-r_lconflim[i,it]) )   
                r_lconflim[i,it]*=1.01


    return sla_pref, r_uconflim, r_lconflim, t