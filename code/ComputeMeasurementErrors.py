#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  2 10:43:50 2024

@author: ashleybellas

The following functions were copied from Thmoas Frederikse's github assocaited with Nerem et al 2022:
    https://github.com/thomasfrederikse/2021_GMSL_acc/tree/main

*I have modified the code to produce only an ensemble in the rate & accel associated with measurement errors. 
*This will be combined with the GIA & serially correlated formal errors in another script.

"""

import numpy as np
import os


def ComputeMeasurementErrors(inputfile,intputdir,outputfile,N):
    settings = def_settings(inputfile,intputdir,N)
    gsl = read_gsl(settings)

    trend_ptb_ensemble_mat = np.zeros([settings['num_ens'],3])
    
    ptb_gsl_alt = compute_gsl_ensemble(gsl,settings)

    trend_ptb_ensemble_mat = compute_trend_acc_ptb_ensemble(gsl, ptb_gsl_alt, settings)
    
    trend_stats = compute_statistics(trend_ptb_ensemble_mat, settings)
    
    # save result to file 
    np.savez('data/'+outputfile, trend_ptb_ensemble_mat=trend_ptb_ensemble_mat)

    return trend_ptb_ensemble_mat, trend_stats


def def_settings(inputfile,inputdir,N):
    settings = {}
    # Directories
    settings["dir_project"] = inputdir
    settings["dir_data"] = settings["dir_project"] + 'data/'

    # Input files
    settings['fn_gmsl'] = settings["dir_data"] + inputfile

    # Define time steps
    data = np.load(settings["fn_gmsl"])
    settings['time_sat'] = data['t']

    # Number of ensembles for Monte Carlo simulation
    settings['num_ens'] = N
    return settings


def read_gsl(settings):
    
    data = np.load(settings["fn_gmsl"])
    gsl  = data['sla_']  
    
    return gsl


def compute_trend_acc_ptb_ensemble(gsl,ptb_gsl_alt,settings):
    # Compute the trend and acceleration and the resulting extrapolation

    # Create an ensemble of GMSL estimates: GMSL = GSL + measurement errror perturbations
    obs_ensemble = gsl[np.newaxis,:] + ptb_gsl_alt

    # Remove 2000-2020 mean
    baseline_idx = (settings['time_sat']>2000)&(settings['time_sat']<2020)
    obs_ensemble -= obs_ensemble[:,baseline_idx].mean(axis=1)[:,np.newaxis]

    # Estimate the trend and acceration and
    trend_ensemble     = np.zeros([settings['num_ens'],3])
    trend_ptb_ensemble = np.zeros([settings['num_ens'],3])
    
    amat = np.ones([len(settings['time_sat']),3])
    amat[:,1] = settings['time_sat'] - np.mean(settings['time_sat'])
    amat[:,2] = 0.5*(settings['time_sat'] - np.mean(settings['time_sat']))**2

    amat_T  = amat.T
    amat_sq = np.linalg.inv(np.dot(amat_T, amat))
    
    trend_gsl = np.dot(amat_sq, np.dot(amat_T, gsl))    # 0: bias, 1: trend; 2: accel; - Ashley
    
    for i in range(settings['num_ens']):
        trend_ensemble[i,:] = np.dot(amat_sq, np.dot(amat_T, obs_ensemble[i,:]))    # 0: bias, 1: trend; 2: accel; - Ashley
        trend_ptb_ensemble[i,:] = trend_ensemble[i,:] - trend_gsl 

    return(trend_ptb_ensemble)



def compute_gsl_ensemble(gsl,settings):
    # Generate an ensemble of global-mean geocentric sea level curve
    # uncertainties around the mean, corrected for ENSO/Pinatubo

    # Term 1: Ablain et al: measurement uncertainties
    t_rand_highf  = gen_autocov(1.7,1.5, 1.2, 2/12, settings)
    t_rand_medf   = gen_autocov(1.3,1.2, 1, 1, settings)
    t_rand_wettr  = gen_autocov(1.1,1.1, 1.1, 5, settings)
    t_rand_lorbit = gen_autocov(1.12,0.5, 0.5, 10, settings)
    t_rand_intmis = gen_randjump(settings)
    t_rand_dorbit = gen_drift(0.1,settings)
    t_rand_tpx    = gen_tpxdrift(settings)
    ptb_gsl_alt   = (t_rand_highf + t_rand_medf + t_rand_wettr + t_rand_lorbit + t_rand_intmis + t_rand_dorbit+t_rand_tpx).astype(np.float32)
    
    return ptb_gsl_alt



def compute_statistics(trend_ensemble, settings):
    # Compute statistics

    trend_stats = {}
    trend_stats['percentiles'] = [5, 17, 50, 83, 95]
    trend_stats["trend"] = np.percentile(trend_ensemble[:, 1], trend_stats['percentiles'])
    trend_stats["accel"] = np.percentile(trend_ensemble[:, 2], trend_stats['percentiles'])

    return(trend_stats)

# -------------------------------------------------------------------
# Ablain et al. (2018) functions
# These functions are used to generate an ensemble of GMSL
# estimates with the uncertainty structure from Ablain et al. (2018)
# -------------------------------------------------------------------
def gen_randjump(settings):
    tjump = [(1999 + 1.5 / 12),(2002 + 4.5 / 12),(2008 + 8.5 / 12),(2016 + 8.5 / 12)]
    hjump = [2,0.5,0.5,0.5]
    t_rand = np.zeros([settings['num_ens'],len(settings['time_sat'])])
    for jump in range(len(tjump)):
        rndjump = np.random.randn(settings['num_ens'])*hjump[jump]
        t_acc= settings['time_sat']>tjump[jump]
        t_rand[:,t_acc] += rndjump[:,np.newaxis]
    return(t_rand)

def gen_drift(drift,settings):
    t_rand = drift * np.random.randn(settings['num_ens'])[:,np.newaxis] * (settings['time_sat']-settings['time_sat'].mean())[np.newaxis, :]
    return(t_rand)

def gen_tpxdrift(settings):
    t_rand = np.zeros([settings['num_ens'], len(settings['time_sat'])])
    t_txa = settings['time_sat'] < (1999 + 1.5 / 12)
    t_txb = (settings['time_sat'] > (1999 + 1.5 / 12)) & (settings['time_sat']<(2002 + 4.5 / 12))
    t_rand[:, t_txa] = (np.random.randn(settings['num_ens']) * 0.7/12)[:,np.newaxis]
    t_rand[:, t_txb] = (np.random.randn(settings['num_ens']) * 0.1/12)[:,np.newaxis]
    t_rand = np.cumsum(t_rand, axis=1)
    return t_rand


def gen_autocov(sig_TPX, sig_J1, sig_J23, l_factor, settings):
    jump_0 = np.argmin(np.abs(settings['time_sat']-2002+4.5/12))
    jump_1 = np.argmin(np.abs(settings['time_sat']-2008+8.5/12))
    sig_array = np.zeros(len(settings['time_sat']))
    sig_array[:jump_0] = sig_TPX
    sig_array[jump_0:jump_1] = sig_J1
    sig_array[jump_1:] = sig_J23
    t_distance = np.abs(settings['time_sat'][:,np.newaxis] - settings['time_sat'][np.newaxis,:])
    covmat = np.exp(-0.5*(t_distance/l_factor)**2) + np.eye(len(settings['time_sat']))*1.0e-10
    covc = np.linalg.cholesky(covmat)

    t_rand = np.zeros([settings['num_ens'], len(settings['time_sat'])])
    for n in range(settings['num_ens']):
        t_rand[n, :] = sig_array*np.matmul(covc, np.random.randn(len(settings['time_sat'])))
    return t_rand


