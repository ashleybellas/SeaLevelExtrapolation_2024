#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:10:53 2024

@author: ashleybellas

"""

#####################################################################################
# Computes the serial correlation inflation factor follow the methods of Maul & Martin 1993 

# INPUTS
    # data:     tells what corrections have been applied to the MEaSUREs data 
    # regions:  predefined regions of interest
#####################################################################################    

import numpy as np 
import statsmodels.api as sm

def ComputeSerialCorrelationScaleFactor(data,regions):


    # load the residuals (data-fit)
    file   = np.load('data/ssh_grids_v2205_'+data+'_rate_accel_'+regions+'.npz')
    sla_ds = file['ds'] 
    Nbox   = sla_ds.shape[0] 
    
    # compute the autocorrelation sequence with lag=  0,1 and take the result for lag=1
    acs  = np.zeros((Nbox))
    for n in range(Nbox):
        acs[n] = sm.tsa.acf(sla_ds[n,:], nlags=1)[1]     
             
    # compute the scale factor following Maul & Martin 1993
    Fi = 1/np.sqrt((1-acs)/(1+acs))
    
    # save results to file 
    np.savez('data/ScaleFMaul_MEaSUREs_'+data+'_'+regions+'.npz', Fi=Fi)
    
    return Fi