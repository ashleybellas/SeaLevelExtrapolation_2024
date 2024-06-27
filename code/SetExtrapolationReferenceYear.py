#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 25 12:37:42 2024

@author: ashleybellas
"""

#############################################################################
# Subtracts the value of the extrapolation at YEAR0 from all entries of the timeseries
# 
#   INPUTS
#   t:      time stamps for the timeseries 
#   ts:     the timeseries with dimensions [region, time]
#   YEAR0:  the desired reference year for the extrapolation 
#############################################################################


import numpy as np

def SetExtrapolationReferenceYear(t,ts,YEAR0):
    
    ind_t     = np.where(np.abs(t-YEAR0) == np.min(np.abs(t-YEAR0)))[0][0]
    ts2D      = np.zeros((ts.shape[0],1))
    ts2D[:,0] = ts[:,ind_t]    
    tsYEAR0   = ts - ts2D
    
    return tsYEAR0
    
    
    