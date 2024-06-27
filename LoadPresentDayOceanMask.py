#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 13 12:19:22 2024

@author: ashleybellas
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def LoadPresentDayOceanMask(grid_lat, grid_lon):
   
    if np.min(grid_lon) >= 0:
        file0 = np.load('data/PresentDayOceanMask_0360.npz')
   
    if np.min(grid_lon) < 0:
        file0 = np.load('data/PresentDayOceanMask_180180.npz')
   
    olat  = file0['lat']
    olon  = file0['lon']
    omask = file0['mask']
   
    # interpolate onto the desired grid
    interp_omask = RegularGridInterpolator((olat, olon), omask)
    omaski       = interp_omask((grid_lat, grid_lon))
    omaski[omaski> 0.5] = 1
    omaski[omaski<=0.5] = 0

    return omaski
