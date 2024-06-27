#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 15:30:19 2024

@author: ashleybellas
"""

########################################################################################

# First coordinate (lon0,lat0) must be the lower left corner 
# List coordinates counter-clockwise

########################################################################################

import numpy as np 

#                IO    NP    TP'    SP     NA'    SA      AC
lon0 = np.array([20,   94,  125,    150,   290,  290,    20])
lat0 = np.array([-45,  20,  -20,    -45,   5,    -45,   -70])

lon1 = np.array([150,  255,  290,   290,   350,  380,   380])
lat1 = np.array([-45,  20,   -20,   -45,   5,    -45,   -70])

lon2 = np.array([85,  280,  290,    290,   380,  380,   380])
lat2 = np.array([30,   70,    5,    -20,   70,     5,   -45])

lon3 = np.array([20,   120,  255,   125,   280,  290,    20])
lat3 = np.array([30,   70,    20,   -20,   70,     5,   -45])

lon4 = np.array([20,   94,   94,   150,    255,  290,    20])
lat4 = np.array([-45,  20,   20,   -45,    20,   -45,   -70])

labels = ['Indian Ocean', 'North Pacific', 'Tropical Pacific', 'South Pacific', 'North Atlantic', 'South Atlantic', 'Antarctic Circumpolar']

# save to file 
np.savez('data/boxsetup1.npz', lat0=lat0,lon0=lon0, lat1=lat1,lon1=lon1,lat2=lat2,lon2=lon2, lat3=lat3,lon3=lon3, lat4=lat4,lon4=lon4, labels=labels)
