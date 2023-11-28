#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 27 19:10:14 2023

@author: u1318104
"""

import numpy as np
import os

#%%

Pers=[24,30,40,50,60,70,80,90,100]

LonMin=-113
LonMax=-107

LatMin=40
LatMax=45

fdir=os.path.dirname(__file__)

dir_in=fdir+'Data_1psi'
dir_out=fdir+'Cut_1psi'

per=60

Lons,Lats,Vphs,A1,Psi1=np.loadtxt(dir_in+'USArr_%ds_1psi_0.txt'%per,unpack=True,usecols=(0,1,2,3,4))

fp=[]
for i in range(len(Lons)):
    if Lons[i]<LonMin or Lons[i]>LonMax or Lats[i]<LatMin or Lats[i]>LatMax:
        continue
    # fp.append('%f %f %f %f %f %f %f'%(Lon[i],Lat[i],CCi[i],AA1[i],Psi1[i],AA2[i],Psi2[i]))
    fp.append('%f %f %f '%(Lon[i],Lat[i],AA1[i]))

