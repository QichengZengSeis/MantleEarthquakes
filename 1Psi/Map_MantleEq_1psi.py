#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 11:23:15 2023

@author: u1318104
"""

#%% Modules and settings

import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from matplotlib.offsetbox import AnchoredText
from cartopy.io.img_tiles import GoogleTiles, Stamen, OSM, MapQuestOSM, MapQuestOpenAerial#,StamenTerrain #StamenTerrain
import matplotlib.cm as cm
import numpy as np
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from obspy.clients.fdsn import Client
from obspy import UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth
import matplotlib
from tqdm import tqdm
import cartopy.feature as cfeature
import os
# from mpl_toolkits.basemap import Basemap
#US
xmax=297
ymax=50 #36.2
xmin=235 #-121.0
ymin=25 #32.2
#WUS
# xmax=260
# ymax=51 #36.2
# xmin=230 #-121.0
# ymin=25 #32.2
fdir=os.path.dirname(__file__)
os.chdir(fdir)
#%% Colormap
def colormap():

    c1=((125./255.),0.,0.)
    c2=((255./255.),0.,0.) #red

    c3=((255./255.),(255./255.),0.) #yellow
    # c4=((230./255.),(230./255.),(230./255.)) #white'
    c4=((255./255.),(255./255.),(255./255.)) #white'

    c5=((150./255.),(255./255.),(150./255.)) #lightgreen
    c6=((100./255.),(185./255.),(255./255.)) #cyan
    c7=((100./255.),0.,(170./255.)) #purple
    
    # colors = [c1,c2,c3,c4,c5,c6,c7]  # black, red, yellow, white, lightgreen, cyan, purple
    colors = [c1,c2,c3,c4,c5,c6,c7]  # black, red, yellow, white, lightgreen, cyan, purple
    
    # colors1=cmap2rgb('jet_r') #new color black without white color in the middle, TBD
    #n_bins = [3, 6, 10, 100]  # Discretizes the interpolation into bins
    cmap_name = 'my_list'
    cmnew = matplotlib.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=500)
    return cmnew
cmnew=colormap()

def colormap1():
    #load colormap
    c0=((255./255.),(255./255.),(255./255.))
    c1=((255./255.),(255./255.),(255./255.))
    c2=((220./255.),(220./255.),(220./255.))

    c3=((0./255.),(255./255.),(255./255.))
    c4=((0./255.),(128./255.),(255./255.))
    c5=((76./255.),(0./255.),(153./255.))
    c6=((76./255.),(0./255.),(153./255.))

    colors = [c0,c1,c2,c3,c4,c5,c6]  # black, red, yellow, white, lightgreen, cyan, purple

    cmap_name = 'my_list'
    cmnew = matplotlib.colors.LinearSegmentedColormap.from_list(cmap_name, colors, N=500)
    return cmnew
cmnew1=colormap1()

#%% Earthquake Catalog
client = Client("IRIS")
startT = UTCDateTime("1800-01-01T00:00:00")
# endT = UTCDateTime("2020-03-18T23:59:59")
endT = UTCDateTime("2023-12-31T23:59:59")

catalog=client.get_events(starttime=startT, endtime=endT, 
                          minlatitude=ymin, maxlatitude=ymax, 
                          minlongitude=xmin if xmin<180 else xmin-360, 
                          maxlongitude=xmin if xmin<180 else xmax-360,
                          mindepth=45, minmagnitude=2)



cat1=[]
for i in range(len(catalog)):
    tmp=catalog[i]
    if any(x in tmp.event_descriptions[0].text for x in ['GULF','OFF','MEXICO','CANADA','VANCOUVER']):
        continue
    
    cat1.append([tmp.origins[0].longitude,tmp.origins[0].latitude,tmp.magnitudes[0].mag])
cat1=np.array(cat1)
tmp_sort=np.argsort(cat1[:,2])[::-1]
cat1=cat1[tmp_sort]

tmpM1=cat1[0,2];tmpM0=cat1[-1,2] #biggest/smallest magnitude
tmpS1=20;tmpS0=2 #biggest/smallest Marker Size

catMS=tmpS0+(tmpS1-tmpS0)*(cat1[:,2]-tmpM0)/(tmpM1-tmpM0)
#%% Geophysical Boundaries

# fdir='/uufs/chpc.utah.edu/common/home/koper-group1/keith/HISTORICAL_PROJECTS/2020.Magna/util/'
# fp=open(fdir+'wasatch_pang.xy','r')
# Faults=fp.read().split('>')

fp=open(fdir+'/US_boundaries/wus_province_II.dat','r')
Physio=fp.read().split('99999 99999')

flag=1 #flag for lon
faultslon=list()
faultslat=list()
for fault in Physio:
    if fault=='':
        continue
    faultlonlat=fault.split()
    tmplon=list()
    tmplat=list()
    for tmp in faultlonlat:
        if flag:
            tmplon.append(float(tmp))
        else:
            tmplat.append(float(tmp))
        flag=1-flag
    faultslon.append(tmplon)
    faultslat.append(tmplat)


fp_rift=open(fdir+'/US_boundaries/rift.dat')
rift_lon=[];rift_lat=[]
flag=0 #0 - beginning of file; 1 - not beginning
for line in fp_rift:
    if '>>' in line:
        if flag:
            rift_lon.append(tmp_lon) #ignore warning sign
            rift_lat.append(tmp_lat) #ignore warning sign
        flag=1
        tmp_lon=[]
        tmp_lat=[]
        continue
    tmp_lon.append(float(line.split()[0]))
    tmp_lat.append(float(line.split()[1]))
rift_lon.append(tmp_lon)
rift_lat.append(tmp_lat)

fp_BD=open(fdir+'/US_boundaries/BD.20.txt')
BD_lon=[];BD_lat=[]
flag=0 #0 - beginning of file; 1 - not beginning
for line in fp_BD:
    if '>>' in line:
        if flag:
            BD_lon.append(tmp_lon) #ignore warning sign
            BD_lat.append(tmp_lat) #ignore warning sign
        flag=1
        tmp_lon=[]
        tmp_lat=[]
        continue
    tmp_lon.append(float(line.split()[0]))
    tmp_lat.append(float(line.split()[1]))
BD_lon.append(tmp_lon)
BD_lat.append(tmp_lat)

fp_GF=open(fdir+'/US_boundaries/Grenville.txt')
GF_lon=[];GF_lat=[]
flag=0 #0 - beginning of file; 1 - not beginning
for line in fp_GF:
    if '>>' in line:
        if flag:
            GF_lon.append(tmp_lon) #ignore warning sign
            GF_lat.append(tmp_lat) #ignore warning sign
        flag=1
        tmp_lon=[]
        tmp_lat=[]
        continue
    tmp_lon.append(float(line.split()[0]))
    tmp_lat.append(float(line.split()[1]))
GF_lon.append(tmp_lon)
GF_lat.append(tmp_lat)

# Physio_Prov_US.txt
fp_PP=open(fdir+'/US_boundaries/Physio_Prov_US.txt')
PP_lon=[];PP_lat=[]
tmp_lon=[]; tmp_lat=[]
flag=0 #0 - beginning of file; 1 - not beginning
for line in fp_PP:
    if '360 NaN' in line:
        if flag:
            PP_lon.append(tmp_lon) #ignore warning sign
            PP_lat.append(tmp_lat) #ignore warning sign
        flag=1
        tmp_lon=[]
        tmp_lat=[]
        continue
    tmp_lon.append(float(line.split()[0]))
    tmp_lat.append(float(line.split()[1]))
PP_lon.append(tmp_lon)
PP_lat.append(tmp_lat)


fp_MCR=open(fdir+'/US_boundaries/1.2_1.1.MCR.gmt')
MCR_lon=[];MCR_lat=[];
flag=0 #0 - beginning of file; 1 - begin of reading for one section; 2 - end of reading for one section
for line in fp_MCR:
    if '# @P' in line:
        flag=1
        tmp_lon=[]
        tmp_lat=[]
        continue
    if '>' in line and flag==1:
        flag=2
        MCR_lon.append(tmp_lon) #ignore warning sign
        MCR_lat.append(tmp_lat) #ignore warning sign
        continue
    if flag==1:
        tmp_lon.append(float(line.split()[0]))
        tmp_lat.append(float(line.split()[1]))
MCR_lon.append(tmp_lon)
MCR_lat.append(tmp_lat)



#%% 1Psi
# dir_EMC='/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi/'
# LatsEMC, LonsEMC, DepsEMC, VsEMC=np.loadtxt(dir_EMC+'US_CrustVs_SLK_GRL_2015.txt',unpack=True)

# dir_Vph='/uufs/chpc.utah.edu/common/home/u1318104/Research/NearFieldXcorr/TOMO_BARMIN_ET_AL/30.0_iso/400_70/'
# Lons,Lats,Vphs=np.loadtxt(dir_Vph+'USArr-iso_400_70_30.0.1',unpack=True)

dir_Vph=fdir+'/Data_1psi/'
# Lons,Lats,Vphs=np.loadtxt(dir_Vph+'60s_iso_ani_v1_scale_1.iso',unpack=True,usecols=(0,1,2))
Lons,Lats,Vphs,A1,Psi1=np.loadtxt(dir_Vph+'USArr_1psi_0.txt',unpack=True,usecols=(0,1,2,3,4))

SubN=3 #Subsample ratio of 1psi arrows
Cri1=A1>0.02#np.where(A1>0.02)#0.0001
Cri2=~np.array(np.round((Lons-np.min(Lons))/0.2)%SubN,dtype=bool)
#US
Cri3=~np.array(np.round((Lats-np.min(Lats))/0.2)%SubN,dtype=bool)
#AK
# Cri3=~np.array(np.round((Lats-np.min(Lats))/0.1)%SubN,dtype=bool)

CriS=np.logical_and.reduce((Cri1,Cri2,Cri3))
# CriS=Cri1

LonsC=Lons[CriS]
LatsC=Lats[CriS]
A1C=A1[CriS]
Psi1C=Psi1[CriS]
#%% Mesh for 2D plot
# from PlotAll import tightgaussiansmooth
import matplotlib 

# Lats60km=LatsEMC[DepsEMC==60]
# Lons60km=LonsEMC[DepsEMC==60]
# Vs60km=VsEMC[DepsEMC==60]

Lon1=np.unique(Lons);Lon1.sort()
Lat1=np.unique(Lats);Lat1.sort()
DATA=np.zeros((len(Lat1),len(Lon1)))
# if scale==[]:
#     scale=[min(data),max(data)]
    
for k in range(len(Vphs)):
    iLon=np.where(Lon1==Lons[k])
    iLat=np.where(Lat1==Lats[k])
    # DATA[iLat,iLon]=Vphs[k]
    DATA[iLat,iLon]=A1[k]
    
DATA=np.ma.array(DATA,mask=~np.array(DATA,dtype=bool))

#%% Plotting setting

ftname=[f.name for f in matplotlib.font_manager.fontManager.afmlist]
import os
from matplotlib import font_manager
import sys
# sys.path.append()
# sys.path.insert(0,'/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi')
# sys.path=list(np.unique(np.array(sys.path)))
# sys.path
import test_obs_1psi as to1


plt.rcParams.update({
    "text.usetex": True,
    "font.family": 'Times New Roman',
    "font.size": 11,
    "figure.autolayout": True}) #'axes.linewidth': 0.8
#%%

#WUS
# lonticks=np.arange(xmin,xmax,5)
# latticks=np.arange(ymin,ymax,5)
#AK
#US
# lonticks=np.arange(xmin,xmax,10)
# latticks=np.arange(ymin,ymax,5)
lonticks=np.arange(-120,-70+1,10)
latticks=np.arange(25,50+1,5)

# request=MapQuestOpenAerial() 
# request=MapQuestOSM()
# request=StamenTerrain()
request=GoogleTiles(style='satellite')
# GoogleTiles(style='satellite') #StamenTerrain()
request.desired_tile_form='L'

MarkerSize=12.5 #50
# fig=plt.figure(figsize=(12, 8)) #, facecolor="none"


ProjectCCRS=ccrs.PlateCarree()#ccrs.LambertConformal()
fig, ax0=plt.subplots(figsize=(8.5, 6),subplot_kw=dict(projection=request.crs)) #,subplot_kw=dict(projection=request.crs)#12 8
ax0.set_facecolor('gainsboro')
## Reset Zoom level ###
# xmine=-112.5; xmaxe=-111.6; ymine=40.4; ymaxe=40.9
xmine=xmin; xmaxe=xmax; ymine=ymin; ymaxe=ymax

ax0.cla()
ax0.set_extent([xmine,xmaxe,ymine,ymaxe], crs=ProjectCCRS)
# ax0.add_image(request,4,cmap='terrain') #,cmap='gray'
# ax0.add_image(request,4,cmap='gray') #,cmap='gray'

# Tick labels etc
ax0.set_xticks(lonticks, crs=ProjectCCRS) #, crs=ProjectCCRS
ax0.set_yticks(latticks, crs=ProjectCCRS) #, crs=ProjectCCRS
lon_formatter = LongitudeFormatter(number_format='.1f', degree_symbol='')  # ,zero_direction_label=True)
lat_formatter = LatitudeFormatter(number_format='.1f', degree_symbol='')
ax0.xaxis.set_major_formatter(lon_formatter)
ax0.yaxis.set_major_formatter(lat_formatter)
ax0.xlabel_style = {'size': 11, 'color': 'gray'}
ax0.ylabel_style = {'size': 11, 'color': 'gray'}  # 'weight':'bold'

#plot faults...
for ii in np.arange(len(faultslon)):
    ax0.plot(faultslon[ii],faultslat[ii],linewidth=2,color='r',zorder=2,transform=ProjectCCRS) #transform=ccrs.Geodetic() #linewidth=2.
for ii in np.arange(len(rift_lon)):
    if ii !=0:
        continue
    ax0.plot(rift_lon[ii],rift_lat[ii],linewidth=2,color='r',zorder=2,transform=ProjectCCRS) #transform=ccrs.Geodetic() #linewidth=2.
for ii in np.arange(len(BD_lon)):
    ax0.plot(BD_lon[ii],BD_lat[ii],linewidth=2,color='r',zorder=2,transform=ProjectCCRS) #transform=ccrs.Geodetic() #linewidth=2.
# for ii in np.arange(len(GF_lon)):
#     ax0.plot(GF_lon[ii],GF_lat[ii],linewidth=2,color='r',zorder=2,transform=ccrs.Geodetic()) #transform=ccrs.Geodetic() #linewidth=2.

# for ii in np.arange(len(PP_lon)):
#     ax0.plot(PP_lon[ii],PP_lat[ii],linewidth=2,color='r',zorder=2,transform=ccrs.Geodetic()) #transform=ccrs.Geodetic() #linewidth=2.


im=ax0.contourf(Lon1,Lat1,DATA*100,np.linspace(0,6,7),cmap=cmnew1,extend='both',zorder=1,alpha=1,transform=ccrs.PlateCarree())    
fig.colorbar(im,orientation='horizontal',ticks=np.linspace(0,6,7),pad=0.06,fraction=0.046,shrink=0.5,format='%.0f',label='1psi (\%)')

# ax0.add_feature(cfeature.STATES.with_scale('1m'),linewidth=0.5) #.with_scale('110m')
ax0.add_feature(cfeature.LAKES)
ax0.add_feature(cfeature.BORDERS)
ax0.add_feature(cfeature.STATES.with_scale('50m'),linewidth=0.3)
ax0.add_feature(cfeature.LAND)
ax0.add_feature(cfeature.OCEAN)
###############events>45km###############
for icat in range(len(cat1)):
    # ax0.plot(cat1[icat,0],cat1[icat,1],'.',MarkerSize=catMS[icat],color='gray',transform=ccrs.PlateCarree())
    ax0.scatter(cat1[icat,0],cat1[icat,1],s=catMS[icat]*2,color='gray',edgecolor='k',transform=ccrs.PlateCarree())
lonSRPtmp0=-112.826942
latSRPtmp0=43.384350
lonSRPtmp1=-114.461650
latSRPtmp1=45.427743

lonSRP=(lonSRPtmp0+lonSRPtmp1)/2
latSRP=(latSRPtmp0+latSRPtmp1)/2
# lonSRP=(lonSRPtmp0*2/3+lonSRPtmp1*1/3)
# latSRP=(latSRPtmp0*2/3+latSRPtmp1*1/3)

# ax0.scatter(lonSRPtmp0,latSRPtmp0,s=50,color='k',marker='^',zorder=3,transform=ProjectCCRS)
# ax0.scatter(lonSRPtmp1,latSRPtmp1,s=50,color='b',marker='^',zorder=3,transform=ProjectCCRS)
# ax0.scatter(lonSRP,latSRP,s=50,color='g',marker='o',zorder=3,transform=ProjectCCRS)

##################################### Inset ####################
left,bottom, width, height=[0.735,0.22,0.22,0.22]

ax2=fig.add_axes([left,bottom,width,height])

ax2.plot(to1.aa1,to1.C_a1,color='g',label='1-psi %.4f \n2-psi %.4f'%(to1.A1,to1.A2),zorder=0)
ax2.errorbar(to1.azi1,to1.vel1,yerr=to1.un_vel1,fmt='.',color='r',capsize=2,elinewidth=1,ms=4)
ax2.text(0,3.89,'1-psi 3.0\%\n 2-psi 0.1\%')

# ax2.legend()
ax2.set_xlabel('Azimuth (degree)')
ax2.set_ylabel('Phase Velocity (km/s)')

ax2.yaxis.tick_right()
ax2.yaxis.set_label_position('right')
ax2.tick_params(axis='y',direction='in',pad=-17)
ax2.tick_params(axis='x',direction='in',pad=-12)

# ax2.set_xlim([-100,370])
ax2.set_xticks([0,100,200,300])
ax2.set_yticks([3.8,3.9])
ax2.set_ylim([3.7,3.95])
# ax2.grid()

TRANSPARENT=False
plt.savefig(fdir+'/US_1psi_Eqs.png',transparent=TRANSPARENT,dpi=300)
plt.savefig(fdir+'/US_1psi_Eqs.pdf',transparent=TRANSPARENT)
