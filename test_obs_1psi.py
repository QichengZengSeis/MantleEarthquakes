#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 11:09:53 2023

@author: u1318104
"""

#%%
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from numpy import cos, sin, pi
from numpy.linalg import inv
import os 

dir_data=os.path.dirname(__file__) #'/uufs/chpc.utah.edu/common/home/u1318104/Research/1Psi/' #60s_AK/
fp=dir_data+'/Data_1psi/60s_iso_ani_v1_scale_0.ani'

#%%
Lon=[]
Lat=[]
Azi=[]
Vel=[]
Un_vel=[]

with open(fp) as f:
    # while True:
    lines=f.readlines()
    
    flag=1 # title
    # for line in tqdm(lines):
    for line in lines:
        if flag:
            n=int(line.split()[2])
            Lon.append(float(line.split()[0]))
            Lat.append(float(line.split()[1]))
            
            azi=[]
            vel=[]
            un_vel=[]
            tmpn=0
            flag=0
            continue
        
        if not flag:
            azi.append(float(line.split()[0]))
            vel.append(float(line.split()[1]))
            un_vel.append(float(line.split()[2]))
            tmpn+=1
            
        if tmpn==n:
            Azi.append(azi)
            Vel.append(vel)
            Un_vel.append(un_vel)
            flag=1
            
#%%
CCi=[]
AA1=[]
AA2=[]
Psi1=[]
Psi2=[]
# for ii in tqdm(range(len(Azi))):
for ii in range(len(Azi)):
    azi=np.asarray(Azi[ii])
    vel=np.asarray(Vel[ii])
    # vel=4*(1+2/2*cos(azi/180*pi-pi/2)+0.1/2*cos(2*(azi/180*pi-pi/6)))
    
    un_vel=np.asarray(Un_vel[ii])
    
    #Am=d, W weighting; WAm=Wd;m: Ci,A1cos(*),A1sin(*),A2cos(2*),A2sin(2*)
    A=np.ones((len(azi),5))
    d=vel.reshape((-1,1))
    A[:,1]=cos(azi/180*pi)/2
    A[:,2]=sin(azi/180*pi)/2
    A[:,3]=cos(2*azi/180*pi)/2
    A[:,4]=sin(2*azi/180*pi)/2
    
    W=np.diag(1/un_vel)
    # W=np.diag(np.ones(un_vel.shape))
    m=inv(A.T @ W @ W @ A) @ (A.T @ W @ W @ d)
    # m=np.matmul(inv(np.matmul(A.T,A)),np.matmul(A.T,d))
    
    Ci=m[0]
    A1=np.sqrt(m[1]**2+m[2]**2)/Ci
    psi1=np.arctan2(m[2],m[1])
    A2=np.sqrt(m[3]**2+m[4]**2)/Ci
    psi2=np.arctan2(m[4],m[3])/2
    
    CCi.append(Ci);AA1.append(A1);AA2.append(A2)
    Psi1.append(psi1);Psi2.append(psi2)
    
    if len(azi)<=5:
        print(len(azi),Lon[ii],Lat[ii])
    # if ii==3308:
    #     break

#%%
lonSRPtmp0=-112.826942
latSRPtmp0=43.384350
lonSRPtmp1=-114.461650
latSRPtmp1=45.427743

lonSRP=(lonSRPtmp0+lonSRPtmp1)/2
latSRP=(latSRPtmp0+latSRPtmp1)/2
# lonSRP=(lonSRPtmp0*2/3+lonSRPtmp1*1/3)
# latSRP=(latSRPtmp0*2/3+latSRPtmp1*1/3)
# lonSRP=(lonSRPtmp0*3/5+lonSRPtmp1*2/5)
# latSRP=(latSRPtmp0*3/5+latSRPtmp1*2/5)

# lonSRP=-113.451797 #TA.I14A
# latSRP=43.9286

tmp=np.argmin((np.array(Lon)-360-lonSRP)**2+(np.array(Lat)-latSRP)**2)
#%%
CCi=np.array(CCi).reshape(-1)
AA1=np.array(AA1).reshape(-1)
AA2=np.array(AA2).reshape(-1)
Psi1=np.array(Psi1).reshape(-1)
Psi2=np.array(Psi2).reshape(-1)
fp=[]
# for i in tqdm(range(len(CCi))):
for i in range(len(CCi)):
    fp.append('%f %f %f %f %f %f %f'%(Lon[i],Lat[i],CCi[i],AA1[i],Psi1[i],AA2[i],Psi2[i]))
    
# np.savetxt(dir_data+'USArr_1psi_0.txt',fp,fmt='%s')

#%%
ii=tmp#3399#3308
A1=AA1[ii]
A2=AA2[ii]
psi1=Psi1[ii]
psi2=Psi2[ii]
Ci=CCi[ii]
azi=Azi[ii]
vel=Vel[ii]
un_vel=Un_vel[ii]
aa=np.linspace(0,360,40)
C_a=Ci*(1+A1/2*cos(aa/180*pi-psi1)+A2/2*cos(2*(aa/180*pi-psi2)))

aa1=np.array(aa)+180
azi1=np.array(azi)+180
for i in range(len(aa1)):
    if aa1[i]>360:
        aa1[i]=aa1[i]-360
for i in range(len(azi1)):
    if azi1[i]>360:
        azi1[i]=azi1[i]-360
C_a1=C_a[np.argsort(aa1)]
aa1=np.sort(aa1)
vel1=np.array(vel)[np.argsort(azi1)]
un_vel1=np.array(un_vel)[np.argsort(azi1)]
azi1=np.sort(azi1)
#%%
if __name__=="__main__":
    # plt.close('all')
    plt.figure()
    plt.plot(aa,C_a,label='lon %.2f lat %.2f\n1-psi %.4f \n2-psi %.4f'%(Lon[ii],Lat[ii],A1,A2))
    plt.errorbar(azi,vel,yerr=un_vel,fmt='.',color='r',capsize=4)
    plt.legend()
    plt.xlabel('back-azimuth (degs)')
    plt.ylabel('Vph (km/s)')
    # plt.savefig('/uufs/chpc.utah.edu/common/home/u1318104/Figures/04272023_USArr_1PsiV3/SRP_1psi',transparent=True)
    # plt.ylim([3.65,4.0])
    
        
        
    
    
    
    
    
    
    
    
    
