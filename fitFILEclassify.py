# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 16:54:15 2021

@author: dingxu
"""

import pandas as pd
import numpy as np
from PyAstronomy.pyasl import foldAt
from PyAstronomy.pyTiming import pyPDM
import matplotlib.pylab as plt
from scipy import interpolate
import os,pickle,time,shutil
import matplotlib.pyplot as plt

def ztf_2(CSV_FILE_PATH,P):
    dfdata = pd.read_csv(CSV_FILE_PATH)
    
    hjd = dfdata['HJD']
    mag = dfdata['mag']
    
    rg = dfdata['band'].value_counts()
    try:
        lenr = rg['r']
    except:
        return [0,0],0
    
    nphjd = np.array(hjd)
    npmag = np.array(mag)
    
    hang = rg['g']
    nphjd = nphjd[hang:]
    npmag = npmag[hang:]-np.mean(npmag[hang:])
    
    phases = foldAt(nphjd, P)
    sortIndi = np.argsort(phases)
    phases = phases[sortIndi]
    resultmag = npmag[sortIndi]
    
    listmag = resultmag.tolist()
    listmag.extend(listmag)
    
    listphrase = phases.tolist()
    listphrase.extend(listphrase+np.max(1)) 
    
    
    nplistmag = np.array(listmag)
    sortmag = np.sort(nplistmag)
    try:
        maxindex = np.median(sortmag[-15:])
    except:
        return [0,0],0
    indexmag = listmag.index(maxindex)
    
    nplistphrase = np.array(listphrase)
    nplistphrase = nplistphrase-nplistphrase[indexmag]
    nplistmag = np.array(listmag)
    
    phasemag = np.vstack((nplistphrase, nplistmag)) #纵向合并矩阵
    phasemag = phasemag.T
    
    phasemag = phasemag[phasemag[:,0]>=0]
    phasemag = phasemag[phasemag[:,0]<=1]
    
    dmag1=np.diff(phasemag[:,1],2).std()/np.sqrt(6)
    
    return phasemag,dmag1



#tot=781602
tot=781602
w=10000
t1w=tot//w
dat=np.genfromtxt('Table2data.txt',dtype=str)

datemp = []
tot=dat.shape[0]
ID=0
for j in range(t1w+1):
    for i in range(w):
        if ID>(tot-1):
            break
        sourceid=dat[ID,1]
        P = float(dat[ID][4])
        gmag=float(dat[ID,8])
        dirnm='Z:/DingXu/ZTF_jkf/alldata/'+str(int(sourceid)//w).zfill(4)
        filename = dirnm+'/'+str(sourceid).zfill(7)+'.csv'
        
        if os.path.getsize(filename)>100:
            if gmag<20:
                print(dat[ID,24].upper())
                print(sourceid)
                
                if (dat[ID,24].upper()=='RR'):
                    pm,d1 = ztf_2(filename, P)
  
                    if d1 <= 0.2:
                        #np.savetxt('H:\\ZTFDATA\\EA\\'+sourceid+'.txt', pm)
                        try:
                            
                            sx1 = np.linspace(0,1,100)
                            sy1 = np.interp(sx1, pm[:,0], pm[:,1])
                        
                            sx1sy1 = np.vstack((sx1, sy1)) #纵向合并矩阵
                            sx1sy1 = sx1sy1.T
                            np.savetxt('H:\\ZTFDATA\\RR\\'+sourceid+'.txt', sx1sy1)
                            
#                            plt.figure(0)
#                            plt.plot(pm[:,0], pm[:,1], '.')
#                            plt.plot(sx1,sy1,'.')
#                            plt.pause(0.001)
#                            plt.clf()
#                            ax1 = plt.gca()
#                            ax1.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
#                            ax1.invert_yaxis() #y轴反向 
                        
                        except:
                            print('it is error!')
                            
                            
                if (dat[ID,24].upper()=='RRC'):
                    pm,d1 = ztf_2(filename, P)
  
                    if d1 <= 0.2:
                        #np.savetxt('H:\\ZTFDATA\\EA\\'+sourceid+'.txt', pm)
                        try:
                            sx1 = np.linspace(0,1,100)
                            sy1 = np.interp(sx1, pm[:,0], pm[:,1])
                            
                            sx1sy1 = np.vstack((sx1, sy1)) #纵向合并矩阵
                            sx1sy1 = sx1sy1.T
                            np.savetxt('H:\\ZTFDATA\\RRC\\'+sourceid+'.txt', sx1sy1)
                            
#                            plt.figure(1)
#                            plt.plot(pm[:,0], pm[:,1], '.')
#                            plt.plot(sx1,sy1,'.')
#                            plt.pause(0.001)
#                            plt.clf()
#                            ax2 = plt.gca()
#                            ax2.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
#                            ax2.invert_yaxis() #y轴反向 
                        
                        except:
                            print('it is error!')
                    
                if (dat[ID,24].upper()=='DSCT'):
                    pm,d1 = ztf_2(filename, P)
  
                    if d1 <= 0.02:
                        #np.savetxt('H:\\ZTFDATA\\EA\\'+sourceid+'.txt', pm)
                        try:
                            sx1 = np.linspace(0,1,100)
                            sy1 = np.interp(sx1, pm[:,0], pm[:,1])
                            
                            sx1sy1 = np.vstack((sx1, sy1)) #纵向合并矩阵
                            sx1sy1 = sx1sy1.T
                            np.savetxt('H:\\ZTFDATA\\DSCT\\'+sourceid+'.txt', sx1sy1)
                            
#                            plt.figure(1)
#                            plt.plot(pm[:,0], pm[:,1], '.')
#                            plt.plot(sx1,sy1,'.')
#                            plt.pause(0.001)
#                            plt.clf()
#                            ax2 = plt.gca()
#                            ax2.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
#                            ax2.invert_yaxis() #y轴反向 
                        
                        except:
                            print('it is error!')
                            
                    
        ID+=1
            