# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 15:56:50 2021

@author: dingxu
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from PyAstronomy.pyasl import foldAt
from PyAstronomy.pyTiming import pyPDM
from astropy.timeseries import LombScargle
from tensorflow.keras.models import load_model


def readfits(fits_file):
    with fits.open(fits_file, mode="readonly") as hdulist:
        tess_bjds = hdulist[1].data['TIME']
        sap_fluxes = hdulist[1].data['SAP_FLUX']
        pdcsap_fluxes = hdulist[1].data['PDCSAP_FLUX']
        print(hdulist[0].header['OBJECT'])
        print(hdulist[0].header['RA_OBJ'], hdulist[0].header['DEC_OBJ'])
        
        indexflux = np.argwhere(pdcsap_fluxes > 0)
#        print(sap_fluxes)
        time = tess_bjds[indexflux]
        time = time.flatten()
        flux = pdcsap_fluxes[indexflux]
        flux =  flux.flatten()
        
        return time, flux
    
    

path = 'I:\\TESSDATA\\section1\\' #tess2018206045859-s0001-0000000025063986-0120-s_lc.fits
file = 'tess2018206045859-s0001-0000000092254719-0120-s_lc.fits'
#file = "https://archive.stsci.edu/missions/tess/tid/s0001/0000/0000/2515/5310/tess2018206045859-s0001-0000000025155310-0120-s_lc.fits"
tbjd, fluxes = readfits(path+file)
mag = -2.5*np.log10(fluxes)
mag = mag-np.mean(mag)

plt.figure(0)
plt.plot(tbjd, fluxes, '.')
plt.xlabel('JD',fontsize=14)
plt.ylabel('FLUX',fontsize=14) 

plt.figure(1)
plt.plot(tbjd, mag, '.')
plt.xlabel('phase',fontsize=14)
plt.ylabel('mag',fontsize=14)
ax1 = plt.gca()
ax1.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
ax1.invert_yaxis() #y轴反向


#magjddata = np.column_stack((tbjd, mag))   
#
#np.savetxt('magjd.txt', magjddata)
