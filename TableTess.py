# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 17:09:46 2021

@author: dingxu
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from PyAstronomy.pyasl import foldAt
from PyAstronomy.pyTiming import pyPDM
from astropy.timeseries import LombScargle

# For the purposes of this tutorial, we just know the MAST URL location of the file we want to examine.
#fits_file = "https://archive.stsci.edu/missions/tess/tid/s0001/0000/0000/2515/5310/tess2018206045859-s0001-0000000025155310-0120-s_lc.fits"
fits_file = 'kplr006852488-2013131215648_llc.fits'

with fits.open(fits_file, mode="readonly") as hdulist:
    tess_bjds = hdulist[1].data['TIME']
    sap_fluxes = hdulist[1].data['SAP_FLUX']
    pdcsap_fluxes = hdulist[1].data['PDCSAP_FLUX']
    print(len(tess_bjds))
    
print(fits.info(fits_file))   

print(hdulist[0].header['RA_OBJ'], hdulist[0].header['DEC_OBJ'])
plt.figure(0)
plt.plot(tess_bjds, sap_fluxes)
#plt.plot(tess_bjds, pdcsap_fluxes)

#pdcsap_fluxes = pdcsap_fluxes[0:500]
#tess_bjds = tess_bjds[0:500]

indexflux = np.argwhere(pdcsap_fluxes > 0)
time = tess_bjds[indexflux]
time = time.flatten()
flux = pdcsap_fluxes[indexflux]
flux =  flux.flatten()
#mag1 = 25 -2.5*np.log10(sap_fluxes)


mag2 = 25-2.5*np.log10(flux)
mag2 = mag2-np.mean(mag2)

ls = LombScargle(time, mag2, normalization='model')
frequency, power = ls.autopower(minimum_frequency=0.01,maximum_frequency=20)
lenfre = len(frequency)
index = np.argmax(power)
plt.figure(2)
plt.plot(frequency, power)  
print('day=', 1/frequency[index])
print('probility=', ls.false_alarm_probability(power.max()) ) 
##mag2 = mag2-np.mean(mag2)
plt.figure(1)
##plt.plot(tess_bjds, mag1)
plt.plot(time, mag2)
ax1 = plt.gca()
ax1.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
ax1.invert_yaxis() #y轴反向

P = 1/frequency[index]
lendata = 2*int((P/27)*len(time))
time = time[0:lendata]
mag2 = mag2[0:lendata]
phases = foldAt(time, P)
#phases = phases[phases<1]
sortIndi = np.argsort(phases)
phases = phases[sortIndi]
resultmag = mag2[sortIndi]

plt.figure(10)
plt.plot(phases, resultmag, '.')
ax1 = plt.gca()
ax1.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
ax1.invert_yaxis() #y轴反向
