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
fits_file = "https://archive.stsci.edu/missions/tess/tid/s0001/0000/0000/2515/5310/tess2018206045859-s0001-0000000025155310-0120-s_lc.fits"
#fits_file = 'kplr006852488-2013131215648_llc.fits'

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
flux = flux.flatten()
#mag1 = 25 -2.5*np.log10(sap_fluxes)


mag2 = 25-2.5*np.log10(flux)

#frequency, power = LombScargle(time.flatten(), flux.flatten()).autopower()
#lenfre = len(frequency)
#index = np.argmax(power[0:int(lenfre/3)])
#plt.figure(2)
#plt.plot(frequency, power)  
#print(1/frequency[index])

##mag2 = mag2-np.mean(mag2)
plt.figure(1)
##plt.plot(tess_bjds, mag1)
plt.plot(time, mag2)
ax1 = plt.gca()
ax1.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
ax1.invert_yaxis() #y轴反向
#
#
S = pyPDM.Scanner(minVal=0.06, maxVal=10, dVal=0.0001, mode="frequency")
P = pyPDM.PyPDM(time, mag2)
#
#f1, t1 = P.pdmEquiBinCover(10, 3, S)
f2, t2 = P.pdmEquiBin(10, S)
plt.figure(2)
plt.plot(f2, t2, 'gp-')
##plt.plot(f1, t1, 'rp-')
plt.xlabel('frequency',fontsize=14)
plt.ylabel('Theta', fontsize=14)
print(f2[np.argmin(t2)])
#
#valuet = np.sort(t2[t2<0.5])
#print(valuet)

#
P = 3.288776 

phases = foldAt(time, P)
#phases = phases[phases<1]
sortIndi = np.argsort(phases)
phases = phases[sortIndi]
resultmag = mag2[sortIndi]

plt.figure(10)
plt.plot(phases, resultmag)
ax1 = plt.gca()
ax1.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
ax1.invert_yaxis() #y轴反向
