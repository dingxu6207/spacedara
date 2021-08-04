# -*- coding: utf-8 -*-
"""
Created on Wed Jun 30 20:35:44 2021

@author: dingxu
"""

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from PyAstronomy.pyasl import foldAt
from PyAstronomy.pyTiming import pyPDM
from astropy.timeseries import LombScargle
from tensorflow.keras.models import load_model

model = load_model('resultztfmodel.hdf5')
def classfiydata(phasemag):
    sx1 = np.linspace(0,1,100)
    sy1 = np.interp(sx1, phasemag[:,0], phasemag[:,1])
    nparraydata = np.reshape(sy1,(1,100))
    prenpdata = model.predict(nparraydata)

    index = np.argmax(prenpdata[0])
    print(index)
    return index
    
    
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


def computeperiod(JDtime, targetflux):
   
    ls = LombScargle(JDtime, targetflux, normalization='model')
    frequency, power = ls.autopower(minimum_frequency=0.01,maximum_frequency=20)
    #frequency, power = ls.autopower()
    index = np.argmax(power)
    maxpower = np.max(power)
    period = 1/frequency[index]
#    maxpower = np.max(power)
#    sortpower = np.sort(power)
#    index = np.argwhere(power == sortpower[-3])
#    index = index[0][0]
#    period = 1/frequency[index]
    
    wrongP = ls.false_alarm_probability(power.max())
    return period*1, wrongP, maxpower

def computePDM(f0, time, fluxes):
    mag = -2.5*np.log10(fluxes)
    mag = mag-np.mean(mag)
    S = pyPDM.Scanner(minVal=f0-0.01, maxVal=f0+0.01, dVal=0.00001, mode="frequency")
    P = pyPDM.PyPDM(time, mag)
    f2, t2 = P.pdmEquiBin(60, S)
    delta = np.min(t2)
    pdmp = 1/f2[np.argmin(t2)]
    return pdmp, delta 
    

def pholddata(per, times, fluxes):
    mags = -2.5*np.log10(fluxes)
    mags = mags-np.mean(mags)
    
    if per<27:
        lendata =  int((per/25)*len(times))
    else:
        lendata = len(times)
        
    time = times[0:lendata]
    mag = mags[0:lendata]
    phases = foldAt(time, per)
    sortIndi = np.argsort(phases)
    phases = phases[sortIndi]
    resultmag = mag[sortIndi]
    return phases, resultmag

def zerophse(phases, resultmag):
    listmag = resultmag.tolist()
    listmag.extend(listmag)
    listphrase = phases.tolist()
    listphrase.extend(listphrase+np.max(1))
    
    nplistmag = np.array(listmag)
    sortmag = np.sort(nplistmag)
    maxindex = np.median(sortmag[-1:])
    indexmag = listmag.index(maxindex)
    nplistphrase = np.array(listphrase)
    nplistphrase = nplistphrase-nplistphrase[indexmag]
    nplistmag = np.array(listmag)
    
    phasemag = np.vstack((nplistphrase, nplistmag)) #纵向合并矩阵
    phasemag = phasemag.T
    phasemag = phasemag[phasemag[:,0]>=0]
    phasemag = phasemag[phasemag[:,0]<=1]
    
    return phasemag

path = 'I:\\TESSDATA\\variable\\section1\\' #tess2018206045859-s0001-0000000025063986-0120-s_lc.fits
file = 'tess2018206045859-s0001-0000000139699256-0120-s_lc.fits'
#file = "https://archive.stsci.edu/missions/tess/tid/s0001/0000/0000/2515/5310/tess2018206045859-s0001-0000000025155310-0120-s_lc.fits"
tbjd, fluxes = readfits(path+file)
comper, wrongP, maxpower = computeperiod(tbjd, fluxes)
pdmp, delta  = computePDM(1/comper, tbjd, fluxes)
pdmp2, delta2  = computePDM(1/(comper*2), tbjd, fluxes)

delm = np.abs(delta-delta2)
if (delta < delta2):
    phases, resultmag = pholddata(comper, tbjd, fluxes)
else:
    phases, resultmag = pholddata(comper*2, tbjd, fluxes)

phasemag = zerophse(phases, resultmag)
index = classfiydata(phasemag)
#t_fit = np.linspace(0, 1)
#y_fit = ls.model(t_fit, 1/comper)
#y_fit = -2.5*np.log10(y_fit)
#y_fit = y_fit-np.mean(y_fit)

plt.figure(0)
plt.plot(tbjd, fluxes, '.')
plt.xlabel('JD',fontsize=14)
plt.ylabel('FLUX',fontsize=14)


plt.figure(1)
plt.plot(phases, resultmag, '.')
plt.xlabel('phase',fontsize=14)
plt.ylabel('mag',fontsize=14)
ax1 = plt.gca()
ax1.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
ax1.invert_yaxis() #y轴反向

plt.figure(2)
plt.plot(phasemag[:,0], phasemag[:,1], '.')
plt.xlabel('phase',fontsize=14)
plt.ylabel('mag',fontsize=14)
ax1 = plt.gca()
ax1.yaxis.set_ticks_position('left') #将y轴的位置设置在右边
ax1.invert_yaxis() #y轴反向

if index == 0:
    plt.title('Prediction is BYDra')
    
if index == 1:
    plt.title('Prediction is DSCT')

if index == 2:
    plt.title('Prediction is EA')

if index == 3:
    plt.title('Prediction is EW')

if index == 4:
    plt.title('Prediction is RR')
    
if index == 5:
    plt.title('Prediction is SR')
